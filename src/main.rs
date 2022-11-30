#![allow(dead_code, unused_variables, unused_imports, unused_mut)]
#![feature(drain_filter)]
//TODO: remove this once I have some working code

mod particles;
mod util;
use atomic_float::AtomicF64;

extern crate nalgebra as na;
use crate::particles::Particles;
use crate::util::*;
use clap::{Arg, Parser};
use crossbeam::channel::{self, unbounded, Receiver, Sender};
use na::coordinates::X;
use num::{Bounded, Integer};
use particles::{Particle, ParticlePointer, ParticleType};
use rand::prelude::*;
use rand::seq::index;
use rayon::prelude::*;
use std::cell;
use std::cmp::max;
use std::f64::consts::PI;
use std::fmt::Display;
use std::fs::File;
use std::io::Write;
use std::ops::ControlFlow;
use std::path::Path;
use std::sync::{Arc, RwLock};

//const CHONKY: usize = 8;
const COLLISION_SIZE: f32 = 1e-28;
//number of molecules each particle represents
const molecules_per_particle: f64 = 1e27;
const sigma_k: f64 = 1.000000003e-28;

fn random_direction() -> na::Vector3<f64> {
    let b = 2. * rand::random::<f64>() - 1.;
    let a: f64 = (1. - b * b).sqrt();
    let theta = rand::random::<f64>() * 2. * PI;
    na::Vector3::<f64>::new(b, a * theta.cos(), a * theta.sin())
}
/// Computes a velocity magnitude distribution that would correspond to
/// thermal equilibrium
fn random_velocity(region_temp: f64) -> f64 {
    return region_temp as f64 * (-f64::max(rand::random::<f64>(), 1e-200).ln()).sqrt();
}

// Information that is used to control the collision probability code
struct CollisionData {
    // Maximum collision rate seen for this cell so far in the simulation
    max_collision_rate: Vec<f64>,
    // Non-integral fraction of collisions that remain to be performed
    // and are carried over into the next timestep
    remainder: Vec<f64>,
    //maxreciver:
}

impl CollisionData {
    fn new(num_cells: usize, region_temp: f64) -> Self {
        Self {
            max_collision_rate: vec![sigma_k * region_temp; num_cells],
            remainder: vec![0.; num_cells],
        }
    }
    //get the collision rate and remainder for a cell
    fn get_data(&self, index: usize) -> (f64, f64) {
        let collision_rate = self.max_collision_rate[index];

        let remainder = self.remainder[index];
        (collision_rate, remainder)
    }

    fn set_remainder(&mut self, index: usize, remainder: f64) {
        //let mut remainders = self.max_collision_rate[index];
        self.remainder[index] = remainder;
    }
}

struct CellSample {
    mean_particles: f64,
    total_velocity: na::Vector3<f64>,
    total_kinetic_energy: f64,
    members: Vec<usize>,
}

impl CellSample {
    fn new() -> Self {
        Self {
            mean_particles: 0.,
            total_velocity: na::Vector3::zeros(),
            total_kinetic_energy: 0.0,
            //Note: balance of reads to writes is about even, so use mutex
            //https://stackoverflow.com/questions/50704279/when-or-why-should-i-use-a-mutex-over-an-rwlock
            members: Vec::new(),
        }
    }
    fn reset(&mut self) {
        self.mean_particles = 0.;
        self.total_kinetic_energy = 0.;
        //https://stackoverflow.com/a/56114369/11019565
        self.total_velocity.x = 0.;
        self.total_velocity.y = 0.;
        self.total_velocity.z = 0.;
    }
}

// Initialize particles at inflow boundaries
/// Create particles at inflow boundary
/// This works by creating particles at a ghost cell just before the boundary
/// any particles that don't make it into the domain are discarded.

fn initialize_boundaries(
    particles: &mut Vec<Particle>,
    num_x: usize,
    num_y: usize,
    num_z: usize,
    mean_velocity: f64,
    region_temp: f64,
    mean_particles_per_cell: i32,
) {
    // TODO: I might be able to avoid resizing the particle containing
    // data structure by running the functions
    // move_particles and remove_outside_particles (or similar function)
    //directly on the generated particles
    let dx = 2. / num_x as f64;
    let dy = 2. / num_y as f64;
    let dz = 2. / num_z as f64;
    particles.extend(
        &(0..num_y)
            .into_par_iter()
            .flat_map(|j| -> Vec<Particle> {
                (0..num_z)
                    .into_par_iter()
                    .flat_map(|k| -> Vec<Particle> {
                        let current_x = -1. - dx;
                        let current_y = -1. + j as f64 * dy;
                        let current_z = -1. + k as f64 * dz;

                        (0..mean_particles_per_cell)
                            .into_par_iter()
                            //.chunks(CHONKY)
                            .map(|i| -> Particle {
                                let position = na::Vector3::<f64>::new(
                                    current_x + rand::random::<f64>() * dx,
                                    current_y + rand::random::<f64>() * dy,
                                    current_z + rand::random::<f64>() * dz,
                                );
                                let velocity = random_velocity(region_temp) * random_direction();
                                Particle::new(position, velocity)
                            })
                            .collect::<Vec<Particle>>()
                    })
                    .collect()
            })
            .collect::<Vec<Particle>>(),
    );
}

/// Move particle for the timestep.  
///Also handle side periodic boundary conditions
/// and specular reflections off of a plate
fn move_particles_with_bcs(
    particles: &mut Vec<Particle>,
    plate_x: f64,
    plate_dy: f64,
    plate_dz: f64,
    delta_t: f64,
) {
    particles.into_iter().for_each(|particle| {
        let tmp = particle.position;
        let mut new_position = tmp + particle.velocity * delta_t;

        if (tmp[0] < plate_x) && (new_position[0] > plate_x)
            || (tmp[0] > plate_x) && (new_position[0] < plate_x)
        {
            let t = (tmp[0] - plate_x) / (tmp[0] - new_position[0]);
            //create a particle representing the potential position
            let tmp = tmp * (1. - t) + ((tmp + delta_t * particle.velocity) * t);
            if (-plate_dy..plate_dy).contains(&tmp[1]) && (-plate_dz..plate_dz).contains(&tmp[2]) {
                new_position.x -= 2. * new_position.x - plate_x;
                particle.velocity[0] = -particle.velocity[0];
                particle.particle_type = ParticleType::Thud;
            }
        }

        if (-1.0..1.).contains(&new_position.y) {
            new_position.y -= 2.0;
        }
        if (-1.0..1.).contains(&new_position.z) {
            new_position.z -= 2.0;
        }
        particle.position = new_position;
    });
}

/// any particles outside of the cells need to be discarded as they cannot be indexed.  
/// Since this can happen only at the x=-1 or x=1 boundaries, we only need to check the x coordinate
fn remove_outside_particles(particles: &mut Vec<Particle>, num_cells: usize) {
    // particles.sort_by_key(|p| (p.parent_cell));
    // if let Some(cutoff) = particles.iter().rposition(|p| p.parent_cell < num_cells) {
    //     particles.truncate(cutoff + 1);
    // }
    let _ = particles.drain_filter(|p| p.parent_cell < num_cells);
}
// Move particles based on their present positions and velocities and the
//time-step
//$x^{n+1}_{p} = x^{n}_{p} + v_{p}^n ∗ ∆t$
// For particles that collide with solid objects, apply specular boundary condition
// Particles that move outside of the mesh are removed
// Locate particles to the containing cell in the mesh
//indexParticles()
// Collect statistics of particles contained in each cell (e.g. number of particles and average velocity and kinetic energy)
//sampleParticles()
// Visit cells and determine chance that particles can collide
//fn collide_particles() {

//}

//collideParticles() : collision is a random process where the result of
//a collision conserves momentum and kinetic energy
// time-step is complete, continues
fn main() {
    //-----------------------------INITIALIZATION-----------------------------
    //------------------------------------------------------------------------
    //
    let num_cores = std::thread::available_parallelism().unwrap();
    let mean_velocity: f64 = 1.; //switched to f32 to avoid repetitive casting
                                 //let (s, r) = channel::unbounded();
                                 //struct destructuring for the win
    let Args {
        two_dim,
        num_x,
        num_y,
        num_z,
        mean_particles_per_cell,
        mach_number,
        density,
        plate_x,
        plate_height,
        plate_width,
        time_step,
    } = Args::parse();

    //temperature of the region, affects velocity #TODO: update description
    let region_temp: f64 = mean_velocity / mach_number;

    let number_of_cells = num_x * num_y * num_z;
    // Compute number of molecules a particle represents
    let cell_vol: f64 = 2. / (number_of_cells) as f64;

    //create channels to share between threads
    let (collision_idx_sender, collision_idx_receiver) = unbounded();
    let (sample_data_sender, sample_data_reciever) = unbounded();

    // Create simulation data structures
    //let mut particles: Particles = Particles::new(collision_idx_receiver, sample_data_sender);
    let mut particles: Vec<Particle> = Vec::new();
    let mut cell_data: &mut Vec<CellSample> = &mut Vec::with_capacity(number_of_cells);
    println!("{:?}", number_of_cells);
    //go ahead and initialize the cell_data
    for i in 0..number_of_cells {
        cell_data.push(CellSample::new());
    }

    let mut collision_data: &mut CollisionData =
        &mut CollisionData::new(number_of_cells, region_temp);
    // Compute reasonable timestep
    let delta_x: f64 = 2. / (max(max(num_x, num_y), num_z) as f64);
    let delta_t: f64 = 0.1 * delta_x / (mean_velocity + region_temp);
    // compute nearest power of 2 timesteps
    let run_time: f64 = if time_step < 0. {
        8. / (mean_velocity + region_temp)
    } else {
        time_step
    } / delta_t;
    let run_time: usize = 1 << (run_time.log2().ceil()) as usize;
    println!("time step {:?}", time_step);
    println!("end of time {:?}", run_time);
    let mut num_sample = 0;
    //
    let sample_reset = run_time / 4;
    println!("starting main loop");
    //---------------------------HOT LOOP ----------------------------------------
    //----------------------------------------------------------------------------
    (0..run_time).for_each(|n| {
        // Add particles at inflow boundaries
        initialize_boundaries(
            &mut particles,
            num_x,
            num_y,
            num_z,
            mean_velocity,
            region_temp,
            mean_particles_per_cell,
        );
        // Move particles
        //particles.update_positions();
        move_particles_with_bcs(&mut particles, plate_x, plate_height, plate_width, delta_t);

        num_sample += 1;
        // If time to reset cell samples, reinitialize data
        if n % sample_reset == 0 {
            initialize_sample(cell_data);
            // Remove any particles that are now outside of boundaries
            remove_outside_particles(&mut particles, number_of_cells);
            num_sample = 0
        }

        // Compute cell index for particles based on their current
        // locations
        //NOTE: these task should
        rayon::join(
            || index_particles(&mut particles, &sample_data_sender, num_x, num_y, num_z),
            || sample_particles(cell_data, &sample_data_reciever),
        );

        rayon::join(
            || {
                calc_collisions(
                    &collision_idx_sender,
                    collision_data,
                    cell_data,
                    num_sample,
                    cell_vol,
                    delta_t,
                )
            },
            || collide_particles(&collision_idx_receiver, &mut particles),
        );
        println!("finished iter i");
    });
    println!("finished");
    //-------------------------WRITE RESULTS--------------------------------------
    //----------------------------------------------------------------------------

    // Write out final particle data
    write_particle_data(String::from("cells.dat"), particles).unwrap();
}

fn calc_collisions(
    collision_data_sender: &CollisionIndexSender,
    collision_data: &mut CollisionData,
    cell_data: &mut Vec<CellSample>,
    num_sample: i32,
    cell_vol: f64,
    delta_t: f64,
) {
    let number_of_cells = cell_data.len();
    //1. Visit cells and determine chance that particles can collide

    for i in 0..number_of_cells {
        // Compute a number of particles that need to be selected for
        // collision tests
        let (mut current_max_rate, mut remainder) = collision_data.get_data(i);
        let num_selections = (cell_data[i].members.len() as f64
            + cell_data[i].mean_particles
            + molecules_per_particle)
            / (cell_vol + remainder);
        let selection_count = num_selections as usize;

        remainder = num_selections - selection_count as f64;
        if selection_count > 0 {
            if cell_data[i].members.len() < 2 {
                collision_data.set_remainder(i, remainder + num_selections);
                continue;
            }
            // Select nselect particles for possible collision
            let mut rng = &mut rand::thread_rng();
            //2. Sample particles for collision and perform collision
            //  collision is a random process where the result of
            //  a collision conserves momentum and kinetic energy
            let collision_indices = (0..selection_count)
                .into_iter()
                .map(|c| -> (usize, usize) {
                    // select two points in the cell
                    let cell_indices = rand::seq::index::sample(rng, cell_data[i].members.len(), 2);
                    let p1 = cell_data[c].members[*(&cell_indices.index(0))];
                    let p2 = cell_data[c].members[*(&cell_indices.index(1))];
                    (p1, p2)
                })
                .collect::<Vec<(usize, usize)>>();

            collision_data_sender
                .send((collision_indices, current_max_rate))
                .unwrap();
            //collision_data.max_collision_rate[i] =
            //collide(current_max_rate, collision_indices, particles, rng);
        }
        collision_data.set_remainder(i, remainder);
    }
}

fn collide_particles(
    collision_index_receiver: &CollisionIndexReceiver,
    particles: &mut [Particle],
) -> f64 {
    let (collision_indices, mut current_max) = collision_index_receiver.recv().unwrap();
    let particle_arr = ParticlePointer(particles.as_mut_ptr());
    current_max = collision_indices
        .into_par_iter()
        .filter_map(move |(p1, p2)| {
            let mut rng = rand::thread_rng();
            let particle_1 = unsafe { &mut *{ particle_arr }.0.add(p1) };
            let particle_2 = unsafe { &mut *{ particle_arr }.0.add(p2) };
            let relative_velocity = particle_1.velocity - particle_2.velocity;
            let relative_velocity_norm = relative_velocity.norm();
            let collision_rate = sigma_k * relative_velocity_norm;
            let threshold = if current_max.lt(&collision_rate) {
                1.
            } else {
                collision_rate / current_max
            };
            let r: f64 = rng.gen();
            if r < threshold {
                // Collision Accepted, adjust particle velocities
                // Compute center of mass velocity, vcm
                let center_of_mass = 0.5 * (particle_1.velocity + particle_2.velocity);
                // Compute random perturbation that conserves momentum
                let perturbation = random_direction() * relative_velocity_norm;
                particle_1.velocity = center_of_mass + 0.5 * perturbation;
                particle_1.velocity = center_of_mass - 0.5 * perturbation;

                if particle_1.particle_type != ParticleType::Inflow
                    || particle_2.particle_type != ParticleType::Inflow
                {
                    let maxtype = particle_1.particle_type.max(particle_2.particle_type);
                    particle_1.particle_type =
                        particle_1.particle_type.max(particle_1.particle_type);
                    particle_2.particle_type =
                        particle_2.particle_type.max(particle_2.particle_type);
                }
            }
            if current_max.lt(&collision_rate) {
                Some(collision_rate)
            } else {
                None
            }
        })
        .collect::<Vec<f64>>()
        .into_par_iter()
        .reduce(|| 0.0, |a, b| (a.max(b)));
    return current_max;
}

// Initialize the sampled cell variables to zero
fn initialize_sample(cell_data: &mut Vec<CellSample>) {
    for sample in cell_data {
        sample.reset();
    }
}

fn sample_particles(cell_data: &mut Vec<CellSample>, update_receiver: &SampleUpdateReceiver) {
    update_receiver.recv().into_iter().for_each(
        |(cell_membership, particle_idx, particle_velocity)| {
            cell_data[cell_membership].members.push(particle_idx);
        },
    );
}

///cell membership goes like this:
/// 1. abs(x) starting from the origin,
///     - 0 and negative numbers are even
///     - positive numbers are odd
/// 2. y just goes up from the bottom,
/// 3. z from front to back, so
/// if you have a 2 x 4 x 3 grid it's members would be
/// 6 14 22    7 15 33
/// 4 12 20    5 13 21
/// 2 10 18    3 11 19
/// 0 8  16    1 9  17
fn index_particles(
    particles: &mut [Particle],
    sample_sender: &SampleUpdateSender,
    num_x: usize,
    num_y: usize,
    num_z: usize,
) {
    let num_cells = num_x * num_y * num_z;
    //assuming number of cells must be even in the x direction
    let half_x = num_x / 2;
    let grid_size = num_y * num_z;
    let z_mult = num_y * 2;

    let dy: f64 = 2. / num_y as f64;
    let dz: f64 = 2. / num_y as f64;

    particles
        .into_par_iter()
        .enumerate()
        .for_each(|(i, particle)| {
            let y_offset = ((particle.position.y + 1.0 / dy).floor() as usize).min(num_y - 1) * 2;
            let z_offset =
                ((particle.position.z + 1.0 * dz).floor() as usize).min(num_z - 1) * (num_y * 2);
            //figure out where a particle belongs based of it's location along the x-axis
            let cell_membership: usize = match particle.position {
                //scenario one: it's left of 0 or 0
                p if (-1.0..=0.0).contains(&p.x) => {
                    let x_offset = (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
                    x_offset + y_offset + z_offset
                }
                //scenario two: it's right of zero or one
                p if (0.0..=1.).contains(&p.x) => {
                    let x_offset = (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
                    x_offset + y_offset + z_offset + 1
                }
                //scenario three: it's no longer relevant
                _ => num_cells,
            };
            particle.parent_cell = cell_membership;
            if cell_membership < num_cells {
                sample_sender
                    .send((cell_membership, i, particle.velocity))
                    .unwrap();
                //cell_data[cell_membership].members.push(i);
            }
        });
}

fn update_cell_sample(
    cell_data: &mut [CellSample],
    cell_membership: usize,
    i: usize,
    particle_velocity: na::Vector3<f64>,
) {
    cell_data[cell_membership].members.push(i);
    cell_data[cell_membership].total_velocity += particle_velocity;
    cell_data[cell_membership].total_kinetic_energy += 0.5 * particle_velocity.magnitude_squared();
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///switch to two dimensional mode (only one cell in z directions)
    #[arg(long = "2d")]
    two_dim: bool,
    ///Number of cells in x direction
    #[arg(long, long = "ni", default_value_t = 32)]
    num_x: usize,
    ///Number of cells in y direction
    #[arg(long, long = "nj", default_value_t = 32)]
    num_y: usize,
    ///Number of cells in z direction
    #[arg(long, long = "nk", default_value_t = 32)]
    num_z: usize,
    ///Mean Particles per Cell in simulation, adjust number of virtual
    /// particles to meet this target for the inflow
    #[arg(long, long = "mppc", default_value_t = 10)]
    mean_particles_per_cell: i32,
    ///Mach number or ratio of mean atom velocity versus mean thermal velocity
    #[arg(long, long = "mach", default_value_t = 20.)]
    mach_number: f64,
    ///Density of the incoming flow
    #[arg(long, long = "density", default_value_t = 1e30)]
    density: f64,
    //x-location of plate
    #[arg(long, long = "px", default_value_t = -0.25)]
    plate_x: f64,
    ///y height of plate
    #[arg(long, long = "platedy", default_value_t = 0.25)]
    plate_height: f64,
    ///z width of plate
    #[arg(long, long = "platedz", default_value_t = 0.5)]
    plate_width: f64,
    ///simulation time step size (usually computed from the above parameters)
    #[arg(long, long = "time", default_value_t = -1.)]
    time_step: f64, //not usize?
}

fn write_particle_data(file_name: String, particles: Vec<Particle>) -> std::io::Result<()> {
    let mut write_buf = Vec::new();
    let ptype_int: i32;
    for particle in particles {
        write!(write_buf, "{}\n", particle).unwrap()
    }
    let file_path = Path::new(&file_name);
    let mut particle_file = File::create(file_name).unwrap();
    particle_file.write(&write_buf).unwrap();
    Ok(())
}
