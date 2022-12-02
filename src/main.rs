//#![feature(drain_filter)]
//TODO: remove this once I have some working code

mod particles;
mod util;

extern crate nalgebra as na;

#[macro_use]
use crate::util::*;
use clap::Parser;
use crossbeam::channel::{self, unbounded, Receiver, Sender};

//use lazy_static::lazy_static;
use once_cell::sync::Lazy;
use particles::{CellPointer, MaxPointer, Particle, ParticlePointer, ParticleType};
use rand::prelude::*;

use rayon::prelude::*;

use std::cmp::max;

use std::fs::File;
use std::io::Write;

use std::path::Path;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Instant;
//const CHONKY: usize = 8;
use once_cell::sync::OnceCell;
const SIGMA_K: f64 = 1.000000003e-28;

///command line arguments
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
    mean_particles_per_cell: usize,
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

fn main() {
    //-----------------------------INITIALIZATION-----------------------------
    //------------------------------------------------------------------------
    let num_cores = std::thread::available_parallelism().unwrap();
    let mean_velocity: f64 = 1.; //switched to f32 to avoid repetitive casting

    let molecules_per_particle: f64 = 1e27;
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
    let molecules_per_particle = density * cell_vol / (mean_particles_per_cell as f64);
    //create channels to share between threads

    // Create simulation data structures
    //let mut particles: Particles = Particles::new(collision_idx_receiver, sample_data_sender);
    let mut particles: Vec<Particle> = Vec::new();
    let mut cell_data: Vec<CellSample> = Vec::with_capacity(number_of_cells);
    println!("{:?}", number_of_cells);
    let mut collision_remainder: Vec<f64> = Vec::with_capacity(number_of_cells);
    let mut collision_max_rate: Vec<AtomicU64> = Vec::with_capacity(number_of_cells);
    //go ahead and initialize the cell_data
    for _ in 0..number_of_cells {
        collision_max_rate.push(AtomicU64::new((SIGMA_K * region_temp).to_bits()));
        collision_remainder.push(rand::random::<f64>());
        cell_data.push(CellSample::new());
    }

    //let mut collision_data: CollisionData = CollisionData::new(number_of_cells, region_temp);

    // store kinetic_
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
    println!("run time {:?}", run_time);
    let mut num_sample = 0;
    //
    let sample_reset = run_time / 4;
    println!("starting main loop");
    //---------------------------HOT LOOP ----------------------------------------
    //----------------------------------------------------------------------------
    (0..run_time).for_each(|n| {
        println!("initializaing boundaries");
        rayon::scope(|s| {
            // Add particles at inflow boundaries
            rayon::join(
                || {
                    initialize_boundaries(
                        &mut particles,
                        num_x,
                        num_y,
                        num_z,
                        mean_velocity,
                        region_temp,
                        mean_particles_per_cell,
                    )
                },
                || clear_cell_members(&mut cell_data),
            );
        });
        // Move particles
        println!("updating positions");
        let start = Instant::now();
        //particles.update_positions();
        move_particles_with_bcs(&mut particles, plate_x, plate_height, plate_width, delta_t);
        println!("update took {:?}", start.elapsed());
        num_sample += 1;
        // If time to reset cell samples, reinitialize data
        if n % sample_reset == 0 {
            initialize_sample(&mut cell_data);
            // Remove any particles that are now outside of boundaries
            remove_outside_particles(&mut particles, number_of_cells);
            num_sample = 0
        }
        println!("starting join index/sample");
        let start = Instant::now();
        // Compute cell index for particles based on their current
        // locations
        //NOTE: these task should
        rayon::scope(|s| {
            let (sample_sender, sample_receiver) = unbounded();
            rayon::join(
                || {
                    index_particles(&mut particles, &sample_sender, num_x, num_y, num_z);
                    drop(sample_sender)
                },
                || update_sample(&mut cell_data, &sample_receiver.clone()),
            )
        });
        println!("index/sample took {:?}", start.elapsed());
        //println!("members {:?}", cell_data);
        println!("starting join collisions");
        let start = Instant::now();
        rayon::scope(|s| {
            let (collision_idx_sender, collision_idx_receiver) = unbounded();
            rayon::join(
                || {
                    calc_collisions(
                        &collision_idx_sender,
                        &mut collision_remainder,
                        &mut cell_data,
                        molecules_per_particle,
                        num_sample,
                        cell_vol,
                        delta_t,
                    );
                    drop(collision_idx_sender)
                },
                || {
                    collide_particles(
                        &collision_idx_receiver,
                        &mut particles,
                        &mut collision_max_rate,
                    )
                },
            );
        });
        println!("collisions took {:?}", start.elapsed());
        println!("finished iter {:?}", n);
    });
    println!("finished");
    //-------------------------WRITE RESULTS--------------------------------------
    //----------------------------------------------------------------------------

    // Write out final particle data
    write_particle_data(String::from("cells.dat"), particles).unwrap();
}

///collection of data for cells, total kinetic energy and velocity have been removed
#[derive(Debug)]
struct CellSample {
    mean_particles: f64,
    //total_velocity: na::Vector3<f64>,
    //total_kinetic_energy: f64,
    members: Vec<usize>,
    //members: Arc<Mutex<Vec<usize>>>,
}

impl CellSample {
    fn new() -> Self {
        Self {
            mean_particles: 0.,
            //Turns out this is unnecessary
            //total_velocity: na::Vector3::zeros(),
            //total_kinetic_energy: 0.0,
            //Note: balance of reads to writes is about even, so use mutex
            //https://stackoverflow.com/questions/50704279/when-or-why-should-i-use-a-mutex-over-an-rwlock
            //members: Arc::new(Mutex::new(Vec::new())),
            members: Vec::new(),
        }
    }
    fn reset(&mut self) {
        self.mean_particles = 0.;
        //self.total_kinetic_energy = 0.;
        //https://stackoverflow.com/a/56114369/11019565
        //self.total_velocity.x = 0.;
        //self.total_velocity.y = 0.;
        //self.total_velocity.z = 0.;
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
    mean_particles_per_cell: usize,
) {
    // TODO: I might be able to avoid resizing the particle containing
    // data structure by running the functions
    // move_particles and remove_outside_particles (or similar function)
    //directly on the generated particles
    let dx = 2. / (num_x as f64);
    let dy = 2. / (num_y as f64);
    let dz = 2. / (num_z as f64);
    let cell: OnceCell<Vec<(usize, usize)>> = OnceCell::new();
    let outer_range = cell.get_or_init(|| gen_outer_range(num_y, num_z));
    //let tmp=Vec<Particle> = Vec::new()
    particles.par_extend(
        outer_range
            .into_par_iter()
            .flat_map(|(j, k)| -> Vec<Particle> {
                let current_x = -1. - dx;
                let current_y = -1. + *j as f64 * dy;
                let current_z = -1. + *k as f64 * dz;
                println!("d {} {} {}", dx, dy, dz);
                println!("curr {} {} {}", current_x, current_y, current_z);
                (0..mean_particles_per_cell)
                    .into_par_iter()
                    //.chunks(CHONKY)
                    .map(|_| -> Particle {
                        let mut rng = rand::thread_rng();
                        let position = na::Vector3::<f64>::new(
                            current_x + rng.gen::<f64>() * dx,
                            current_y + rng.gen::<f64>() * dy,
                            current_z + rng.gen::<f64>() * dz,
                        );
                        println!("{}", position);
                        let mut velocity =
                            random_velocity(&mut rng, region_temp) * random_direction(&mut rng);
                        velocity.x += mean_velocity;
                        Particle::new(position, velocity)
                    })
                    .collect()
            })
            //.filter(|p| (-1.0..=1.0).contains(&p.velocity.x))W
            .collect::<Vec<Particle>>(),
    );
    println!("particle count: {:?}", particles.len());
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
    particles.par_sort_unstable_by_key(|p| (p.parent_cell));
    if let Some(cutoff) = particles.iter().rposition(|p| p.parent_cell < num_cells) {
        particles.truncate(cutoff + 1);
    }
    //let _ = particles.drain_filter(|p| p.parent_cell < num_cells);
}

//collideParticles() : collision is a random process where the result of
//a collision conserves momentum and kinetic energy
// time-step is complete, continues

fn calc_collisions(
    collision_data_sender: &CollisionIndexSender,
    collision_remainder: &mut Vec<f64>,
    cell_data: &mut Vec<CellSample>,
    molecules_per_particle: f64,
    num_sample: i32,
    cell_vol: f64,
    delta_t: f64,
) {
    let number_of_cells = cell_data.len();
    //1. Visit cells and determine chance that particles can collide
    println!("starting calc collisions");
    collision_remainder
        .into_par_iter()
        .enumerate()
        .for_each(|(i, remainder)| {
            // Compute a number of particles that need to be selected for
            // collision tests

            //let members = cell_data[i].members.lock().unwrap();
            let members = &cell_data[i].members;

            let num_selections =
                (members.len() as f64 + cell_data[i].mean_particles + molecules_per_particle)
                    / (cell_vol + *remainder);
            let selection_count = num_selections as usize;
            println!("selection count {:?}", selection_count);
            *remainder = num_selections - selection_count as f64;
            if selection_count > 0 {
                if members.len() < 2 {
                    *remainder += num_selections;
                } else {
                    // Select nselect particles for possible collision

                    //2. Sample particles for collision and perform collision
                    //  collision is a random process where the result of
                    //  a collision conserves momentum and kinetic energy
                    let collision_indices = (0..selection_count)
                        .into_par_iter()
                        .map(|c| -> (usize, usize) {
                            let mut rng = &mut rand::thread_rng();
                            // select two points in the cell
                            let cell_indices = rand::seq::index::sample(rng, members.len(), 2);
                            let p1 = members[*(&cell_indices.index(0))];
                            let p2 = members[*(&cell_indices.index(1))];
                            (p1, p2)
                        })
                        .chunks(8)
                        .for_each(|chunk| collision_data_sender.send(chunk).unwrap());

                    //collision_data.max_collision_rate[i] =
                    //collide(current_max_rate, collision_indices, particles, rng);
                }
            }
        });
    println!("end_reached");
}

fn collide_particles(
    collision_index_receiver: &CollisionIndexReceiver,
    particles: &mut [Particle],
    collision_max_rate: &mut [AtomicU64],
) {
    let particle_arr = ParticlePointer(particles.as_mut_ptr());
    let current_max_pointer = MaxPointer(collision_max_rate.as_mut_ptr());
    let num_cells = collision_max_rate.len();
    collision_index_receiver
        .into_iter()
        .for_each(|collision_indices| {
            //println!("got collision indices: {:?}", collision_indices);

            let current_max_atomic = get_max!(collision_indices, particle_arr, current_max_pointer);
            let current_max = f64::from_bits(current_max_atomic.load(Ordering::Relaxed));
            collision_indices.into_par_iter().for_each(move |(p1, p2)| {
                let mut rng = rand::thread_rng();
                let particle_1 = unsafe { &mut *{ particle_arr }.0.add(p1) };
                let particle_2 = unsafe { &mut *{ particle_arr }.0.add(p2) };
                //let cell_index = particle_1.parent_cell;
                let relative_velocity = particle_1.velocity - particle_2.velocity;
                let relative_velocity_norm = relative_velocity.norm();
                let collision_rate: &mut f64 = &mut (SIGMA_K * relative_velocity_norm);
                let threshold = if current_max.lt(&collision_rate) {
                    1.
                } else {
                    *collision_rate / current_max
                };
                let r: f64 = rng.gen();
                if r < threshold {
                    // Collision Accepted, adjust particle velocities
                    // Compute center of mass velocity, vcm
                    let center_of_mass = 0.5 * (particle_1.velocity + particle_2.velocity);
                    // Compute random perturbation that conserves momentum
                    let perturbation = random_direction(&mut rng) * relative_velocity_norm;
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
                    current_max_atomic.store(collision_rate.to_bits(), Ordering::Relaxed)
                }
            });
        });
}

// Initialize the sampled cell variables to zero
fn initialize_sample(cell_data: &mut Vec<CellSample>) {
    for sample in cell_data {
        sample.reset();
    }
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

    println!("indexing particles start");
    particles
        .into_par_iter()
        .enumerate()
        .map(|(i, particle)| {
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
            //println!("send it");
            particle.parent_cell = cell_membership;
            (cell_membership, i)
        })
        .chunks(8)
        .for_each(|chunk| sample_sender.send(chunk).unwrap());
    //sample_sender.send((num_cells, num_cells)).unwrap();
}

fn update_sample(cell_data: &mut Vec<CellSample>, sample_reciever: &SampleUpdateReceiver) {
    let cell_ptr = CellPointer(cell_data.as_mut_ptr());
    let num_cells = cell_data.len();
    sample_reciever.into_iter().for_each(|chunk| {
        //chunk.into_par_iter().chunks(4).for_each(|chunk| {
        chunk.into_iter().for_each(|(cell_idx, particle_idx)| {
            println!("receiving {:?}", cell_idx);
            if cell_idx < num_cells {
                // let cell = unsafe { &*{ cell_ptr }.0.add(cell_idx) };
                // let guard = Arc::clone(&cell.members);

                // let mut members = guard.lock().unwrap();
                if cell_data[cell_idx].members.len() != 0 {
                    println!("yup")
                }
                // members.push(particle_idx);
                cell_data[cell_idx].members.push(particle_idx)
                //println!("members {:?}", members);
            }
        })
        //});
        //});
    });
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

fn clear_cell_members(cell_data: &mut Vec<CellSample>) {
    cell_data.into_par_iter().for_each(|cell| {
        //let mut members = cell.members.lock().unwrap();
        //members.clear();
        cell.members.clear()
    });
}
