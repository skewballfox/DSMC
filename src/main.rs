#![allow(dead_code, unused_variables, unused_imports, unused_mut)]
//TODO: remove this once I have some working code

extern crate nalgebra as na;
use clap::{Arg, Parser};
use na::coordinates::X;
use num::{Bounded, Integer};
use rand::thread_rng;
use std::cmp::max;
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use std::path::Path;

const COLLISION_SIZE: f32 = 1e-28;

enum ParticleType {
    Inflow,
    Ricochet,
    Thud,
}

#[derive()]
struct Particle {
    velocity: na::Vector3<f64>,
    position: na::Vector3<f64>,
    particle_type: ParticleType,
    parent_cell: usize,
}

impl Particle {
    fn new(velocity: na::Vector3<f64>, position: na::Vector3<f64>) -> Self {
        Self {
            velocity,
            position,
            particle_type: ParticleType::Inflow,
            parent_cell: 0,
        }
    }
}

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
struct CollisionInfo {
    // Maximum collision rate seen for this cell so far in the simulation
    max_collision_rate: f64,
    // Non-integral fraction of collisions that remain to be performed
    // and are carried over into the next timestep
    collision_remainder: f64,
}

struct CellSample {
    number_of_particles: i32,
    total_velocity: na::Vector3<f64>,
    total_kinetic_energy: f64,
    start: usize,
    stop: usize,
}

impl CellSample {
    fn new() -> Self {
        Self {
            number_of_particles: 0,
            total_velocity: na::Vector3::zeros(),
            total_kinetic_energy: 0.0,
            start: 0,
            stop: 0,
        }
    }
    fn reset(&mut self) {
        self.number_of_particles = 0;
        self.total_kinetic_energy = 0.;
        //https://stackoverflow.com/a/56114369/11019565
        self.total_velocity.x = 0.;
        self.total_velocity.y = 0.;
        self.total_velocity.z = 0.;
    }
}

pub struct Cell {
    start: usize,
    stop: usize,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CellIndex {
    index: [u32; 3],
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
    for j in 0..num_y {
        for k in 0..num_z {
            let current_x = -1. - dx;
            let current_y = -1. + j as f64 * dy;
            let current_z = -1. + k as f64 * dz;

            for i in 0..mean_particles_per_cell {
                let position = na::Vector3::<f64>::new(
                    current_x + rand::random::<f64>() * dx,
                    current_y + rand::random::<f64>() * dy,
                    current_z + rand::random::<f64>() * dz,
                );
                let velocity = random_velocity(region_temp) * random_direction();
                let particle = Particle::new(position, velocity);
                particles.push(particle)
            }
        }
    }
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
    for particle in particles {
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
    }
}

/// any particles outside of the cells need to be discarded as they cannot be indexed.  
/// Since this can happen only at the x=-1 or x=1 boundaries, we only need to check the x coordinate
fn remove_outside_particles(particles: &mut Vec<Particle>, num_cells: usize) {
    particles.sort_by_key(|p| (p.parent_cell));
    if let Some(cutoff) = particles.iter().rposition(|p| p.parent_cell < num_cells) {
        particles.truncate(cutoff + 1);
    }
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
    let mean_velocity: f64 = 1.; //switched to f32 to avoid repetitive casting
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
    //number of molecules each particle represents
    let molecules_per_particle = 1e27;
    let number_of_cells = num_x * num_y * num_z;
    // Compute number of molecules a particle represents
    let cell_vol: f64 = 2. / (number_of_cells) as f64;
    // Create simulation data structures
    let mut particles: &mut Vec<Particle> = &mut Vec::new();
    let mut cell_data: &mut Vec<CellSample> = &mut Vec::with_capacity(number_of_cells);
    //go ahead and initialize the cell_data
    for i in 0..number_of_cells {
        cell_data.push(CellSample::new());
    }

    let mut collision_data: &mut Vec<CollisionInfo> = &mut Vec::with_capacity(number_of_cells);
    // Compute reasonable timestep
    let delta_x: f64 = 2. / (max(max(num_x, num_y), num_z)) as f64;
    let delta_t: f64 = 0.1 * delta_x / (mean_velocity + region_temp);
    // compute nearest power of 2 timesteps
    let end_of_time: f64 = if time_step < 0. {
        8. / (mean_velocity + region_temp)
    } else {
        time_step
    } / delta_t;
    let end_of_time: usize = 1 << ((time_step).log2().ceil()) as usize;

    let mut num_sample = 0;
    //
    let sample_reset = end_of_time / 4;
    //---------------------------HOT LOOP ----------------------------------------
    //----------------------------------------------------------------------------
    (0..end_of_time).for_each(|n| {
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
        move_particles_with_bcs(&mut particles, plate_x, plate_height, plate_width, delta_t);

        // Compute cell index for particles based on their current
        // locations
        index_particles(particles, num_x, num_y, num_z);
        // Remove any particles that are now outside of boundaries
        remove_outside_particles(particles, number_of_cells);

        num_sample += 1;
        // If time to reset cell samples, reinitialize data
        if n % sample_reset == 0 {
            initialize_sample(cell_data);
            num_sample = 0
        }

        sample_particles(cell_data, particles);
        collide_particles(
            particles,
            collision_data,
            cell_data,
            num_sample,
            cell_vol,
            delta_t,
        )
    });
    //-------------------------WRITE RESULTS--------------------------------------
    //----------------------------------------------------------------------------

    // Write out final particle data
}

fn collide_particles(
    particles: &[Particle],
    collision_data: &mut Vec<CollisionInfo>,
    cell_data: &mut Vec<CellSample>,
    num_sample: i32,
    cell_vol: f64,
    delta_t: f64,
) {
    let number_of_cells = cell_data.len();
    //1. Visit cells and determine chance that particles can collide
    //2. Sample particles for collision and perform collision
    //  collision is a random process where the result of
    //  a collision conserves momentum and kinetic energy
    for i in 0..number_of_cells {
        let n_mean = cell_data[i].number_of_particles as f64 / num_sample as f64;
    }
}

// Initialize the sampled cell variables to zero
fn initialize_sample(cell_data: &mut Vec<CellSample>) {
    for sample in cell_data {
        sample.reset();
    }
}

fn sample_particles(cell_data: &mut Vec<CellSample>, particles: &[Particle]) {
    for particle in particles {
        let sample = &mut cell_data[particle.parent_cell];
        sample.number_of_particles += 1;
        sample.total_velocity += particle.velocity;
        sample.total_kinetic_energy += 0.5 * particle.velocity.magnitude_squared();
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
fn index_particles(particles: &mut Vec<Particle>, num_x: usize, num_y: usize, num_z: usize) {
    let num_cells = num_x * num_y * num_z;
    //assuming number of cells must be even in the x direction
    let half_x = num_x / 2;
    let grid_size = num_y * num_z;
    let z_mult = num_y * 2;

    let dy: f64 = 2. / num_y as f64;
    let dz: f64 = 2. / num_y as f64;

    for particle in particles {
        let y_offset = ((particle.position.y + 1.0 / dy).floor() as usize).min(num_y - 1) * 2;
        let z_offset =
            ((particle.position.z + 1.0 * dz).floor() as usize).min(num_z - 1) * (num_y * 2);

        let cell_membership = match particle.position {
            p if (-1.0..=0.0).contains(&p.x) => {
                let x_offset = (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
                x_offset + y_offset + z_offset
            }
            p if (0.0..=1.).contains(&p.x) => {
                let x_offset = (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
                x_offset + y_offset + z_offset + 1
            }
            _ => num_cells,
        };
        particle.parent_cell = cell_membership
    }
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
    #[arg(long, long = "nk", default_value_t = 0.)]
    time_step: f64, //not usize?
}

fn write_particle_data(file_name: String, particles: Vec<Particle>) -> std::io::Result<()> {
    let mut write_buf = Vec::new();
    let ptype_int: i32;
    for particle in particles {
        let ptype_int = match particle.particle_type {
            ParticleType::Inflow => 0,
            ParticleType::Ricochet => 1,
            ParticleType::Thud => 2,
        };
        write!(
            &mut write_buf,
            "{} {} {} {}\n",
            particle.position[0], particle.position[1], particle.position[2], ptype_int
        )?;
    }
    let file_path = Path::new(&file_name);
    let mut particle_file = File::create(file_name);

    Ok(())
}
