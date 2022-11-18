#![allow(dead_code, unused_variables, unused_imports, unused_mut)]
//TODO: remove this once I have some working code

extern crate nalgebra as na;
use clap::{Arg, Parser};
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
struct Particle {
    velocity: na::Vector3<f64>,
    position: na::Vector3<f64>,
    particle_type: ParticleType,
    parent_cell: CellIndex,
}

impl Particle {
    fn new(velocity: na::Vector3<f64>, position: na::Vector3<f64>) -> Self {
        Self {
            velocity,
            position,
            particle_type: ParticleType::Inflow,
            parent_cell: CellIndex { index: [0, 0, 0] },
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
fn random_velocity(region_temp: f32) -> f64 {
    return region_temp as f64 * (-f64::max(rand::random::<f64>(), 1e-200).ln()).sqrt();
}

/// Unique identifier for a cell on a grid, represented by an index triplet on the 3D cartesian grid
///
/// The local indexing of points and edges are given as follows:
/// ```text
///  Points:                Edges:
///          7 ________ 6           _____6__
///          /|       /|         7/|       /|
///        /  |     /  |        /  |     /5 |
///    4 /_______ /    |      /__4____ /    10
///     |     |  |5    |     |    11  |     |
///     |    3|__|_____|2    |     |__|__2__|
///     |    /   |    /      8   3/   9    /
///     |  /     |  /        |  /     |  /1
///     |/_______|/          |/___0___|/
///    0          1
/// ```
/// Point with index 0 is closest to the coordinate origin.
//pub trait Index: Copy + Integer + Bounded {}
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
}
pub struct Cell {
    index: CellIndex,
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
fn initializeBoundaries(
    particles: Vec<Particle>,
    num_x: usize,
    num_y: usize,
    num_z: usize,
    mean_velocity: f32,
    region_temp: f32,
    mean_particles_per_cell: i32,
) {
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
fn moveParticlesWithBCs(particles: Vec<Particle>, delta_T: f32) {}
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
fn collide_particles() {
    //1. Visit cells and determine chance that particles can collide
    //2. Sample particles for collision and perform collision
    //  collision is a random process where the result of
    //  a collision conserves momentum and kinetic energy
}

//collideParticles() : collision is a random process where the result of
//a collision conserves momentum and kinetic energy
// time-step is complete, continues
fn main() {
    //-----------------------------INITIALIZATION-----------------------------
    //------------------------------------------------------------------------
    let mean_velocity: f32 = 1.; //switched to f32 to avoid repetitive casting
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
    let region_temp = mean_velocity / mach_number;
    //number of molecules each particle represents
    let molecules_per_particle = 1e27;
    let number_of_cells = num_x * num_y * num_z;
    // Compute number of molecules a particle represents
    let cell_vol: f64 = 2. / (number_of_cells) as f64;
    // Create simulation data structures
    let mut particles: Vec<Particle> = Vec::new();
    let mut cell_data: Vec<CellSample> = Vec::with_capacity(number_of_cells);
    let mut collision_data: Vec<CollisionInfo> = Vec::with_capacity(number_of_cells);
    // Compute reasonable timestep
    let delta_x = 2. / (max(max(num_x, num_y), num_z)) as f32;
    let delta_t = 0.1 * delta_x / (mean_velocity + region_temp);
    // compute nearest power of 2 timesteps
    let end_of_time: f32 = if time_step < 0. {
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
        initializeBoundaries(
            particles,
            num_x,
            num_y,
            num_z,
            mean_velocity,
            region_temp,
            mean_particles_per_cell,
        );
        // Move particles
        moveParticlesWithBCs(particleList, deltaT);
        // Remove any particles that are now outside of boundaries
        removeOutsideParticles(particleList);
        // Compute cell index for particles based on their current
        // locations
        indexParticles(particleList, ni, nj, nk);
        // If time to reset cell samples, reinitialize data
        if n % sample_reset == 0 {
            num_sample = 0
        }
    });
    //-------------------------WRITE RESULTS--------------------------------------
    //----------------------------------------------------------------------------

    // Write out final particle data
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
    mach_number: f32,
    ///Density of the incoming flow
    #[arg(long, long = "density", default_value_t = 1e30)]
    density: f32,
    //x-location of plate
    #[arg(long, long = "px", default_value_t = -0.25)]
    plate_x: f32,
    ///y height of plate
    #[arg(long, long = "platedy", default_value_t = 0.25)]
    plate_height: f32,
    ///z width of plate
    #[arg(long, long = "platedz", default_value_t = 0.5)]
    plate_width: f32,
    ///simulation time step size (usually computed from the above parameters)
    #[arg(long, long = "nk", default_value_t = 0.)]
    time_step: f32, //not usize?
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
    let file_path = Path::new(&file_name);
    Ok(())
}
