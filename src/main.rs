#![allow(dead_code, unused_variables, unused_imports)]
//TODO: remove this once I have some working code

extern crate nalgebra as na;
use clap::{Arg, Parser};
use num::{Bounded, Integer};

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
use rand::thread_rng;

const COLLISION_SIZE: f32 = 1e-28;
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
    maxCollisionRate: f64,
    // Non-integral fraction of collisions that remain to be performed
    // and are carried over into the next timestep
    collisionRemainder: f64,
}

struct cellSample {
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
fn collideParticles() {
    //1. Visit cells and determine chance that particles can collide
    //2. Sample particles for collision and perform collision
    //  collision is a random process where the result of
    //  a collision conserves momentum and kinetic energy
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///switch to two dimensional mode (only one cell in z directions)
    #[arg(long = "2d")]
    two_dim: bool,
    ///Number of cells in x direction
    #[arg(long, long = "ni", default_value_t = 32)]
    num_x: i32,
    ///Number of cells in y direction
    #[arg(long, long = "nj", default_value_t = 32)]
    num_y: i32,
    ///Number of cells in z direction
    #[arg(long, long = "nk", default_value_t = 32)]
    num_z: i32,
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
    #[arg(long, long = "nk", default_value_t = 0)]
    time_step: usize,
}
//collideParticles() : collision is a random process where the result of
//a collision conserves momentum and kinetic energy
// time-step is complete, continues
fn main() {
    let particles_per_particle = 1e27;
    let mean_velocity = 1;
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
    let region_temp = mean_velocity as f32 / mach_number;
    // Compute number of molecules a particle represents
    // Create simulation data structures
    loop {}
}
