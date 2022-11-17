extern crate nalgebra as na;
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

//collideParticles() : collision is a random process where the result of
//a collision conserves momentum and kinetic energy
// time-step is complete, continues
fn main() {
    loop {}
}
