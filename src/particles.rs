use std::{
    fmt::{Display, Formatter},
    ops::{Deref, DerefMut},
    sync::Arc,
};

use crossbeam::channel::{unbounded, Receiver, Sender};
use nalgebra::Vector3;
use rand::prelude::*;

use crate::{random_direction, random_velocity, util::CollisionIndexReceiver, CellSample, SIGMA_K};

#[derive(PartialOrd, Ord, PartialEq, Eq, Copy, Clone, Debug)]
pub(crate) enum ParticleType {
    Inflow,
    Ricochet,
    Thud,
}

impl ParticleType {
    fn to_num(self) -> i32 {
        match self {
            ParticleType::Inflow => 0,
            ParticleType::Ricochet => 1,
            ParticleType::Thud => 2,
        }
    }
}

//see https://stackoverflow.com/a/70848420/11019565
//thin wrapper over pointer to make it Send/Sync
#[derive(Copy, Clone)]
pub(crate) struct ParticlePointer(pub(crate) *mut Particle);
unsafe impl Send for ParticlePointer {}
unsafe impl Sync for ParticlePointer {}

#[derive(Copy, Clone)]
pub(crate) struct MaxPointer(pub(crate) *mut f64);
unsafe impl Send for MaxPointer {}
unsafe impl Sync for MaxPointer {}

#[derive(Copy, Clone)]
pub(crate) struct CellPointer(pub(crate) *const CellSample);
unsafe impl Send for CellPointer {}
unsafe impl Sync for CellPointer {}

// let indices = [1, 4, 7, 8];
// let mut arr = [1u32, 2, 3, 4, 5, 6, 7, 8, 9, 10];
// let arr_ptr = Pointer(arr.as_mut_ptr());

// indices.into_par_iter().for_each(move |x| {
//     // safety:
//     // * `indices` must be unique and point inside `arr`
//     // * `place` must not leak outside the closure
//     // * no element of `array` that is in `indices` may be accessed by
//     //   some other thread while this is running
//     let place = unsafe { &mut *{arr_ptr}.0.add(x) };
//     *place = some_function(x);
// });

#[derive(Copy, Clone, Debug)]
pub struct Particle {
    pub(crate) velocity: na::Vector3<f64>,
    pub(crate) position: na::Vector3<f64>,
    pub(crate) particle_type: ParticleType,
    pub(crate) parent_cell: usize,
}

impl Display for Particle {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {} {} {}\n",
            self.position[0],
            self.position[1],
            self.position[2],
            self.particle_type.to_num()
        )
    }
}

impl Particle {
    pub(crate) fn new(velocity: na::Vector3<f64>, position: na::Vector3<f64>) -> Self {
        Self {
            velocity,
            position,
            particle_type: ParticleType::Inflow,
            parent_cell: 0,
        }
    }
}

pub(crate) struct Particles {
    //receiver: ThingBuf,
    collision_index_receiver: CollisionIndexReceiver,
    sample_sender: Sender<(usize, usize, Vector3<f64>)>,
    // velocities: Vec<na::Vector3>,
    // pub(crate) positions: Vec<na::Vector3>,
    // types: Vec<ParticleType>,
    // parent_cells: Vec<usize>,
    vec: Vec<Particle>,
}

// impl Particles {
//     pub(crate) fn new(
//         collision_index_receiver: Receiver<(Vec<(usize, usize)>, f64)>,
//         sample_sender: Sender<(usize, usize, Vector3<f64>)>,
//     ) -> Self {
//         Self {
//             collision_index_receiver,
//             sample_sender,
//             vec: Vec::new(),
//         }
//     }

//     fn collide(&mut self) {
//         let mut rng = rand::thread_rng();
//         let vec = self.vec;
//         let particles = ParticlePointer(self.vec.as_mut_ptr());
//         let (collision_indices, mut current_max) = self.collision_index_receiver.recv().unwrap();
//         collision_indices.into_iter().for_each(|(p1, p2)| {
//             let particle_1 = unsafe { &mut *{ particles }.0.add(p1) };
//             let particle_2 = unsafe { &mut *{ particles }.0.add(p2) };

//             let relative_velocity = particle_1.velocity - particle_2.velocity;
//             let relative_velocity_norm = relative_velocity.norm();
//             let collision_rate = sigma_k * relative_velocity_norm;
//             let threshold = if collision_rate > current_max {
//                 //TODO: if different thread send update here
//                 current_max = collision_rate;
//                 1.
//             } else {
//                 collision_rate / current_max
//             };
//             let r: f64 = rng.gen();
//             if r < threshold {
//                 // Collision Accepted, adjust particle velocities
//                 // Compute center of mass velocity, vcm
//                 let center_of_mass = 0.5 * (particle_1.velocity + particle_2.velocity);
//                 // Compute random perturbation that conserves momentum
//                 let perturbation = random_direction() * relative_velocity_norm;
//                 particle_1.velocity = center_of_mass + 0.5 * perturbation;
//                 particle_1.velocity = center_of_mass - 0.5 * perturbation;
//                 //seems weird but I think this involves fewest comparisons
//                 if particle_1.particle_type == ParticleType::Inflow
//                     && particle_2.particle_type == ParticleType::Inflow
//                 {
//                     return;
//                 } else {
//                     particle_1.particle_type = ParticleType::Ricochet;
//                     particle_2.particle_type = ParticleType::Ricochet;
//                 }
//             }
//         });
//     }

//     // fn new(
//     //     collision_index_receiver: unbounded::Receiver<(Vec<(usize, usize)>, f64)>,
//     //     sample_sender: unbounded::Sender<(usize, usize, Vector3)>,
//     // ) -> Self {
//     //     Self {
//     //         collision_index_receiver,
//     //         sample_sender,
//     //         velocities: Vec::new(),
//     //         positions: Vec::new(),
//     //         types: Vec::new(),
//     //         parent_cells: Vec::new(),
//     //     }
//     // }
//     fn initialize_boundaries(
//         &mut self,
//         num_x: usize,
//         num_y: usize,
//         num_z: usize,
//         mean_velocity: f64,
//         region_temp: f64,
//         mean_particles_per_cell: i32,
//     ) {
//         let dx = 2. / num_x as f64;
//         let dy = 2. / num_y as f64;
//         let dz = 2. / num_z as f64;
//         self.vec
//             .extend((0..num_y).into_iter().flat_map(|j| -> Vec<Particle> {
//                 (0..num_z)
//                     .into_iter()
//                     .flat_map(|k| -> Vec<Particle> {
//                         let current_x = -1. - dx;
//                         let current_y = -1. + j as f64 * dy;
//                         let current_z = -1. + k as f64 * dz;

//                         (0..mean_particles_per_cell)
//                             .into_iter()
//                             .map(|i| -> Particle {
//                                 let position = na::Vector3::<f64>::new(
//                                     current_x + rand::random::<f64>() * dx,
//                                     current_y + rand::random::<f64>() * dy,
//                                     current_z + rand::random::<f64>() * dz,
//                                 );
//                                 let velocity = random_velocity(region_temp) * random_direction();
//                                 Particle::new(position, velocity)
//                             })
//                             .collect::<Vec<Particle>>()
//                     })
//                     .collect()
//             }));
//     }

//     fn update_postions(&mut self, delta_t: f64, plate_x: f64, plate_dy: f64, plate_dz: f64) {
//         self.vec.iter().for_each(|mut particle| {
//             let current_position = particle.position;
//             let mut new_position = current_position + particle.velocity * delta_t;

//             if (current_position[0] < plate_x) && (new_position[0] > plate_x)
//                 || (current_position[0] > plate_x) && (new_position[0] < plate_x)
//             {
//                 let t = (current_position[0] - plate_x) / (current_position[0] - new_position[0]);
//                 //create a particle representing the potential position
//                 let tmp = current_position * (1. - t)
//                     + ((current_position + delta_t * particle.velocity) * t);
//                 if (-plate_dy..plate_dy).contains(&tmp[1])
//                     && (-plate_dz..plate_dz).contains(&tmp[2])
//                 {
//                     new_position.x -= 2. * new_position.x - plate_x;
//                     particle.velocity[0] = -particle.velocity[0];
//                     particle.particle_type = ParticleType::Thud;
//                 }
//             }

//             if (-1.0..1.).contains(&new_position.y) {
//                 new_position.y -= 2.0;
//             }
//             if (-1.0..1.).contains(&new_position.z) {
//                 new_position.z -= 2.0;
//             }
//             particle.position = new_position;
//         });
//     }

//     /// any particles outside of the bounds of the grid need to be discarded as they cannot be indexed.
//     /// Since this can happen only at the x=-1 or x=1 boundaries, we only need to check the x coordinate
//     pub(crate) fn filter_particles(&mut self, num_cells: usize) {
//         // self.vec.sort_by_key(|p| (p.parent_cell));
//         // if let Some(cutoff) = particles.iter().rposition(|p| p.parent_cell < num_cells) {
//         //     particles.truncate(cutoff + 1);
//         // }
//         let _ = self.vec.drain_filter(|p| p.parent_cell < num_cells);
//     }

//     ///cell membership goes like this:
//     /// 1. abs(x) starting from the origin,
//     ///     - 0 and negative numbers are even
//     ///     - positive numbers are odd
//     /// 2. y just goes up from the bottom,
//     /// 3. z from front to back, so
//     /// if you have a 2 x 4 x 3 grid it's members would be
//     /// 6 14 22    7 15 33
//     /// 4 12 20    5 13 21
//     /// 2 10 18    3 11 19
//     /// 0 8  16    1 9  17
//     fn index(&mut self, num_x: usize, num_y: usize, num_z: usize) {
//         let num_cells = num_x * num_y * num_z;
//         //assuming number of cells must be even in the x direction
//         let half_x = num_x / 2;
//         let grid_size = num_y * num_z;
//         let z_mult = num_y * 2;

//         let dy: f64 = 2. / num_y as f64;
//         let dz: f64 = 2. / num_y as f64;

//         for (i, particle) in self.vec.iter_mut().enumerate() {
//             let y_offset = ((particle.position.y + 1.0 / dy).floor() as usize).min(num_y - 1) * 2;
//             let z_offset =
//                 ((particle.position.z + 1.0 * dz).floor() as usize).min(num_z - 1) * (num_y * 2);
//             //figure out where a particle belongs based of it's location along the x-axis
//             let cell_membership: usize = match particle.position {
//                 //scenario one: it's left of 0 or 0
//                 p if (-1.0..=0.0).contains(&p.x) => {
//                     let x_offset = (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
//                     x_offset + y_offset + z_offset
//                 }
//                 //scenario two: it's right of zero or one
//                 p if (0.0..=1.).contains(&p.x) => {
//                     let x_offset = (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
//                     x_offset + y_offset + z_offset + 1
//                 }
//                 //scenario three: it's no longer relevant
//                 _ => num_cells,
//             };
//             particle.parent_cell = cell_membership;
//             if cell_membership < num_cells {
//                 self.sample_sender
//                     .send((cell_membership, i, particle.velocity));
//             }
//         }
//     }
// }

// impl IntoIterator for Particles {
//     type Item = Particle;
//     type IntoIter = <Vec<Particle> as IntoIterator>::IntoIter;
//     fn into_iter(self) -> Self::IntoIter {
//         self.vec.into_iter()
//     }
// }

// We deref to slice so that we can reuse the slice impls
// impl Deref for Particles {
//     type Target = [Particle];

//     fn deref(&self) -> &[Particle] {
//         &self.vec[..]
//     }
// }
// impl DerefMut for Particles {
//     // type Output = [Particle];

//     fn deref_mut(&mut self) -> &mut [Particle] {
//         &mut self.vec[..]
//     }
// }

// impl std::ops::Deref for Particles {
//     type Target = Vec<Particle>;

//     // fn deref(&self) -> &Self::Target {
//     //     &self.particles
//     // }
// }
