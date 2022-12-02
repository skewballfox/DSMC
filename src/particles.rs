use std::{
    fmt::{Display, Formatter},
    ops::{Deref, DerefMut},
    sync::{
        atomic::{AtomicBool, AtomicU64},
        Arc, Mutex,
    },
};

use crate::{
    random_direction, random_velocity,
    util::{
        create_parking_lot, gen_outer_range, grow_parking_lot, BoolPointer, CollisionIndexReceiver,
        ImmutableParentCellPointer, ParentCellPointer, ParticleTypePointer, SampleUpdateSender,
        VectorPointer,
    },
    CellSample, SIGMA_K,
};
use crossbeam::{
    channel::{unbounded, Receiver, Sender},
    sync::Parker,
};
use nalgebra::{Matrix, Vector3};
use once_cell::sync::OnceCell;
use rand::prelude::*;
use rayon::prelude::*;

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

pub struct Particles {
    //receiver: ThingBuf,
    velocities: Vec<na::Vector3<f64>>,
    positions: Vec<na::Vector3<f64>>,
    types: Vec<ParticleType>,
    parent_cells: Vec<usize>,
    keepers: Vec<bool>,
}

impl Particles {
    pub fn new() -> Self {
        Self {
            velocities: Vec::new(),
            positions: Vec::new(),
            types: Vec::new(),
            parent_cells: Vec::new(),
            keepers: Vec::new(),
        }
    }

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

    fn collide(&mut self, p1: usize, p2: usize, current_max_rate: f64, rng: &mut ThreadRng) -> f64 {
        let Particles {
            velocities,
            positions: _,
            types,
            parent_cells: _,
            keepers: _,
        } = self;
        let mut cell: OnceCell<Vec<Arc<Parker>>> = OnceCell::new();
        let mut plock = cell.get_or_init(|| create_parking_lot(self.len()));
        if plock.len() < self.len() {
            grow_parking_lot(plock, lot_size)
        }
        //println!("got collision indices: {:?}", collision_indices);
        let lock1 = plock[p1].clone();
        let lock2 = plock[p2].clone();
        lock1.park();
        lock2.park();
        let relative_velocity = velocities[p1] - velocities[p2];
        let relative_velocity_norm = relative_velocity.norm();
        let collision_rate: &mut f64 = &mut (SIGMA_K * relative_velocity_norm);
        let threshold = if current_max_rate.lt(&collision_rate) {
            1.
        } else {
            *collision_rate / current_max_rate
        };
        let r: f64 = rng.gen();
        if r < threshold {
            // Collision Accepted, adjust particle velocities
            // Compute center of mass velocity, vcm
            let center_of_mass = 0.5 * (velocities[p1] + velocities[p2]);
            // Compute random perturbation that conserves momentum
            let perturbation = random_direction(rng) * relative_velocity_norm;

            // Adjust particle velocities to reflect collision
            velocities[p1] = center_of_mass + 0.5 * perturbation;
            velocities[p2] = center_of_mass - 0.5 * perturbation;

            // Bookkeeping to track particle interactions
            if types[p1] != ParticleType::Inflow || types[p2] != ParticleType::Inflow {
                let maxtype = types[p1].max(types[p2]);
                types[p1] = maxtype;
                types[p2] = maxtype;
            }
        }

        return current_max_rate.max(*collision_rate);
    }
    fn update_postions(&mut self, delta_t: f64, plate_x: f64, plate_dy: f64, plate_dz: f64) {
        let Particles {
            velocities,
            positions,
            types: _,
            parent_cells: _,
            keepers: _,
        } = self;
        let particle_count = velocities.len();

        let position_ptr = VectorPointer(positions.as_mut_ptr());
        let velocity_ptr = VectorPointer(velocities.as_mut_ptr());
        let particle_type_ptr = ParticleTypePointer(self.types.as_mut_ptr());

        (0..particle_count)
            .into_par_iter()
            .chunks(4)
            .for_each(|chunk| {
                chunk.iter().for_each(|i| {
                    let current_position = unsafe { &mut *{ position_ptr }.0.add(*i) };
                    let current_velocity = unsafe { &mut *{ velocity_ptr }.0.add(*i) };
                    let mut new_position = (*current_position) + (*current_velocity) * delta_t;

                    if (current_position.x < plate_x) && (new_position.x > plate_x)
                        || (current_position.x > plate_x) && (new_position.x < plate_x)
                    {
                        let t = (current_position[0] - plate_x)
                            / (current_position[0] - new_position[0]);
                        //create a particle representing the potential position
                        let tmp = *current_position * (1. - t)
                            + ((*current_position + delta_t * (*current_velocity)) * t);
                        if (-plate_dy..plate_dy).contains(&tmp[1])
                            && (-plate_dz..plate_dz).contains(&tmp[2])
                        {
                            new_position.x -= 2. * new_position.x - plate_x;
                            current_velocity.x = -current_velocity.x;
                            let ptype = unsafe { &mut *{ particle_type_ptr }.0.add(*i) };
                            *ptype = ParticleType::Thud;
                        }
                    }

                    if (-1.0..1.).contains(&new_position.y) {
                        new_position.y -= 2.0;
                    }
                    if (-1.0..1.).contains(&new_position.z) {
                        new_position.z -= 2.0;
                    }
                    *current_position = new_position;
                })
            });
    }

    /// any particles outside of the bounds of the grid need to be discarded as they cannot be indexed.
    /// Since this can happen only at the x=-1 or x=1 boundaries, we only need to check the x coordinate
    pub(crate) fn filter_particles(&mut self, num_cells: usize) {
        // self.vec.sort_by_key(|p| (p.parent_cell));
        // if let Some(cutoff) = particles.iter().rposition(|p| p.parent_cell < num_cells) {
        //     particles.truncate(cutoff + 1);
        // }
        let Particles {
            velocities,
            positions,
            types,
            parent_cells,
            keepers,
        } = self;
        let cell_ptr = ParentCellPointer(parent_cells.as_mut_ptr());
        //let mut keep=keepers.iter();
        rayon::join(
            || {
                rayon::join(
                    || {
                        let mut keep = keepers.iter();
                        velocities.retain(|_| *keep.next().unwrap());
                    },
                    || {
                        let mut keep = keepers.iter();
                        positions.retain(|_| *keep.next().unwrap());
                    },
                )
            },
            || {
                rayon::join(
                    || {
                        let mut keep = keepers.iter();
                        parent_cells.retain(|_| *keep.next().unwrap());
                    },
                    || {
                        let mut keep = keepers.iter();
                        types.retain(|_| *keep.next().unwrap());
                    },
                )
            },
        );
        keepers.clear();
        keepers.resize(velocities.len(), true);
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
    fn index(
        &mut self,
        sample_sender: SampleUpdateSender,
        num_x: usize,
        num_y: usize,
        num_z: usize,
    ) {
        let num_cells = num_x * num_y * num_z;
        //assuming number of cells must be even in the x direction
        let half_x = num_x / 2;

        let dy: f64 = 2. / num_y as f64;
        let dz: f64 = 2. / num_y as f64;

        let Particles {
            velocities,
            positions,
            types: _,
            parent_cells,
            keepers,
        } = self;

        let particle_count = velocities.len();
        let cell_ptr = ParentCellPointer(parent_cells.as_mut_ptr());
        let keep_ptr = BoolPointer(keepers.as_mut_ptr());
        (0..particle_count)
            .into_par_iter()
            .chunks(4)
            .map(|chunk| {
                chunk //iterate syncronously for chunk size
                    .iter()
                    .map(|i| {
                        let y_offset =
                            ((positions[*i].y + 1.0 / dy).floor() as usize).min(num_y - 1) * 2;
                        let z_offset = ((positions[*i].z + 1.0 * dz).floor() as usize)
                            .min(num_z - 1)
                            * (num_y * 2);
                        //figure out where a particle belongs based of it's location along the x-axis
                        let cell_membership: usize = match positions[*i] {
                            //scenario one: it's left of 0 or 0
                            p if (-1.0..=0.0).contains(&p.x) => {
                                let x_offset =
                                    (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
                                x_offset + y_offset + z_offset
                            }
                            //scenario two: it's right of zero or one
                            p if (0.0..=1.).contains(&p.x) => {
                                let x_offset =
                                    (((p.x).abs() * half_x as f64).floor() as usize).min(half_x);
                                x_offset + y_offset + z_offset + 1
                            }
                            //scenario three: it's no longer relevant
                            _ => num_cells,
                        };
                        let parent = unsafe { &mut *{ cell_ptr }.0.add(*i) };
                        *parent = cell_membership;
                        if cell_membership == num_cells {
                            let keep = unsafe { &mut *{ keep_ptr }.0.add(*i) };
                            *keep = false;
                        }
                        (cell_membership, *i)
                    })
                    .collect::<Vec<(usize, usize)>>()
            }) //send each chunk of indices to sample update function
            .for_each(|chunk| sample_sender.send(chunk).unwrap());
    }

    pub fn get_velocity(&self, index: usize) -> &Vector3<f64> {
        &(self.velocities[index])
    }

    pub fn set_velocity(&mut self, index: usize, velocity: Vector3<f64>) {
        self.velocities[index] = velocity;
    }
    pub fn get_position(&self, index: usize) -> &Vector3<f64> {
        &(self.positions[index])
    }

    pub fn set_position(&mut self, index: usize, position: Vector3<f64>) {
        self.positions[index] = position;
    }

    pub fn len(self) -> usize {
        self.parent_cells.len()
    }

    pub fn reserve(&mut self, capacity: usize) {
        let Particles {
            velocities,
            positions,
            types,
            parent_cells,
            keepers,
        } = self;
        rayon::join(
            || {
                rayon::join(
                    || {
                        rayon::join(
                            || velocities.reserve(capacity),
                            || positions.reserve(capacity),
                        )
                    },
                    || {
                        rayon::join(
                            || types.reserve(capacity),
                            || parent_cells.reserve(capacity),
                        )
                    },
                )
            },
            || keepers.reserve(capacity),
        );
        self.velocities.reserve(capacity);
        self.positions.reserve(capacity);
        self.parent_cells.reserve(capacity);
        self.types.reserve(capacity);
    }

    pub fn par_extend(&mut self, particles: Particles) {
        if self.velocities.capacity() - self.velocities.len() < particles.velocities.len() {
            self.reserve(particles.velocities.len())
        }
        let Particles {
            velocities,
            positions,
            types,
            parent_cells,
            keepers,
        } = self;

        let Particles {
            velocities: other_velocities,
            positions: other_positions,
            types: other_types,
            parent_cells: other_parent_cells,
            keepers: other_keepers,
        } = particles;

        rayon::join(
            || {
                rayon::join(
                    || {
                        rayon::join(
                            || velocities.par_extend(other_velocities),
                            || positions.par_extend(other_positions),
                        )
                    },
                    || {
                        rayon::join(
                            || types.par_extend(other_types),
                            || parent_cells.par_extend(other_parent_cells),
                        )
                    },
                )
            },
            || keepers.par_extend(other_keepers),
        );
    }
}

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
