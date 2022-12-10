use std::{
    collections::HashMap,
    fmt::{Display, Formatter},
    hint,
    ops::{Deref, DerefMut},
    sync::{
        atomic::{AtomicBool, AtomicU64, Ordering},
        Arc, Mutex,
    },
};

use crate::{
    random_direction, random_velocity,
    util::{
        create_parking_lot, gen_outer_range, grow_parking_lot, BoolPointer, CollisionIndexReceiver,
        ExtendTuple2, ExtendTuple3, ImmutableParentCellPointer, MaxPointer, ParentCellPointer,
        ParticleTypePointer, PlockPointer, SampleUpdateSender, VectorPointer,
    },
    CellSample, SIGMA_K,
};
use crossbeam::{
    atomic::AtomicCell,
    channel::{unbounded, Receiver, Sender},
    sync::Parker,
};
use dashmap::DashMap;
use na::default_allocator;
use nalgebra::{Matrix, Vector3};
use once_cell::sync::OnceCell;
use rand::prelude::*;
use rayon::prelude::*;

#[derive(PartialOrd, Ord, PartialEq, Eq, Copy, Clone, Debug, Default)]
pub enum ParticleType {
    #[default]
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

#[derive(Debug)]
pub struct Particles {
    //receiver: ThingBuf,
    pub velocities: Vec<na::Vector3<f64>>,
    pub positions: Vec<na::Vector3<f64>>,
    pub types: Vec<ParticleType>,
}

impl Particles {
    pub fn new() -> Self {
        Self {
            velocities: Vec::new(),
            positions: Vec::new(),
            types: Vec::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            velocities: Vec::with_capacity(capacity),
            positions: Vec::with_capacity(capacity),
            types: Vec::new(),
        }
    }

    // pub fn collide(
    //     &mut self,
    //     collision_idx_receiver: &CollisionIndexReceiver,
    //     collision_max_rate: Arc<DashMap<usize, AtomicU64>>,
    // ) {
    //     let mut cell: OnceCell<Vec<Arc<AtomicBool>>> = OnceCell::new();
    //     {
    //         cell.get_or_init(|| create_parking_lot(self.len()));
    //     }
    //     let mut plock = cell.get_mut().unwrap();

    //     grow_parking_lot(plock, self.len());

    //     //let current_max_ptr = MaxPointer(collision_max_rate.as_mut_ptr());
    //     let cell_ptr = ParentCellPointer(self.parent_cells.as_mut_ptr());
    //     let velocity_ptr = VectorPointer(self.velocities.as_mut_ptr());
    //     let particle_type_ptr = ParticleTypePointer(self.types.as_mut_ptr());
    //     let plock_ptr = PlockPointer(unsafe { plock.as_ptr() });

    //     collision_idx_receiver
    //         .into_iter()
    //         .par_bridge()
    //         .for_each(|big_chunk| {
    //             //println!("got chunk: {:?}", big_chunk);

    //             //used first particle index to figure out what cell the particle belongs to
    //             let cell_idx = unsafe { &mut *{ cell_ptr }.0.add(big_chunk[0].0) };
    //             //use the cell index to get the current max
    //             //let current_max_atomic = unsafe { &*{ current_max_ptr }.0.add(*cell_idx) }; //collision_max_rate[*cell_idx].clone(); //
    //             let current_max = f64::from_bits(
    //                 collision_max_rate
    //                     .get(cell_idx)
    //                     .unwrap()
    //                     .load(Ordering::Relaxed),
    //             );

    //             big_chunk
    //                 .into_par_iter() //.chunks(8).for_each(|lil_chunk| {
    //                 //lil_chunk.into_iter()
    //                 .for_each(|(p1, p2)| {
    //                     let mut rng = rand::thread_rng();
    //                     //println!("got collision indices: {:?}", collision_indices);
    //                     let lock1 = unsafe { &*{ plock_ptr }.0.add(p1) }.clone();
    //                     let lock2 = unsafe { &*{ plock_ptr }.0.add(p2) }.clone();

    //                     while !lock1.load(Ordering::Acquire) || !lock2.load(Ordering::Acquire) {
    //                         println!("idling");
    //                         hint::spin_loop();
    //                     }
    //                     //println!("locking particles {:#?}{:#?}", lock1, lock2);
    //                     lock1.store(false, Ordering::Release);
    //                     lock2.store(false, Ordering::Release);
    //                     //println!("parked thread");
    //                     //get the pointers for setting the fields
    //                     let velocity1 = unsafe { &mut *{ velocity_ptr }.0.add(p1) };
    //                     let ptype1 = unsafe { &mut *{ particle_type_ptr }.0.add(p1) };
    //                     let velocity2 = unsafe { &mut *{ velocity_ptr }.0.add(p2) };
    //                     let ptype2 = unsafe { &mut *{ particle_type_ptr }.0.add(p2) };
    //                     //println!("got the data");
    //                     let relative_velocity = *velocity1 - *velocity2;
    //                     let relative_velocity_norm = relative_velocity.norm();
    //                     let collision_rate = (SIGMA_K * relative_velocity_norm);
    //                     let threshold = if current_max.lt(&collision_rate) {
    //                         1.
    //                     } else {
    //                         collision_rate / current_max
    //                     };
    //                     let r: f64 = rng.gen();
    //                     if r < threshold {
    //                         // Collision Accepted, adjust particle velocities
    //                         // Compute center of mass velocity, vcm
    //                         let center_of_mass = 0.5 * (*velocity1 + *velocity2);
    //                         // Compute random perturbation that conserves momentum
    //                         let perturbation = random_direction(&mut rng) * relative_velocity_norm;
    //                         //make sure current

    //                         // Adjust particle velocities to reflect collision
    //                         *velocity1 = center_of_mass + 0.5 * perturbation;
    //                         *velocity2 = center_of_mass - 0.5 * perturbation;

    //                         // Bookkeeping to track particle interactions
    //                         if *ptype1 != ParticleType::Inflow || *ptype1 != ParticleType::Inflow {
    //                             let maxtype = *ptype1.max(ptype2);
    //                             *ptype1 = maxtype;
    //                             *ptype2 = maxtype;
    //                         }
    //                     }
    //                     if current_max.lt(&collision_rate) {
    //                         collision_max_rate
    //                             .get(cell_idx)
    //                             .unwrap()
    //                             .store(collision_rate.to_bits(), Ordering::Relaxed);
    //                     }
    //                     lock1.store(true, Ordering::Release);
    //                     lock2.store(true, Ordering::Release);
    //                     //println!("released locks")
    //                     //});
    //                 })
    //         });

    //     return;
    // }

    // Initialize particles at inflow boundaries
    /// Create particles at inflow boundary
    /// This works by creating particles at a ghost cell just before the boundary
    /// any particles that don't make it into the domain are discarded.
    pub fn initialize_boundaries(
        &mut self,
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

        let cell: OnceCell<Vec<(f64, f64, f64)>> = OnceCell::new();
        let outer_range =
            cell.get_or_init(|| gen_outer_range(num_x, num_y, num_z, mean_particles_per_cell));
        let growth_amount = outer_range.len() * mean_particles_per_cell;
        self.reserve(growth_amount);
        let old_len = self.len();

        let Particles {
            velocities,
            positions,
            types,
        } = self;

        //writing directly to uninitialized memory, do not try this at home
        //let position_ptr = VectorPointer(positions.as_mut_ptr());
        //let velocity_ptr = VectorPointer(velocities.as_mut_ptr());
        // let particle_type_ptr = ParticleTypePointer(self.types.as_mut_ptr());
        // let parent_cell_ptr = ParentCellPointer(self.parent_cells.as_mut_ptr());
        // let keep_ptr = BoolPointer(self.keepers.as_mut_ptr());

        rayon::join(
            || {
                ExtendTuple2::new((velocities, positions)).par_extend(
                    outer_range.into_par_iter().map(|(x, y, z)| {
                        let mut rng = rand::thread_rng();
                        let position = na::Vector3::<f64>::new(
                            *x + rng.gen::<f64>(),
                            *y + rng.gen::<f64>(),
                            *z + rng.gen::<f64>(),
                        );
                        //println!("{}", position);
                        let mut velocity =
                            random_velocity(&mut rng, region_temp) * random_direction(&mut rng);
                        velocity.x += mean_velocity;
                        //start writing to array here

                        // }
                        (velocity, position)
                    }),
                );

                //velocities.par_extend(vel_buf.into_par_iter());
                //positions.par_extend(pos_buf.into_par_iter());
            },
            || {
                types.par_extend(
                    (0..growth_amount)
                        .into_par_iter()
                        .map(|_| ParticleType::default()),
                );
            },
        );

        //println!("particle count: {:?}", self.len());
    }
    pub fn update_positions(&mut self, delta_t: f64, plate_x: f64, plate_dy: f64, plate_dz: f64) {
        let Particles {
            velocities,
            positions,
            types,
        } = self;
        let particle_count = velocities.len();

        (velocities, positions, types)
            .into_par_iter()
            //.chunks(4)
            //.for_each(|chunk| {
            //    chunk
            //        .iter()
            .for_each(move |(current_velocity, current_position, ptype)| {
                let mut new_position = (*current_position) + (*current_velocity) * delta_t;

                if (current_position.x < plate_x) && (new_position.x > plate_x)
                    || (current_position.x > plate_x) && (new_position.x < plate_x)
                {
                    let t =
                        (current_position[0] - plate_x) / (current_position[0] - new_position[0]);
                    //create a particle representing the potential position
                    let tmp = *current_position * (1. - t)
                        + ((*current_position + delta_t * (*current_velocity)) * t);
                    if (-plate_dy..plate_dy).contains(&tmp[1])
                        && (-plate_dz..plate_dz).contains(&tmp[2])
                    {
                        new_position.x -= 2. * new_position.x - plate_x;
                        current_velocity.x = -current_velocity.x;
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
                //        })
            });
    }

    /// any particles outside of the bounds of the grid need to be discarded as they cannot be indexed.
    /// Since this can happen only at the x=-1 or x=1 boundaries, we only need to check the x coordinate
    pub(crate) fn filter_out_of_scope(&mut self, off_grid: &Vec<usize>) {
        // self.vec.sort_by_key(|p| (p.parent_cell));
        // if let Some(cutoff) = particles.iter().rposition(|p| p.parent_cell < num_cells) {
        //     particles.truncate(cutoff + 1);
        // }
        let Particles {
            velocities,
            positions,
            types,
        } = self;
        //off_grid;
        let mut keep: Vec<bool> = vec![true; velocities.len()];
        let kptr: BoolPointer = BoolPointer(keep.as_mut_ptr());
        off_grid.into_par_iter().for_each(move |i| {
            let k = unsafe { &mut *(kptr.clone()).0.add(*i) };
            *k = false;
            println!("k: {}", *k);
        });
        let mut tmp_vel = Vec::with_capacity(velocities.len());
        let mut tmp_pos = Vec::with_capacity(velocities.len());
        let mut tmp_typ = Vec::with_capacity(velocities.len());
        //let cell_ptr = ParentCellPointer(parent_cells.as_mut_ptr());
        println!("Starting filter");
        //println!("keep: {:?}", keep);
        ExtendTuple3::new((&mut tmp_vel, &mut tmp_pos, &mut tmp_typ)).par_extend(
            (velocities, positions, types, keep)
                .into_par_iter()
                .filter_map(
                    |(vel, pos, typ, keep)| {
                        if keep {
                            Some((vel, pos, typ))
                        } else {
                            None
                        }
                    },
                ),
        );
        Self {
            velocities: tmp_vel,
            positions: tmp_pos,
            types: tmp_typ,
        };
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
    pub fn index(
        &mut self,

        num_x: usize,
        num_y: usize,
        num_z: usize,
    ) -> Arc<HashMap<usize, Vec<usize>>> {
        let num_cells = num_x * num_y * num_z;
        //assuming number of cells must be even in the x direction
        let half_x = num_x / 2;

        let dy: f64 = 2. / num_y as f64;
        let dz: f64 = 2. / num_y as f64;

        let Particles {
            velocities,
            positions,
            types: _,
        } = self;

        let particle_count = velocities.len();

        let memberships: HashMap<usize, Vec<usize>> = positions
            .into_par_iter()
            .enumerate()
            //.chunks(4)
            .map(
                //|chunk| {
                //chunk //iterate syncronously for chunk size
                //.iter()
                //.map(
                |(i, position)| -> (usize, usize) {
                    let y_offset = ((position.y + 1.0 / dy).floor() as usize).min(num_y - 1) * 2;
                    let z_offset =
                        ((position.z + 1.0 * dz).floor() as usize).min(num_z - 1) * (num_y * 2);
                    //figure out where a particle belongs based of it's location along the x-axis
                    let cell_membership: usize = match position {
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

                    (cell_membership, i)
                },
            )
            .fold(
                || HashMap::new(),
                |mut hm, (cell_index, particle_index): (usize, usize)| {
                    hm.entry(cell_index).or_insert(vec![]).push(particle_index);
                    hm
                },
            )
            .reduce(
                || HashMap::new(),
                |m1, m2| {
                    m2.iter().fold(m1, |mut acc, (cell_index, members)| {
                        acc.entry(*cell_index).or_insert(vec![]).extend(members);
                        acc
                    })
                },
            );
        Arc::new(memberships)
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

    pub fn len(&self) -> usize {
        self.velocities.len()
    }

    pub fn reserve(&mut self, capacity: usize) {
        let Particles {
            velocities,
            positions,
            types,
        } = self;

        rayon::join(
            || {
                rayon::join(
                    || velocities.reserve(capacity),
                    || positions.reserve(capacity),
                )
            },
            || {
                types.reserve(capacity);
            },
        );
    }

    pub fn par_extend(&mut self, particles: Particles) {
        if self.velocities.capacity() - self.velocities.len() < particles.velocities.len() {
            self.reserve(particles.velocities.len())
        }
        let Particles {
            velocities,
            positions,
            types,
        } = self;

        let Particles {
            velocities: other_velocities,
            positions: other_positions,
            types: other_types,
        } = particles;

        rayon::join(
            || {
                rayon::join(
                    || velocities.par_extend(other_velocities),
                    || positions.par_extend(other_positions),
                )
            },
            || types.par_extend(other_types),
        );
    }
}

impl Display for Particles {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        (0..self.len())
            .into_iter()
            .try_for_each(|i| -> Result<(), std::fmt::Error> {
                write!(
                    f,
                    "{} {} {} {}\n",
                    self.positions[i].x,
                    self.positions[i].y,
                    self.positions[i].z,
                    self.types[i].to_num()
                )
            })
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
