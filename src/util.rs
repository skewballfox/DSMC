use std::{
    f64::consts::PI,
    ops::DerefMut,
    sync::{
        atomic::{AtomicBool, AtomicU64},
        Arc,
    },
};

use crossbeam::{
    channel::{Receiver, Sender},
    sync::Parker,
};
use na::Vector3;
use rand::prelude::*;
use rayon::prelude::*;

//types for sending the indices to the particle struct
pub(crate) type CollisionIndexSender = Sender<Vec<(usize, usize)>>;
pub(crate) type CollisionIndexReceiver = Receiver<Vec<(usize, usize)>>;

pub(crate) type SampleUpdateSender = Sender<Vec<(usize, usize)>>;
pub(crate) type SampleUpdateReceiver = Receiver<Vec<(usize, usize)>>;

/// Computes a unit vector with a random orientation and uniform distribution
pub(crate) fn random_direction(rng: &mut ThreadRng) -> na::Vector3<f64> {
    let b = 2. * rng.gen::<f64>() - 1.;
    let a: f64 = (1. - b * b).sqrt();
    let theta = rand::random::<f64>() * 2. * PI;
    na::Vector3::<f64>::new(b, a * theta.cos(), a * theta.sin())
}
/// Computes a velocity magnitude distribution that would correspond to
/// thermal equilibrium
pub(crate) fn random_velocity(rng: &mut ThreadRng, region_temp: f64) -> f64 {
    return region_temp * (-f64::max(rng.gen::<f64>(), 1e-200).ln()).sqrt();
}

//see https://stackoverflow.com/a/70848420/11019565
//thin wrapper over pointer to make it Send/Sync
#[derive(Copy, Clone)]
pub(crate) struct ParticlePointer(pub(crate) *mut Particle);
unsafe impl Send for ParticlePointer {}
unsafe impl Sync for ParticlePointer {}

#[derive(Copy, Clone)]
pub(crate) struct MaxPointer(pub(crate) *mut AtomicU64);
unsafe impl Send for MaxPointer {}
unsafe impl Sync for MaxPointer {}

#[derive(Copy, Clone)]
pub(crate) struct CellPointer(pub(crate) *mut CellSample);
unsafe impl Send for CellPointer {}
unsafe impl Sync for CellPointer {}

#[derive(Copy, Clone)]
pub(crate) struct ParentCellPointer(pub(crate) *mut usize);
unsafe impl Send for ParentCellPointer {}
unsafe impl Sync for ParentCellPointer {}

#[derive(Copy, Clone)]
pub(crate) struct CellMeanPointer(pub(crate) *mut f64);
unsafe impl Send for CellMeanPointer {}
unsafe impl Sync for CellMeanPointer {}

#[derive(Copy, Clone)]
pub(crate) struct ImmutableParentCellPointer(pub(crate) *const usize);
unsafe impl Send for ImmutableParentCellPointer {}
unsafe impl Sync for ImmutableParentCellPointer {}

#[derive(Copy, Clone)]
pub(crate) struct VectorPointer(pub(crate) *mut Vector3<f64>);
unsafe impl Send for VectorPointer {}
unsafe impl Sync for VectorPointer {}
#[derive(Copy, Clone)]
pub(crate) struct ImmutableVectorPointer(pub(crate) *const Vector3<f64>);
unsafe impl Send for ImmutableVectorPointer {}
unsafe impl Sync for ImmutableVectorPointer {}
#[derive(Copy, Clone)]
pub(crate) struct ParticleTypePointer(pub(crate) *mut ParticleType);
unsafe impl Send for ParticleTypePointer {}
unsafe impl Sync for ParticleTypePointer {}

#[derive(Copy, Clone)]
pub(crate) struct BoolPointer(pub(crate) *mut bool);
unsafe impl Send for BoolPointer {}
unsafe impl Sync for BoolPointer {}

#[derive(Copy, Clone)]
pub(crate) struct PlockPointer(pub(crate) *const Arc<AtomicBool>);
unsafe impl Send for PlockPointer {}
unsafe impl Sync for PlockPointer {}
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

pub(crate) fn gen_outer_range(
    num_x: usize,
    num_y: usize,
    num_z: usize,
    mppc: usize,
) -> Vec<(f64, f64, f64)> {
    let dx = 2. / (num_x as f64);
    let dy = 2. / (num_y as f64);
    let dz = 2. / (num_z as f64);
    //let mut outer_range = Vec::with_capacity(num_y * num_z * mppc);

    (0..num_y)
        .into_par_iter()
        .map(|i| -> Vec<(f64, f64, f64)> {
            (0..num_z)
                .into_par_iter()
                .map(|j| -> Vec<(f64, f64, f64)> {
                    (0..mppc)
                        .into_par_iter()
                        .map(|k| -> (f64, f64, f64) {
                            let current_x = (-1. - dx) * dx;
                            let current_y = (-1. + i as f64 * dy) * dy;
                            let current_z = (-1. + j as f64 * dz) * dz;
                            (current_x, current_y, current_z)
                        })
                        .collect()
                })
                .flatten()
                .collect()
        })
        .flatten()
        .collect::<Vec<(f64, f64, f64)>>()
}

///Ironically enough, I defined this to make my code more legible, takes a pointer to the particle array,
/// pointer to max array, and gets the first index out of the collision indices to figure out the parent cell's index
/// then uses that to get the current max
macro_rules! get_max {
    ($coll_idx:ident, $particle_arr:ident,$curr_max_p: ident) => {
        unsafe {
            &*{ $curr_max_p }
                .0
                .add(unsafe { &mut *{ $particle_arr }.0.add($coll_idx[0].0) }.parent_cell)
        }
    };
}
pub(crate) use get_max;

use crate::{
    particles::{Particle, ParticleType},
    CellSample,
};

pub(crate) fn create_parking_lot(lot_size: usize) -> Vec<Arc<AtomicBool>> {
    let mut parking_lot: Vec<Arc<AtomicBool>> = Vec::with_capacity(lot_size);
    for i in 0..lot_size {
        parking_lot.push(Arc::new(AtomicBool::new(true)))
    }
    return parking_lot;
}

pub(crate) fn grow_parking_lot(parking_lot: &mut Vec<Arc<AtomicBool>>, lot_size: usize) {
    for i in parking_lot.len() - 1..lot_size {
        parking_lot.push(Arc::new(AtomicBool::new(true)))
    }
}

//see https://github.com/rayon-rs/rayon/issues/721#issuecomment-585842377
macro_rules! extend_tuple {
    ( $name:ident, $( ( $fields:tt $types:ident ) ),+ ) => {
        pub struct $name<'a, $($types),+> {
            tuple: ($(&'a mut Vec<$types>),+),
        }

        impl<'a, $($types),+> $name<'a, $($types),+> {
            pub fn new(tuple: ($(&'a mut Vec<$types>),+)) -> Self {
                Self { tuple }
            }
        }

        impl<$($types),+> ::rayon::iter::ParallelExtend<($($types),+)> for $name<'_, $($types),+>
        where
            $(
                $types: Send,
            )+
        {
            fn par_extend<PI>(&mut self, par_iter: PI)
            where
                PI: ::rayon::iter::IntoParallelIterator<Item = ($($types),+)>,
            {
                use ::std::{
                    collections::LinkedList, ptr, slice, sync::atomic::{AtomicUsize, Ordering},
                };

                use ::rayon::{
                    iter::plumbing::{Consumer, Folder, Reducer, UnindexedConsumer},
                    prelude::*,
                };

                struct NoopReducer;

                impl Reducer<()> for NoopReducer {
                    fn reduce(self, _left: (), _right: ()) {}
                }

                struct CollectTupleConsumer<'c, $($types: Send),+> {
                    writes: &'c AtomicUsize,
                    targets: ($(&'c mut [$types]),+),
                }

                struct CollectTupleFolder<'c, $($types: Send),+> {
                    global_writes: &'c AtomicUsize,
                    local_writes: usize,
                    targets: ($(slice::IterMut<'c, $types>),+),
                }

                impl<'c, $($types: Send + 'c),+> Consumer<($($types),+)>
                for CollectTupleConsumer<'c, $($types),+>
                {
                    type Folder = CollectTupleFolder<'c, $($types),+>;
                    type Reducer = NoopReducer;
                    type Result = ();

                    fn split_at(self, index: usize) -> (Self, Self, NoopReducer) {
                        let CollectTupleConsumer { writes, targets } = self;

                        let splits = (
                            $(
                                targets.$fields.split_at_mut(index),
                            )+
                        );

                        (
                            CollectTupleConsumer {
                                writes,
                                targets: (
                                    $(
                                        splits.$fields.0,
                                    )+
                                ),
                            },
                            CollectTupleConsumer {
                                writes,
                                targets: (
                                    $(
                                        splits.$fields.1,
                                    )+
                                ),
                            },
                            NoopReducer,
                        )
                    }

                    fn into_folder(self) -> CollectTupleFolder<'c, $($types),+> {
                        CollectTupleFolder {
                            global_writes: self.writes,
                            local_writes: 0,
                            targets: (
                                $(
                                    self.targets.$fields.iter_mut(),
                                )+
                            ),
                        }
                    }

                    fn full(&self) -> bool {
                        false
                    }
                }

                impl<'c, $($types: Send + 'c),+> Folder<($($types),+)>
                for CollectTupleFolder<'c, $($types),+>
                {
                     type Result = ();

                    fn consume(
                        mut self,
                        item: ($($types),+),
                    ) -> CollectTupleFolder<'c, $($types),+> {
                        $(
                            let head = self
                                .targets
                                .$fields
                                .next()
                                .expect("too many values pushed to consumer");
                            unsafe {
                                ptr::write(head, item.$fields);
                            }
                        )+

                        self.local_writes += 1;
                        self
                    }

                    fn complete(self) {
                        self.global_writes.fetch_add(self.local_writes, Ordering::Relaxed);
                    }

                    fn full(&self) -> bool {
                        false
                    }
                }

                impl<'c, $($types: Send + 'c),+> UnindexedConsumer<($($types),+)>
                for CollectTupleConsumer<'c, $($types),+>
                {
                     fn split_off_left(&self) -> Self {
                        unreachable!("CollectTupleConsumer must be indexed!")
                    }
                    fn to_reducer(&self) -> Self::Reducer {
                        NoopReducer
                    }
                }

                struct CollectTuple<'c, $($types: Send),+> {
                    writes: AtomicUsize,
                    tuple: ($(&'c mut Vec<$types>),+),
                    len: usize,
                }

                impl<'c, $($types: Send),+> CollectTuple<'c, $($types),+> {
                    pub fn new(tuple: ($(&'c mut Vec<$types>),+), len: usize) -> Self {
                        Self {
                            writes: AtomicUsize::new(0),
                            tuple,
                            len,
                        }
                    }

                    pub fn as_consumer(&mut self) -> CollectTupleConsumer<'_, $($types),+> {
                        $(
                            self.tuple.$fields.reserve(self.len);
                        )+

                        CollectTupleConsumer {
                            writes: &self.writes,
                            targets: (
                                $(
                                    {
                                        let vec = &mut self.tuple.$fields;
                                        let start = vec.len();
                                        let slice = &mut vec[start..];
                                        unsafe {
                                            slice::from_raw_parts_mut(
                                                slice.as_mut_ptr(),
                                                self.len,
                                            )
                                        }
                                    }
                                ),+
                            ),
                        }
                    }

                    pub fn complete(mut self) {
                        unsafe {
                            let actual_writes = self.writes.load(Ordering::Relaxed);
                            assert!(
                                actual_writes == self.len,
                                "expected {} total writes, but got {}",
                                self.len,
                                actual_writes
                            );

                            $(
                                let vec = &mut self.tuple.$fields;
                                let new_len = vec.len() + self.len;
                                vec.set_len(new_len);
                            )+
                        }
                    }
                }

                let par_iter = par_iter.into_par_iter();
                match par_iter.opt_len() {
                    Some(len) => {
                        let mut collect = CollectTuple::new(($(self.tuple.$fields),+), len);
                        par_iter.drive_unindexed(collect.as_consumer());
                        collect.complete()
                    }
                    None => {
                        let list = par_iter
                            .into_par_iter()
                            .fold(|| ($(Vec::<$types>::new()),+), |mut vecs, elem| {
                                $(
                                    vecs.$fields.push(elem.$fields);
                                )+
                                vecs
                            })
                            .map(|item| {
                                let mut list = LinkedList::new();
                                list.push_back(item);
                                list
                            })
                            .reduce(LinkedList::new, |mut list1, mut list2| {
                                list1.append(&mut list2);
                                list1
                            });
                        let len = list.iter().map(|vecs| vecs.0.len()).sum();

                        $(
                            self.tuple.$fields.reserve(len);
                        )+
                        for mut vecs in list {
                            $(
                                self.tuple.$fields.append(&mut vecs.$fields);
                            )+
                        }
                    }
                }
            }
        }
    };
}

extend_tuple!(ExtendTuple2, (0 A), (1 B));
extend_tuple!(ExtendTuple3, (0 A), (1 B), (2 C));
extend_tuple!(ExtendTuple9, (0 A), (1 B), (2 C), (3 D), (4 E), (5 F), (6 G), (7 H), (8 I));

#[cfg(test)]
mod tests {
    use super::*;

    use rayon::prelude::*;

    #[test]
    fn tuple2() {
        let mut vec0 = vec![];
        let mut vec1 = vec![];

        ExtendTuple2::new((&mut vec0, &mut vec1))
            .par_extend((0..3).into_par_iter().map(|i| (i, i * 2)));

        assert_eq!(vec0, [0, 1, 2]);
        assert_eq!(vec1, [0, 2, 4]);
    }

    #[test]
    fn tuple9() {
        let mut vec0 = vec![];
        let mut vec1 = vec![];
        let mut vec2 = vec![];
        let mut vec3 = vec![];
        let mut vec4 = vec![];
        let mut vec5 = vec![];
        let mut vec6 = vec![];
        let mut vec7 = vec![];
        let mut vec8 = vec![];

        ExtendTuple9::new((
            &mut vec0, &mut vec1, &mut vec2, &mut vec3, &mut vec4, &mut vec5, &mut vec6, &mut vec7,
            &mut vec8,
        ))
        .par_extend(
            (0..3)
                .into_par_iter()
                .map(|i| (i, i + 1, i + 2, i + 3, i + 4, i + 5, i + 6, i + 7, i + 8)),
        );

        assert_eq!(vec0, [0, 1, 2]);
        assert_eq!(vec1, [1, 2, 3]);
        assert_eq!(vec2, [2, 3, 4]);
        assert_eq!(vec3, [3, 4, 5]);
        assert_eq!(vec4, [4, 5, 6]);
        assert_eq!(vec5, [5, 6, 7]);
        assert_eq!(vec6, [6, 7, 8]);
        assert_eq!(vec7, [7, 8, 9]);
        assert_eq!(vec8, [8, 9, 10]);
    }
}
