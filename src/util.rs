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
pub(crate) struct ImmutableParentCellPointer(pub(crate) *const usize);
unsafe impl Send for ImmutableParentCellPointer {}
unsafe impl Sync for ImmutableParentCellPointer {}

#[derive(Copy, Clone)]
pub(crate) struct VectorPointer(pub(crate) *mut Vector3<f64>);
unsafe impl Send for VectorPointer {}
unsafe impl Sync for VectorPointer {}

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
