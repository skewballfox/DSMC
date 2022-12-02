use std::f64::consts::PI;

use crossbeam::channel::{Receiver, Sender};
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

pub(crate) fn gen_outer_range(num_y: usize, num_z: usize) -> Vec<(usize, usize)> {
    let mut outer_range = Vec::with_capacity(num_y * num_z);
    outer_range.par_extend(
        (0..num_y)
            .into_par_iter()
            .flat_map(|i| -> Vec<(usize, usize)> {
                let mut tmp = vec![(0, 0); num_z];
                (0..num_z)
                    .into_par_iter()
                    .map(|j| -> (usize, usize) { (i, j) })
                    .collect_into_vec(&mut tmp);
                tmp
            }),
    );
    outer_range
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
