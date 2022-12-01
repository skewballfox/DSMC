use std::f64::consts::PI;

use crossbeam::channel::{Receiver, Sender};
use na::Vector3;
use rand::prelude::*;

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
