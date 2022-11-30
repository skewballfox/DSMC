use crossbeam::channel::{Receiver, Sender};
use na::Vector3;

//types for sending the indices to the particle struct
pub(crate) type CollisionIndexSender = Sender<(Vec<(usize, usize)>, f64)>;
pub(crate) type CollisionIndexReceiver = Receiver<(Vec<(usize, usize)>, f64)>;

pub(crate) type SampleUpdateSender = Sender<(usize, usize, Vector3<f64>)>;
pub(crate) type SampleUpdateReceiver = Receiver<(usize, usize, Vector3<f64>)>;
