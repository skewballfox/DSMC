use crossbeam::channel::{Receiver, Sender};
use na::Vector3;

//types for sending the indices to the particle struct
pub(crate) type CollisionIndexSender = Sender<Vec<(usize, usize)>>;
pub(crate) type CollisionIndexReceiver = Receiver<Vec<(usize, usize)>>;

pub(crate) type SampleUpdateSender = Sender<Option<(usize, usize, Vector3<f64>)>>;
pub(crate) type SampleUpdateReceiver = Receiver<Option<(usize, usize, Vector3<f64>)>>;
