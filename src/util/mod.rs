mod linalg;
mod ntheory;

pub use self::linalg::{ModPMatrix, ModPSubspace, composition_series};
pub use self::ntheory::{factorize, is_prime, mod_inverse};
