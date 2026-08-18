mod linalg;
pub use self::linalg::{ModularMatrix, ModularSubspace, composition_series};

mod ntheory;
pub use self::ntheory::{factorize, is_prime, mod_inverse};
