mod cyclic;
pub use self::cyclic::cyclic;

mod abelian;
pub use self::abelian::abelian;

mod dihedral;
pub use self::dihedral::dihedral;

mod quaternion;
pub use self::quaternion::quaternion;

mod affine1d;
pub use self::affine1d::affine1d;

mod unipotent;
pub use self::unipotent::unipotent;

mod upper_triangular;
pub use self::upper_triangular::upper_triangular;

mod outliers;
pub use self::outliers::{a4, s3, s4, sl2_3};
