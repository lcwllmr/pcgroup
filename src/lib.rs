pub mod util;
pub mod zoo;

mod word;
pub use crate::word::{Term, Word};

mod presentation;
pub use crate::presentation::Presentation;

mod builder;
pub use crate::builder::{Builder, BuilderError};

mod collector;
pub use crate::collector::Collector;

mod element;
pub use crate::element::Element;

mod consistency;
pub use crate::consistency::{ConsistencyError, verify_consistency};

mod generating_sequence;
pub use crate::generating_sequence::GeneratingSequence;

mod series;
pub use crate::series::{
    chief_series, commutator_subgroup, is_abelian, is_abelian_subgroup, is_nilpotent,
    is_nilpotent_subgroup, is_supersolvable, lower_central_series, nilpotency_class,
    nilpotency_class_subgroup,
};
