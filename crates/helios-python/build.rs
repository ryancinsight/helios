//! Configure standalone Cargo builds of the Python extension.

fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
