
# Prerequisites

To install this package, the following is needed:
- `C++ compiler:` Intel >= 18.0
- `Fortran compiler`
- `System environment:` gcc>=6.3.0, python3, scipy, numpy, mpi4py
- `Library:` Intel MKL

# Complilation

## Building DFT codes
If you want to do DFT+DMFT calculations based on the ABACUS, you need to build the ABACUS at first (Refer to [Building ABACUS](#building-abacus)). If you want to do DFT+DMFT calculations based on the FHI-aims, you need to build the FHI-aims at first (Refer to [Building FHI-aims](#building-fhi-aims)). We suggest you to use the same compliers and libraries when building the DFT codes and this DFT+DMFT package.

### Building ABACUS

- `Download:`ABACUS is an open-source package, and you can download the latest version (a safe suggestion) from the [github](https://github.com/abacusmodeling/abacus-develop).

- `Compilation:` Compile the ABACUS following the introduction in its manual.


### Building FHI-aims

- `Download:`FHI-aims is a full-potential all-electrons LCNAO DFT package, and you can get the latest version (a safe suggestion) from its [offical website](https://fhi-aims.org/).