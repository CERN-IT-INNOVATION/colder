[![Made at CERN!](https://img.shields.io/badge/CERN-CERN%20openlab-blue)](https://openlab.cern/) 
[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1367--2630%2Fad313e-blue)](https://doi.org/10.1088/1367-2630/ad313e)

Read the Docs Documentation available [here](https://colder.readthedocs.io/en/latest/index.html).

---

# Counterdiabatic Optimized Local Driving annealER

This package provides tool for spin system simulation and optimization through a combination of Quantum Optimal Control and Counterdiabatic Driving techniques, as in this [reference paper](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.010312).
The framework allows to prompt the system through a simplified interface and automatically solve the minimization of the Adiabatic Gauge Potential (AGP) for a given local ansatz.

*Package creator*: Barone Francesco Pio, openlab summer 2023 intern

*Supervisors*: Oriel Kiss, Antonio Mandarino, Michele Grossi, Sofia Vallecorsa


## Installation

To install this package, clone this repo and install through pip:
```bash
git clone https://github.com/CERN-IT-INNOVATION/colder
cd colder
pip install .
```

### structure of the package

The colder package consists of the following modules:
- **backend**: implements the time evolution routines in different backends
- **core**: provides the I/O of spin model and the algebraic routines used for the symbolical optimization
- **gauge**: class that optimizes symbolically the AGP ansatz
- **models**: provides useful functions to specify the spin model interactions
- **quantum**: implements basic functions to manipulate quantum states, regardless of the backend
- **simulation**: implements the time evolution interface and the COLD object


## Getting started

The Documentation of this package provides all the required information to use the **physics interface** and run **annealing simulation & optimization**.

A set of examples is provided in the [examples](examples) folder.

## How to cite

If you used this package for your research, please cite:
```text
@article{Barone_2024,
    doi = {10.1088/1367-2630/ad313e},
    url = {https://dx.doi.org/10.1088/1367-2630/ad313e},
    year = {2024},
    month = {mar},
    publisher = {IOP Publishing},
    volume = {26},
    number = {3},
    pages = {033031},
    author = {Francesco Pio Barone and Oriel Kiss and Michele Grossi and Sofia Vallecorsa and Antonio Mandarino},
    title = {Counterdiabatic optimized driving in quantum phase sensitive models},
    journal = {New Journal of Physics}
}
```

## Bibliography

- **Counterdiabatic Optimized Local Driving** - Ieva and Polkovnikov, Anatoli and Daley, Andrew J. and Duncan, Callum W. ([PRXQuantum 4 010312](https://doi.org/10.1103/PRXQuantum.4.010312))
