[![Made at CERN!](https://img.shields.io/badge/CERN-CERN%20openlab-blue)](https://openlab.cern/) 

Read the Docs Documentation: *soon available*

---

# Counterdiabatic Optimized Local Driving annealER

This package provides tool for Spin system simulation and optimization through a combination of Quantum Optimal Control and Counterdiabatic Driving techniques, as inspired by this [reference paper](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.010312).
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

## Getting started

The Documentation of this package provides all the required information to use the **physics interface** and run **annealing simulation & optimization**.

A set of examples is provided in the [examples](examples) folder.


## Bibliography

- **Counterdiabatic Optimized Local Driving** - Ieva and Polkovnikov, Anatoli and Daley, Andrew J. and Duncan, Callum W. ([PRXQuantum 4 010312](https://doi.org/10.1103/PRXQuantum.4.010312))
