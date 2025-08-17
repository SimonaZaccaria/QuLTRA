## QuLTRA:  Quantum hybrid Lumped and TRansmission lines circuits Analyzer

QuLTRA is an open-source Python library for the analysis of superconducting circuits that include both lumped elements (capacitors, inductors, Josephson junctions, resistors) and distributed components (such as coplanar waveguides and CPW couplers).

The library provides tools to compute:

-Eigenfrequencies of the circuit modes,

-Dissipation rates, including coupling to external lines,

-The cross-Kerr matrix, describing the nonlinear interactions between modes.

QuLTRA is designed to support fast, layout-independent modeling in the early design stages, before full electromagnetic simulations.

## Installation

You can install QuLTRA directly from PyPI using pip:

```bash
pip install qultra
```

For more details, visit the PyPI project page: [https://pypi.org/project/qultra/](https://pypi.org/project/qultra/)

You can also install QuLTRA by directly cloning the repository

```bash
git clone https://github.com/SimonaZaccaria/QuLTRA.git
cd QuLTRA
pip install -e .
```
## Documentation
For the complete documentation, please visit the [QuLTRA web page](https://simonazaccaria.github.io/QuLTRA/).

## Future developments
Future improvements of QuLTRA may involve the implementation of a graphical interface, 
the possibility to export the full Hamiltonian of quantum circuits as a Qobj() 
for seamless integration with QuTiP, the calculation of zero-point 
fluctuation elements via energy participation coefficients that account for higher-order 
terms in the Taylor expansion, and the extension of the method to circuits with strong 
anharmonicity.

## Contact
For inquiries or to contribute to the project, please contact Simona at simona.zaccaria4@unibo.it.