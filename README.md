[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16897531.svg)](https://doi.org/10.5281/zenodo.16897531)
[![PyPI version](https://img.shields.io/pypi/v/qultra.svg)](https://pypi.org/project/qultra/)
![License](https://img.shields.io/github/license/SimonaZaccaria/QuLTRA)
[![Docs](https://img.shields.io/badge/docs-online-blue.svg)](https://simonazaccaria.github.io/QuLTRA/)
[![Build Status](https://github.com/SimonaZaccaria/QuLTRA/actions/workflows/tests.yml/badge.svg)](https://github.com/SimonaZaccaria/QuLTRA/actions)





## QuLTRA:  Quantum hybrid Lumped and TRansmission lines circuits Analyzer

QuLTRA is an open-source Python library for the analysis of superconducting circuits that include both lumped elements (capacitors, inductors, Josephson junctions, resistors) and distributed components (such as coplanar waveguides and CPW couplers).

The library provides tools to compute:

-Eigenfrequencies of the circuit modes,

-Dissipation rates, including coupling to external lines,

-The cross-Kerr matrix, describing the nonlinear interactions between modes.

QuLTRA is designed to support fast, layout-independent modeling in the early design stages, before full electromagnetic simulations.

## References

 QuLTRA Python package, release v1.0.1, archived on Zenodo: DOI [10.5281/zenodo.16897531](https://doi.org/10.5281/zenodo.16897531). The preprint of the article is on arXiv:[https://arxiv.org/abs/2509.03651](https://arxiv.org/abs/2509.03651)
 

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
and the extension of the method to circuits with strong 
anharmonicity.

## Contact
For inquiries or to contribute to the project, please contact Simona at simona.zaccaria4@unibo.it.