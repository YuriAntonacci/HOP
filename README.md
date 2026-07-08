# HOP – High-Order Analysis of Random Processes

**HOP** is a MATLAB toolbox for the **High-Order analysis of random Processes**.

The toolbox provides methods for the analysis of high-order interactions in multivariate dynamical systems, with a specific focus on the decomposition of directed and causal interactions in networks of oscillatory processes.

HOP includes tools for the assessment of **unique**, **redundant**, and **synergistic** contributions carried by multiple source processes toward a target process, extending classical pairwise descriptions of information transfer and Granger-causal interactions.

---

## Overview

Granger causality is a widely used statistical framework for the inference of directional influences between time series recorded from complex networks. However, in multivariate systems, the interaction between a set of source processes and a target process may involve high-order effects that cannot be fully described by pairwise interactions alone.

In this context, HOP implements methods designed to dissect multivariate interactions into interpretable high-order components, allowing the user to quantify whether the influence of multiple sources on a target is:

- **unique**, when it is carried exclusively by one source;
- **redundant**, when it is shared by two or more sources;
- **synergistic**, when it emerges only from the joint observation of multiple sources.

These measures are particularly useful for the analysis of complex physiological, neural, and oscillatory systems, where collective dynamics may arise from non-trivial interactions among multiple subsystems.

---

## Main methods

The toolbox includes implementations of the following approaches.

### Partial Decomposition of Granger Causality

**Partial Decomposition of Granger Causality** allows the decomposition of multivariate Granger-causal effects from a set of driver processes to a target process into unique, redundant, and synergistic components.

This approach combines the framework of Granger causality with the principles of partial information decomposition, providing a way to characterize how different source processes contribute to the predictability of a target process.

### Partial Information Rate Decomposition

**Partial Information Rate Decomposition** allows the decomposition of the information rate shared between multiple source processes and a target process into high-order informational components.

This framework provides a complementary information-theoretic description of high-order dependencies in multivariate stochastic processes.

## Citate our work
If you use this code, please cite:

```bibtex
@article{faes2026dissecting,
  title={Dissecting Spectral Granger Causality through Partial Information Decomposition},
  author={Faes, Luca and Mijatovic, Gorana and Pernice, Riccardo and Marinazzo, Daniele and Stramaglia, Sebastiano and Antonacci, Yuri},
  journal={arXiv preprint arXiv:2603.07634},
  year={2026}
}
```

## Contact

For questions, please contact: yuri.antonacci@unipa.it

**Yuri Antonacci**
Department of Engineering, University of Palermo, Palermo, Italy
GitHub: [@YuriAntonacci](https://github.com/YuriAntonacci)


## License

Copyright (c) 2026 Yuri Antonacci

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files, to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

