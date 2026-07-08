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
