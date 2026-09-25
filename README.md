# Moment--SOS Risk Quantification for Multistage Colorectal Carcinogenesis

This repository contains the computational implementation and supplementary
materials associated with the manuscript:

**"Moment--SOS Risk Quantification for Multistage Colorectal Carcinogenesis:
A Semidefinite Optimization Framework for Stochastic Cancer Progression"**

The repository provides the MATLAB implementation of the proposed
Moment--SOS computational framework, numerical experiments, validation
procedures, and supplementary parameter tables required to reproduce the
reported results.

---

## 1. Overview

This work develops a mathematical framework for quantifying stochastic
carcinogenesis risk using a continuous-time Markov jump process (CTMC) combined
with Moment--SOS semidefinite relaxations.

The computational pipeline includes:

- Construction of the multistage APC--KRAS--TP53 carcinogenesis model
- Definition of stochastic transition events and propensity functions
- Generator-based moment formulation
- Construction of Moment--SOS relaxations
- Semidefinite optimization for lower and upper risk bounds
- First-passage risk analysis
- Gillespie stochastic simulation benchmarks
- Sensitivity and convergence analyses

The framework provides deterministic uncertainty bounds for stochastic cancer
progression models without requiring repeated sampling for risk estimation.

---

# 2. Repository Structure
