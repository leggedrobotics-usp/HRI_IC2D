HRI_IC2D
Human–Robot Interaction Force Estimation with the IC2D Test Bench

This repository provides experimental data and MATLAB implementations of three estimators for human–robot interaction (HRI) forces, validated on the IC2D (Impedance Control in 2 Dimensions) platform. The goal is to evaluate methods for improving transparency in wearable robots when direct force sensing is unavailable or limited.

The three main approaches are:
- **Kalman Filter (KF) with automatic noise tuning via Particle Swarm Optimization**
- **Adaptive Extended Kalman Filter (EKF) with real-time noise covariance adaptation and smoothing**
- **Monte Carlo (MC) method for parameter estimation**

Sure! Here's a clean, clear rewrite of your repository structure section for the README:

---

## Repository Structure

```
HRI_IC2D/
│
├── data/                        # Raw data files from dSPACE MicroLabBox (MF4 format)
│   ├── spring2k_1.mf4           
│   ├── spring2k_2.mf4           
│   ├── spring5k_1.mf4           
│   ├── spring5k_2.mf4           
│   ├── spring10k_1.mf4          
│   ├── spring10k_2.mf4          
│   ├── spring20k_1.mf4          
│   └── spring20k_2.mf4          
│
├── scripts/                     # MATLAB scripts for data processing and estimation
│   ├── read_mf4_to_mat.m        # Script to convert .mf4 raw data files into .mat files
│   ├── KF_tuner.m               # Kalman Filter tuning using Particle Swarm Optimization
│   ├── estimator_KF.m           # Linear Kalman Filter estimator implementation
│   ├── estimator_EKF.m          # Adaptive Extended Kalman Filter estimator
│   └── estimator_MC.m           # Monte Carlo based force estimator
│
├── Spring2k_1.mat               # Processed MATLAB data files (.mat format)
├── Spring2k_2.mat
├── Spring5k_1.mat
├── Spring5k_2.mat
├── Spring10k_1.mat
├── Spring10k_2.mat
├── Spring20k_1.mat
├── Spring20k_2.mat
└── README.md                    # Project overview and instructions
```

---

If you want, I can help you integrate this directly into your existing README!



## Table of Contents

- [Background](#background)
- [Data](#data)
- [Installation](#installation)
- [Usage](#usage)
  - [Kalman Filter Tuner](#kalman-filter-tuning-via-particle-swarm-optimization)
  - [Adaptive EKF Estimator](#adaptive-ekf-with-tighter-tuning-and-real-time-smoothing)
  - [Monte Carlo Estimator](#monte-carlo-estimator)
- [Performance Metrics](#performance-metrics)
- [Functions Description](#functions-description)
- [References](#references)

---

## Background

Estimating the interaction force in wearable robotics is crucial for enhancing transparency and safety in human-robot interaction. This repository provides tested MATLAB codes to estimate this force using state estimation techniques based on position, velocity, and acceleration data measured from an experimental platform.

---

## Data

The raw experimental data is originally stored in MF4 format (.mf4), which contains detailed sensor and CAN bus recordings. To facilitate analysis and estimator development, these .mf4 files must be converted to MATLAB .mat files containing the required variables (t, x_r, v_r, x_h, v_h, f_r, f_int, a_r, a_h).

This repository includes a script `read_mf4_to_mat.m` that handles this conversion. Use it to import your raw .mf4 files and export preprocessed .mat files compatible with the estimators.

The estimators and examples provided operate directly on the .mat files, ensuring fast loading and streamlined processing.

---

All algorithms use experimental data files such as:

- `Spring2k_1.mat` — dataset with signals of robot and human positions, velocities, accelerations, and forces.

Make sure these `.mat` files are placed in the working directory before running the scripts.

---

## Installation

1. Clone or download this repository.
2. Open MATLAB and set the working directory to this folder.
3. Ensure required data files are present (`Spring2k_1.mat`, `Spring20k_1.mat`, etc.)
4. Run scripts as needed.

---

## Usage

### Kalman Filter Tuning via Particle Swarm Optimization

Script: `KF_tuner.m`

This script automatically tunes the KF process and measurement noise covariance matrices `Q` and `R` by minimizing the root mean square error (RMSE) of the interaction force estimate using MATLAB's particle swarm optimizer.

**Features:**
- Loads experimental data.
- Sets physical parameters.
- Defines optimization variables and bounds in logarithmic scale.
- Runs PSO to find optimal diagonal elements of `Q` and `R`.
- Runs the KF with optimized parameters.
- Plots interaction force estimates against ground truth.
- Prints error metrics.

**Run with:**
```matlab
KF_tuner
````

---

### Adaptive EKF with Tighter Tuning and Real-Time Smoothing

Script: `EKF_adaptive.m`

Implements an Extended Kalman Filter that estimates the state vector including stiffness and damping parameters (`k` and `c`) adaptively. It updates process and measurement noise covariances online using innovation and residual statistics, and applies a real-time low-pass filter to smooth the estimated interaction force.

**Key points:**

* Uses Runge-Kutta 4 integration in the prediction step.
* Adapts `Q` and `R` every fixed number of steps.
* Provides smoothed interaction force output.
* Plots results and computes error metrics.

**Run with:**

```matlab
EKF_adaptive
```

---

### Monte Carlo Estimator

Script: `MC_estimator.m`

Uses a Monte Carlo approach to estimate stiffness `k` and damping `b` parameters by sampling candidate values, computing errors between predicted and measured interaction forces, and filtering the samples to adaptively refine estimates over time.

**Highlights:**

* Maintains a sliding window of recent measurements.
* Filters parameter samples based on error thresholds.
* Outputs estimates with confidence intervals.
* Compares estimated force against ground truth.

**Run with:**

```matlab
MC_estimator
```

---

## Performance Metrics

For each estimator, the following error metrics are computed and printed:

* **Bias:** Mean error
* **MAE:** Mean absolute error
* **RMSE:** Root mean squared error
* **MaxE:** Maximum absolute error
* **MinE:** Minimum absolute error
* **Statistical test:** Two-sample t-test comparing estimates to ground truth

Plots compare estimated interaction force to measured force over time.

---

## Functions Description

* `kfRmseObj(p, data, m1, m2, k, c)`: Objective function for PSO tuning, returns RMSE of KF force estimates.
* `runKF(Q, R, data, m1, m2, k, c)`: Runs KF with given noise matrices and data, returns estimated interaction force.
* `fProcess_RK4_KF(...)`: Runge-Kutta integrator for KF state prediction.
* `processJacobian_KF(...)`: Computes the continuous-time Jacobian matrix of the system dynamics.
* `measurementJacobian_KF(...)`: Computes measurement Jacobian matrix.
* `hMeas_KF(...)`: Computes predicted measurement from state.
* Equivalent functions are defined for the adaptive EKF with state augmentation.
* Monte Carlo estimator script maintains parameter sampling and error filtering routines.

---

## References

* Boaventura et al., "Acceleration-based Transparent Control Framework for Wearable Robots," *IEEE Transactions on Robotics*, 2016.
* MATLAB Documentation for Particle Swarm Optimization and Kalman Filtering.

---

If you find this repository useful or have questions, feel free to open an issue or contact the author.

---

**Author:** Elisa G. Vergamini and André Vecchione
**Date:** July 2025

---
