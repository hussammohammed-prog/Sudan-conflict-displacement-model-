# Sudan Conflict-Related Displacement Model

This repository contains the MATLAB code associated with a dynamical model of conflict-related displacement between Khartoum and Port Sudan during May 2023–April 2024.

## Model Overview

The model represents displacement dynamics using two interconnected compartments:

- Khartoum
- Port Sudan

It incorporates:

- Time-dependent conflict forcing
- Bidirectional population movement
- Destination capacity
- External source terms
- Removal processes
- Parameter estimation
- Numerical simulation
- Local sensitivity analysis
- Global sensitivity analysis
- Displacement Pressure Index (DPI)

## Data

The model uses available monthly displacement observations and conflict-event information for the study period.

Missing observations are not imputed and are excluded from the corresponding calibration and performance calculations.

Direct Port Sudan observations are limited, and broader Red Sea data are treated as contextual information rather than direct Port Sudan observations.

## Numerical Methods

The model is implemented in MATLAB.

The code performs:

1. Data preparation
2. Conflict forcing construction
3. Parameter estimation
4. Numerical solution of the dynamical system
5. Model-performance evaluation
6. Sensitivity analysis
7. DPI calculation
8. Generation of numerical results and figures

## Reproducibility

The code is provided to facilitate reproduction of the numerical results reported in the associated manuscript.

Results depend on the model structure, input data, parameter bounds, and calibration settings.

## Repository Contents

- `README.md` — Description of the model and repository
- `Sudan_Displacement_Model.m` — Main MATLAB implementation
- `data/` — Input data
- `results/` — Numerical results
- `figures/` — Generated figures

## Software

- MATLAB
