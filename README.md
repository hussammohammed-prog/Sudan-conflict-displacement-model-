# Sudan Conflict-Related Displacement Model

This repository contains the MATLAB code for a dynamical model of conflict-related displacement between Khartoum and Port Sudan, Sudan.

## Study Period

May 2023 – April 2024

## Data

The model uses monthly displacement and conflict-event data for the study period.

The available displacement data include:

- Khartoum monthly displacement observations.
- Limited direct Port Sudan observations.
- Broader Red Sea displacement data used as contextual information.
- Monthly conflict-event information used to construct the conflict forcing.

Missing observations are not imputed and are excluded from the corresponding calculations.

## Model

The model represents displacement dynamics between Khartoum and Port Sudan and includes:

- Time-dependent conflict forcing
- Population movement between the two locations
- Destination capacity
- External source terms
- Removal processes

## Numerical Analysis

The MATLAB code performs:

- Parameter estimation
- Numerical simulation
- Model performance evaluation
- Local sensitivity analysis
- Global sensitivity analysis
- Displacement Pressure Index (DPI)

## Main File

`Sudan_Displacement_Model.m`

This file contains the main MATLAB implementation of the model and numerical analysis.

## Reproducibility

The code is provided to support reproducibility of the numerical results reported in the associated manuscript.

The results depend on the available data, model assumptions, parameter bounds, and calibration settings.

## Software

MATLAB
