# Simulation Workflow

The simulation codes are written in **MATLAB**. The workflow is organized into sequential files named **`mainstep1` – `mainstep10`**.  

## File Descriptions

- **`mainprepare.m`**  
  Generates the synthetic data $\mathbf{Z}$ and $\mathbf{Z}$.

- **`mainstep1.m`**  
  Performs regression with either fixed $\lambda$ or cross-validation.

- **`mainstep2.m`**  
  Computes the debiased matrix.

- **`mainstep3.m`**  
  Computes the debiased estimator.

- **`mainstep4.m`**  
  Reports the results of our proposed methods.

- **`mainstep5.m`**  
  Runs OLS regression with the support set known (oracle case).

## Usage

Run the files in the following order:

1. `mainprepare.m`  
2. `mainstep1.m` → `mainstep5.m` (sequentially)

This will generate synthetic data, run the regressions, compute debiased estimators, and produce the final simulation results.
