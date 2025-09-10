# Simulation Workflow

The simulation codes are written in **MATLAB**. The workflow is organized into sequential files named **`mainstep1` – `mainstep10`**.  

## File Descriptions

- **`mainprepare.m`**  
  Generates the synthetic data $\mathbf{Z}$ and $\mathbf{U}$.

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

- **`mainstep6.m`**  
  Performs for Aim 2 and gives the test for the four cases.

- **`mainstep7.m`**  
  Reports the results for Aim 2.

- **`mainstep8.m`**  
  A small experiments for the comparison of optimization time.

- **`mainstep9.m`**  
  Gives the comparison in hypothesis testing in Section S3.2.

- **`mainstep10.m`**  
  Reports the results in Section S3.2 .

- **`mainunknownprepare.m`**  
  Generates the synthetic data $\mathbf{Z}$ and $\mathbf{U}$.

- **`mainunknownstep1.m`**  
  Performs regression to estimate matrix $\boldsymbol{\Gamma}$.

- **`mainunknownstep2.m`**
  Performs regression.
 
- **`mainunknownstep3.m`**
  Computes the debiased matrix.

- **`mainunknownstep4.m`**  
  Computes the debiased estimator.

- **`mainunknownstep5.m`**  
  Reports the results.

## Usage

Run the files in the following order:

1. `mainprepare.m`  `mainunknownprepare.m` 
2. `mainstep1.m` → `mainstep10.m` (sequentially)
3. `mainunknownstep1.m` → `mainunknownstep5.m` (sequentially)

This will generate synthetic data, run the regressions, compute debiased estimators, and produce the final simulation results. The corresponding .sh file helps us to submit the job to High Performance Computation Platform.
