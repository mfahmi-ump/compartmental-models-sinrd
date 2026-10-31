# COVID-19 Compartmental Models: Estimation and Simulation (Malaysia)

MATLAB code for fitting and simulating SIR-family epidemic models on Malaysian COVID-19 data. It covers the full path from raw daily data to fitted parameters to simulated trajectories:

1. **Estimation**: Bayesian (DRAM-MCMC) parameter estimation for an SI4RD model on automatically detected time segments.
2. **Simulation**: five models of increasing complexity (SIR, SIRD, SI3RD, SI4RD, SI4RDV), each solved with three numerical schemes and compared with the data.

## Models

| Script | Model | Compartments |
|---|---|---|
| `sir_Malaysia.m` | SIR | S, I, R (R includes deaths) |
| `sird_Malaysia.m` | SIRD | S, I, R, D |
| `si3rd_Malaysia.m` | SI3RD | S, I1-I3, R, D |
| `si4rd_Malaysia.m` | SI4RD | S, I1-I4, R, D |
| `si4rdv_Malaysia.m` | SI4RDV | SI4RD plus a time-varying vaccinated fraction `pvax` and natural birth/death (`mu = 2.1e-5`) |

In the multi-class models the infectious compartments I1, I2, I3, I4 are severity stages (infected, hospitalised, ICU, ventilated). Each infectious class has its own transmission rate (`beta`), recovery rate (`r`) and fatality rate (`d`), and a share `p` of new infections. Parameters change day by day and are read from the estimate files produced by the estimation script.

## Numerical methods

Each simulation script solves the model three ways and reports the RMSE of each against the observed data:

- **RK4**: classical fourth-order Runge-Kutta.
- **Adams-Bashforth-Moulton**: fourth-order predictor-corrector.
- **SRK4**: four-stage stochastic Runge-Kutta with multiplicative noise (`sigma = 0.05`), run over 500 paths. The 25th, 50th and 75th percentiles give the median and uncertainty band.

`si3rd_Malaysia.m` runs RK4 and Adams-Bashforth-Moulton only (its SRK4 section is commented out). `sir_Malaysia.m` evaluates the SIR model over sliding windows (lengths 3-90 days, every start day) and stores the RMSE of each, instead of producing a single fit.

## Parameter estimation

`si4rd_parameter_estimation.m` fits the SI4RD model:

1. Detects change points from peaks in the smoothed second derivative of I1 and I3, splitting the series into segments.
2. Builds data-driven starting values for the 18 parameters of each segment.
3. Refines them with `fminsearch` (least squares).
4. Samples the posterior with DRAM (`mcmcrun`) and takes the posterior means.
5. Simulates with the posterior means and plots data against model.
6. Stores the per-day estimates in `paramestimate` (18 x T).

The 18 parameters are `p1-p4`, `beta1-beta4`, `r1-r4`, `d1-d4`, `gamma` and `pvax`.

## Repository structure

```
.
├── si4rd_parameter_estimation.m   # MCMC parameter estimation (SI4RD)
├── si4rdss.m                      # sum-of-squares objective used by the estimation
├── si4rdsrk4.m                    # RK4 solver used by the estimation
├── sir_Malaysia.m                 # SIR simulation, sliding-window RMSE
├── sird_Malaysia.m                # SIRD simulation
├── si3rd_Malaysia.m               # SI3RD simulation
├── si4rd_Malaysia.m               # SI4RD simulation
├── si4rdv_Malaysia.m              # SI4RD with vaccination
├── data/                          # input CSV files (see below)
└── README.md
```

## Requirements

- MATLAB R2019b or later (`readmatrix`, `datetime`, `tiledlayout`).
- Estimation script only: Signal Processing Toolbox (`findpeaks`), [MCMC toolbox](https://mjlaine.github.io/mcmcstat/) by Marko Laine (`mcmcrun`, `chainstats`), and [Panel](https://www.mathworks.com/matlabcentral/fileexchange/20003-panel) by Ben Mitch.
- The simulation scripts need no extra toolboxes.

## Input data

Scripts read CSV files from the current folder. Column layouts, as used by the code:

| File | Columns |
|---|---|
| `sir_data_Malaysia.csv`, `sird_data_Malaysia.csv` | `t, R, I, D, S` |
| `si3rd_data_Malaysia.csv` | `t, R, I1, D, I2, I3, S` |
| `si4rd_data_Malaysia.csv`, `si4rdv_data_Malaysia.csv`, `si4rd_data.csv` | `t, I1, R, D, I2, I3, I4, S, H, ..., ES (col 11)` |

| Parameter file | Columns |
|---|---|
| `Malaysia-SIR-param-estimate.csv` | `t, beta, r` |
| `Malaysia-SIRD-param-estimate.csv` | `t, beta, r, d, gamma` |
| `Malaysia-SI3RD-param-estimate-det.csv` | `t, beta1-3, r1-3, d1-3, p1-3, gamma` |
| `Malaysia-SI4RD-param-estimate.csv` | `t, p1-4, beta1-4, r1-4, d1-4, gamma` |
| `Malaysia-SI4RDV-param-estimate.csv` | as SI4RD, with `pvax` in column 20 |

Population is set to `Ptotal = 3.153e7`. Plot dates start on 24 March 2020 (SI3RD, SI4RD, SI4RDV), 26 March 2020 (SIR, SIRD) and 24 February 2021 (estimation script).

## Usage

```matlab
% 1. (optional) estimate parameters
si4rd_parameter_estimation

% 2. simulate with fitted parameters
si4rd_Malaysia
```

Each simulation script produces a figure comparing data with the three solvers and an `RMSE` matrix (rows: SRK4, RK4, Adams-Bashforth-Moulton).

## License

MIT (see `LICENSE`).

## Citation

If you use this code, please cite the associated paper: *[add reference / DOI here]*.
