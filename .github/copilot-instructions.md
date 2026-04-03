# copilot-instructions for computeSOSFunnels

## 1) Big picture: pipeline + data flow
- Entry point: `main.m`. It runs four stages in order:
  1. `Step1_computeNominalTrajectory.m` -> solves direct collocation for nominal path
  2. `Step2_FeedbackControllerSynthesis.m` -> TVLQR gains + P matrices
  3. `Step3_getDeviationDynamics.m` -> polynomialize dynamics around nominal trajectory
  4. `Step4_computeTimeSampledInvarianceCertificates.m` -> SOS funnel synthesis
- Shared data in `./precomputedData` (`nominalTrajectory.mat`, `LQRGainsAndCostMatrices.mat`, `deviationDynamics.mat`, `taylorApprox_degree3.mat`). Use these saves to avoid repeated expensive solves.
- Modular dynamics:
  - default unicycle in `main.m` with `dynamicsFnHandle = @(x,u) unicycle_dynamics(x,u)`
  - system-specific branches: `cartPole`, `quadrotor` (same pipeline structure, different dynamics/parameters)
- Visualization in `utils/plottingScript.m` and class `lib/plottingFnsClass.m`.

## 2) Characters of the codebase: project conventions
- Scripts are stateful and rely on `exist('var','var')` patterns from wrapper; avoid re-definition across scripts.
- `main.m` supplies overrides such as `initialState`, `finalState`, `maxTimeHorizon`, `numTimeSteps`, `Q`, `R`, `terminalRegionScaling`, `upsamplingFactor`, `maxIter`, `rhoInitialGuessConstant`, `rhoInitialGuessExpCoeff`.
- Non-functional side effects are normal, e.g., `clearvars` at end of Step scripts, then saving artifacts to `.mat`.
- SOS variables use `pvar` (`createPvarVector`) and symbolic conversion through `taylor_approx` chain in `Step3_getDeviationDynamics.m`.

## 3) Workflow commands and debugging
- From MATLAB command line:
  - `cd /home/khalid/git_localRepos/computeSOSFunnels`
  - `main` (runs full pipeline)
  - alternate: `run('Step1_computeNominalTrajectory.m')`, etc.
- To do partial re-runs quickly, load precomputed files:
  - `load('./precomputedData/nominalTrajectory.mat')`
  - `load('./precomputedData/LQRGainsAndCostMatrices.mat')`
- For SOS tuning: set `drawFlag = 0/1`, `maxIter`, and `usageMode` in `main.m` or `Step4`.
- For testing a system branch:
  - `git checkout cartPole` / `git checkout quadrotor`, then run same `main.m` pattern.

## 4) Toolchain and dependencies
- MATLAB R2021b+.
- YALMIP (`yalmip('clear')` at start of Step1)
- IPOPT for direct collocation solver in `lib/getNominalTrajectory_using_DirectCollocation.m`.
- SOSTOOLS + SDP solver (MOSEK recommended, SeDuMi available).
- Symbolic Toolbox used for the Jacobian/Taylor expansions in Step2 and Step3.
- Parallelization: `parfor` in Step3; optional `parpool` in `main.m` (commented-out).

## 5) Source-file map for quick code navigation
- high-level script orchestration: `main.m`
- nominal path: `Step1_computeNominalTrajectory.m` + `lib/getNominalTrajectory_using_DirectCollocation.m`
- TVLQR: `Step2_FeedbackControllerSynthesis.m` + `lib/compute_tvlqr_gains.m`
- deviation dynamics: `Step3_getDeviationDynamics.m`
- SOS certificates: `Step4_computeTimeSampledInvarianceCertificates.m`
- plots: `utils/plottingScript.m`, `lib/plottingFnsClass.m`
- precomputed constants: `precomputedData/*.mat`

## 6) Agent guiding heuristics (do these first)
- Respect the existing `exist('var','var')` override pattern.
- Don’t rewrite to full functional style in one PR; keep backward compatibility with current script sequencing.
- For stability improvements, adjust `Step2` cost weights first (`Q`, `R`, `terminalRegionScaling`) before changing `Step4` SOS init (`rho0`, `c`).
- Keep `addpath('./lib/')` line at top of scripts and avoid relative path breakage.

## 7) Review ask
- If a section (e.g., SOS parameter meaning in `Step4`) is unclear, ask the user for the preferred level of detail and whether to include specific matrix dimensions for each supported branch.
