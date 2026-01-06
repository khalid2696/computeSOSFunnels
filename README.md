# computeSOSFunnels

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18047554.svg)](https://doi.org/10.5281/zenodo.18047554)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-orange.svg)](https://www.mathworks.com/)

**A MATLAB Framework for Invariant Funnel Synthesis using Sum-of-Squares Optimization.**

This repository implements a modular pipeline to compute invariant funnels around nominal trajectories for nonlinear dynamical systems. It integrates trajectory optimization, time-varying LQR (TVLQR) stabilization, and Sum-of-Squares (SOS) programming to generate rigorous Lyapunov certificates of stability.

---

## 🚀 Key Features
* **Trajectory Optimization:** Direct collocation (using IPOPT) to compute nominal trajectories and feedforward inputs.
* **Feedback Control:** Computation of TVLQR gains and solution to the Riccati differential equation (useful in initializing the SOS program).
* **Taylor Approximation:** Automated symbolic expansion of nonlinear dynamics into polynomial form (required for SOS verification).
* **Funnel Synthesis:** Bilinear alternation scheme to compute invariant funnels via SOS programming.
* **Multi-System Support:** Ready-to-use implementations for:
  * Unicycle (Default)
  * Cart-Pole (Branch: `cartPole`)
  * Quadrotor (Branch: `quadrotor`)

## 📄 Documentation
A comprehensive technical report and user tutorial is available:
**[Download the Tutorial PDF](https://github.com/khalid2696/computeSOSFunnels/releases/latest/download/computeSOSFunnels_Tutorial.pdf)**

## 🛠️ Prerequisites
* **MATLAB** (R2021b or later recommended)
* **[YALMIP](https://yalmip.github.io/)** (Optimization interface)
* **[IPOPT](https://coin-or.github.io/Ipopt/)** (Nonlinear solver)
* **[SOSTOOLS](https://www.cds.caltech.edu/sostools/)** (SOS parsing)
* **SDP Solver:** [MOSEK](https://www.mosek.com/) (Recommended) or SeDuMi.

## 📦 Installation
```bash
git clone https://github.com/khalid2696/computeSOSFunnels.git
cd computeSOSFunnels
```
Open MATLAB and run `main.m`. The script handles path setup automatically.

## 🏃 Quick Start
1.  Open `main.m` in MATLAB.
2.  Run the script to execute the pipeline:
    * Step 1: Compute Nominal Trajectory
    * Step 2: Synthesize TVLQR Controller
    * Step 3: Polynomialize Deviation Dynamics
    * Step 4: Compute SOS Funnels
3.  Visualize results using `utils/plottingScript.m`.

## 🤖 Supported Dynamical Systems

This framework is designed to be modular, versatile, and system-agnostic. It is currently implemented and tested on three distinct nonlinear systems, each maintained on a separate branch to serve as a reference template.

| System | Branch Name | Complexity | Description |
| :--- | :--- | :--- | :--- |
| **Unicycle** | `main` | 3 States, 2 Inputs | A non-holonomic mobile robot model (Default). |
| **Cart-Pole** | `cartPole` | 4 States, 1 Input | A classic underactuated benchmark system. |
| **Quadrotor** | `quadrotor` | 12 States, 4 Inputs | A high-dimensional system with complex nonlinear dynamics. |

### Switching Systems
To use a specific system, simply checkout the corresponding branch. For example, to use the quadrotor system:

```bash
git fetch origin
git checkout quadrotor
```

## 📧 Contact & Support

If you have questions about the codebase, installation issues, or encounter bugs, please **open a GitHub Issue** &mdash; this would allow the community to benefit from the discussions.

For research inquiries, collaboration proposals, or theoretical questions regarding the Funnel Synthesis framework, feel free to reach out:

> **Mohamed Khalid M Jaffar**  
> 🌐 [khalid2696.github.io](https://khalid2696.github.io)

## 🖊️ Citation
If you find this code useful in your research, please cite the software:

**BibTeX:**
```bibtex
@software{jaffar_2025_computeSOSFunnels,
  author       = {M Jaffar, Mohamed Khalid},
  title        = {computeSOSFunnels: A MATLAB Framework for Invariant Funnel Synthesis using Sum-of-Squares Optimization},
  version      = {1.0.0},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18047554},
  url          = {https://doi.org/10.5281/zenodo.18047554},
  year         = {2025}
}
```

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
