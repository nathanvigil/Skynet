# 1D Multi-Region Reactor Core Spatial Diffusion Simulator

## Project Overview
This repository contains a computational physics simulator that models the spatial distribution of neutron flux across a heterogeneous, multi-region 1D nuclear reactor core. While standard point kinetics are useful for time transients, they fail to account for spatial variations across different material zones. This simulator bridges that gap by resolving the steady-state neutron diffusion equation numerically using the **Finite Difference Method (FDM)** and solving a generalized eigenvalue system to determine the core's effective multiplication factor ($k_{\text{eff}}$).

### Key Features
* **Multi-Region Geometry:** Configurable domains mapping an enriched fuel core flanked by low-absorption material reflectors.
* **Numerical Matrix Engine:** Leverages `scipy.linalg.eig` to solve the generalized eigenvalue problem ($A\vec{\phi} = \frac{1}{k_{\text{eff}}}B\vec{\phi}$).
* **Automated Data Filtering:** Features algorithmic masking to filter out non-physical harmonic modes and automatically extract the fundamental physical profile.
* **Engineering Dashboard:** Generates side-by-side visualizations comparing spatial material cross-sections directly against the resulting continuous neutron flux profile.

---

## Mathematical & Physical Framework

The simulator is governed by the 1-group thermal neutron diffusion equation, which serves as a localized neutral-particle balance sheet:

$$-D \frac{d^2\phi(x)}{dx^2} + \Sigma_a \phi(x) = \frac{1}{k_{\text{eff}}} \nu \Sigma_f \phi(x)$$

Where:
* $-D \frac{d^2\phi(x)}{dx^2}$ represents **Net Neutron Leakage** (spatial divergence).
* $\Sigma_a \phi(x)$ represents **Macroscopic Absorption Losses**.
* $\nu \Sigma_f \phi(x)$ represents the **Physical Production Rate** of neutrons via fission.
* $k_{\text{eff}}$ acts as the eigenvalue scaling factor required to force the system into a perfect steady-state mathematical balance.

### Boundary Conditions
The core boundaries are simulated using a **Vacuum Boundary Condition**, where the flux is forced to zero at the absolute outer edges of the reflector zones:
$$\phi(0) = 0, \quad \phi(x_{\text{end}}) = 0$$

---

## Computational Methodology

Because macroscopic cross-sections ($D, \Sigma_a, \Sigma_f$) drop or jump abruptly at the fuel-reflector interfaces, analytical pen-and-paper solutions are highly impractical. This simulator resolves the problem via spatial discretization:

1.  **Grid Assembly:** Slices the continuous 1D spatial domain into a discrete grid of $N$ nodes separated by a spatial step $\Delta x$.
2.  **Central Difference Stencil:** Replaces the continuous second-derivative leakage term with an algebraic finite difference approximation:
    $$\frac{d^2\phi}{dx^2} \approx \frac{\phi_{i+1} - 2\phi_i + \phi_{i-1}}{(\Delta x)^2}$$
3.  **Matrix Formulation:** Maps individual node balances into a global tridiagonal loss matrix $A$ and a diagonal fission production matrix $B$, solving for the fundamental eigenvector ($\vec{\phi}$) and dominant eigenvalue ($k_{\text{eff}}$).

---

## Results & Visualization

When executed with a $60\text{ cm}$ fuel core enveloped by $20\text{ cm}$ water reflectors, the simulator successfully converges on a supercritical system multiplication factor:

$$\text{Calculated } k_{\text{eff}} \approx 1.9134$$

<img width="1183" height="495" alt="Screenshot 2026-06-02 at 8 39 06 PM" src="https://github.com/user-attachments/assets/3563c295-23c3-477b-9b09-514df64a4ecd" />


### Core Spatial Profiles
The script generates an engineering dashboard that contrasts the structural setup against the subatomic behavior:

* **Left Plot (Cross-Sections):** Showcases the abrupt step-functions marking the exact boundaries where the fuel assembly transitions into the pure water reflector zones.
* **Right Plot (Neutron Flux):** Illustrates the smooth, physical response of the neutron population—cresting in the center where fission production dominates, and leaking down toward the boundaries.



---

## How to Run the Code

### Prerequisites
Ensure you have Python installed alongside the following scientific computing packages:
```bash
pip install numpy scipy matplotlib
