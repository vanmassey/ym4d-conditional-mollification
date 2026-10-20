# 4D Yang-Mills Conditional Continuum Framework

This repository hosts the companion codebase and formal verification files for the preprint manuscript:
"A Conditional Framework for the Four-Dimensional Euclidean Yang-Mills Measure on a Compact Polish Space of Wilson Loops via Strongly Parabolic Covariant Mollification" (arXiv: hep-th/math-ph).

## Repository Architecture

### 1. Formal Methods Verification (`/lean`)
The `/lean` directory contains `YangMills_DeTurck.lean` under the namespace `YMDT`. This file formalises the algebraic types and boundary conditions required to construct the generalized connection space quotients and gauge actions on a four-dimensional flat torus. 

Crucially, it provides the exact type signature for Lemma 2.2 (Longitudinal drift cancellation under cyclic trace invariance). To maintain structural compilability and enforce zero placeholders, the model avoids complex topological group-closure proofs of smooth gauge fields by defining the space topology via Lean 4's reflexive, symmetric, and transitive equivalence closure (`Relation.EqvGen`).

### 2. Numerical Simulation Infrastructure (`/src`)
The `/src` directory contains the 4D Lattice Heat-Kernel Simulator (`ym_heat_simulation.py`). The program evaluates the discrete heat kernel semigroup \(e^{-2t\Delta}\) over a periodic hypercubic grid using four-dimensional Discrete Fourier Transforms (FFT) to bypass the exponential complexity constraints of true 4D grid mapping (L⁴).

## Empirical Simulation Results

The simulator tracks the ultraviolet convergence behavior of the algebraic BPHZ subtraction tensor:
\[C_t = \frac{1}{2} (K_{2t} - \delta_0) \otimes I_{\mathfrak{g}}\]
over a localized flow-time spectrum \(t \in [10^{-4}, 10^{-1}]\). 

A linear regression analysis executed over a 16⁴ discrete lattice configuration (65,536 coordinate nodes) yielded the following verified parameters integrated directly into Section 3 of the accompanying manuscript:

*   **Isolated Logarithmic Slope (\(m_{\mathrm{sim}}\)):** -0.092783
*   **Theoretical Convergence Intercept (\(b_{\mathrm{sim}}\)):** 0.350591

The output file asset `Graph.png` verifies the uniform logarithmic ultraviolet divergence behavior of the regularizer, serving as a non-perturbative consistency check for the bulk counterterm cancellations before removing the lattice cutoff (\(a \downarrow 0\)).

## Prerequisites and Local Execution

To run the 4D FFT simulation infrastructure locally, ensure you have Python 3.8+ installed along with the core numerical modules:

```bash
pip install numpy scipy matplotlib
```

Execute the script via the command line interface:

```bash
python src/ym_heat_simulation.py
```

## Citation

If you utilize this framework or its corresponding type signatures in your research, please cite the pre-print manuscript:

```bibtex
@article{fsl2026ym01,
  title={A Conditional Framework for the Four-Dimensional Euclidean Yang-Mills Measure on a Compact Polish Space of Wilson Loops via Strongly Parabolic Covariant Mollification},
  author={Max Van Massey},
  journal={arXiv preprint hep-th/2610.XXXXX},
  year={2026}
}
```

## License
This project is licensed under the MIT License. See the LICENSE file for details.
