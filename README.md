# Solar Spectropolarimetry & Optical Deconvolution
### Analysis of Diffraction Effects on Solar Magnetic Field Estimation via Weak Field Approximation (WFA) from SUNRISE III / TuMag Data

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-3776ab?style=flat-square&logo=python)](https://www.python.org/)
[![Astropy](https://img.shields.io/badge/Astropy-Solar%20Astrophysics-orange?style=flat-square&logo=python)](https://www.astropy.org/)
[![SciPy](https://img.shields.io/badge/SciPy-Numerical%20Optimization-blue?style=flat-square&logo=scipy)](https://scipy.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](LICENSE)
[![UGR](https://img.shields.io/badge/Institution-Universidad%20de%20Granada-red?style=flat-square)](https://www.ugr.es/)
[![IAA-CSIC](https://img.shields.io/badge/Collaborator-IAA--CSIC-blue?style=flat-square)](https://www.iaa.csic.es/)

> **Bachelor's Thesis in Physics (Trabajo de Fin de Grado)** — *Grade: 9.4 / 10*  
> **Author**: [Pablo Orellana Chornyak](https://github.com/porellanaphysics) ([LinkedIn](https://linkedin.com/in/p-orellana) · [Email](mailto:pablo.orellana.physics@gmail.com))  
> **Academic Profile**: Graduado en Física (Universidad de Granada) | Estudiante del Máster en Ciencias y Tecnologías desde el Espacio (Universidad de Alcalá)  
> **Collaborating Institution**: Instituto de Astrofísica de Andalucía (IAA-CSIC), Solar Physics Department  
> 📄 **Thesis Monograph**: [`docs/main.pdf`](docs/main.pdf) (13 MB)  
> 📊 **Defense Slides**: [`docs/PresentacionFinal.pdf`](docs/PresentacionFinal.pdf) (24 MB)  

---

## Executive Summary

Solar spectropolarimetry allows astronomers to reconstruct the magnetic topology of the solar photosphere. However, even with balloon-borne telescopes flying in the stratosphere above 99% of the atmosphere—such as the **SUNRISE III** observatory and its **TuMag** (Tunable Magnetograph) instrument—diffraction and residual optical aberrations introduce spatial blur, degrading fine magnetic flux tubes in the **"Quiet Sun"**.

This project investigates the impact of optical degradation on solar magnetic field inferences using the **Weak Field Approximation (WFA)** of the Radiative Transfer Equation (RTE). 

Because classical image deconvolution methods (e.g. Richardson-Lucy) require non-negative signals ($I \ge 0$), they **fail catastrophically when applied directly to circular polarization Stokes $V$** due to its inherently bipolar (positive/negative) nature and high noise sensitivity.

To resolve this limitation, this thesis implements a **spatially coupled forward modeling inversion** (building upon the theoretical framework of van Noort, 2012). By coupling telescope optical modeling (Point Spread Function, PSF) with radiative transfer physics, and solving the resulting inverse problem via a **matrix-free Conjugate Gradient (CG)** solver, this algorithm reliably reconcentrates dispersed magnetic flux without noise explosion.

---

## Visual Gallery & Results

### 1. Spatial Image Degradation Model
$$\text{Observed Image } I_{\mathrm{obs}} = I_{\mathrm{true}} \ast \mathrm{PSF} + \mathcal{N}(\sigma)$$

![Diffraction & PSF Convolution](docs/img/ejemplo_degradacion.png)
*Spatial degradation pipeline: Ground-truth scene convolved with optical PSF and corrupted by noise.*

---

### 2. Magnetic Field Reconstruction Benchmark ($B_\parallel$)

| Ground Truth vs. Degraded vs. Coupled Restoration |
| :---: |
| ![Magnetic Field Zoom Comparison](docs/img/exp4_3_campoB_zoom.png) |
| *Zoom-in on fine magnetic structures: Comparison between true $B_\parallel$, blurred observation, classical uncoupled deconvolution, and our coupled forward model.* |

---

## Mathematical & Physical Foundations

### 1. Stokes Parameters & The Weak Field Approximation (WFA)
Polarized light from solar spectral lines is described by the Stokes vector $\mathbf{S} = [I, Q, U, V]^T$:
- $I$: Total intensity.
- $Q, U$: Linear polarization states.
- $V$: Circular polarization (sensitive to longitudinal magnetic field $B_\parallel$).

Under the **Weak Field Approximation (WFA)**, when Zeeman splitting $\Delta \lambda_B = 4.67 \times 10^{-13} \, g_{\mathrm{eff}} \, \lambda_0^2 \, B$ is much smaller than the Doppler line width ($\Delta \lambda_B \ll \Delta \lambda_D$), circular polarization $V$ is directly proportional to the line-of-sight magnetic field $B_\parallel$ and the spectral derivative of intensity:

$$V(\lambda) = -C \, g_{\mathrm{eff}} \, \lambda_0^2 \, B_{\parallel} \, \frac{\partial I(\lambda)}{\partial \lambda}$$

where:
- $C = 4.67 \times 10^{-13}\text{ \AA}^{-1}\text{ G}^{-1}$ (gravitational/atomic constant).
- $g_{\mathrm{eff}}$: Effective Landé factor of the spectral transition (e.g., Fe I 5250.2 \AA, Fe I 6173 \AA, Mg I b2 5173 \AA).
- $\lambda_0$: Central line wavelength.
- $B_\parallel$: Longitudinal magnetic field strength in Gauss.

---

### 2. The Optical Inverse Problem & Why Standard Deconvolution Fails
Telescope optics degrade spatial information via convolution with the instrumental Point Spread Function (PSF):

$$I_{\mathrm{obs}}(x, y, \lambda) = [I \ast \mathrm{PSF}](x, y, \lambda) + \mathcal{N}$$

$$V_{\mathrm{obs}}(x, y, \lambda) = [V \ast \mathrm{PSF}](x, y, \lambda) + \mathcal{N}$$

Applying classic deconvolution (Richardson-Lucy, Maximum Likelihood) independently to $V_{\mathrm{obs}}$ fails because:
1. Richardson-Lucy is strictly formulated for Poisson processes and positive distributions; Stokes $V$ is signed ($V \in \mathbb{R}$).
2. The ratio $\frac{V}{\partial I / \partial \lambda}$ causes catastrophic noise amplification near spectral line cores and continuum regions where $|\partial I / \partial \lambda| \to 0$.

---

### 3. Spatially Coupled Inversion (Van Noort Scheme)
Instead of deconvolving $V$ directly, we formulate the forward model:

$$V_{\mathrm{model}}(x, y, \lambda) = \mathbf{P} \left[ -C \, g_{\mathrm{eff}} \, \lambda_0^2 \, B_{\parallel}(x, y) \, \frac{\partial I_{\mathrm{deconv}}(x, y, \lambda)}{\partial \lambda} \right]$$

where $\mathbf{P}$ represents the spatial convolution operator with the optical PSF.

We minimize the regularized quadratic loss functional:

$$\chi^2(B_\parallel) = \sum_{\lambda, x, y} \left| V_{\mathrm{obs}} - V_{\mathrm{model}}(B_\parallel) \right|^2 + \gamma \, \mathcal{R}(B_\parallel)$$

where $\mathcal{R}$ is a Tikhonov / spatial smoothness regularizer.

#### Matrix-Free Conjugate Gradient Optimization
Storing the full Jacobian matrix $\mathbf{J}$ for a $1600 \times 1600$ pixel detector across multiple spectral wavelengths would require terabytes of RAM. 

To overcome this, we implement a **matrix-free linear operator** $\mathbf{A} = \mathbf{J}^T \mathbf{J} + \gamma \mathbf{I}$ via `scipy.sparse.linalg.LinearOperator`. At each iteration of the Conjugate Gradient (CG) algorithm:
1. Forward evaluation: compute $\mathbf{P} [ K(x, y, \lambda) \cdot B(x, y) ]$ using fast 2D FFTs.
2. Adjoint evaluation: compute the adjoint operator $\mathbf{P}^T$ using reflected padding and conjugated FFT kernels.

This allows full-detector inversions in seconds with $\mathcal{O}(N)$ memory footprint.

---

## Repository Architecture

```
.
├── src/                            # Core scientific and numerical algorithms
│   ├── magnetismo.py               # WFA calculation & coupled van Noort inversion
│   ├── deconvolucion.py            # Airy PSF generator, 3D Richardson-Lucy, Wiener filter
│   ├── visualizacion.py            # Spectropolarimetric plotting & magnetogram visualizer
│   └── probarSis.py                # Synthetic experiment pipeline & metric evaluation
├── tests/                          # Unit and verification tests
│   ├── test_adjoint_operator.py    # Adjoint property verification: <Px, y> == <x, P^T y>
│   └── test_sunspot_norms.py       # Linear operator norm and gradient consistency tests
├── docs/                           # Academic thesis and presentation
│   ├── main.pdf                    # Complete Bachelor's Thesis Monograph (PDF)
│   ├── PresentacionFinal.pdf       # Thesis Defense Slides (PDF)
│   ├── main.tex                    # LaTeX monograph master document
│   ├── content/                    # Individual monograph chapters (01 to 09)
│   ├── img/                        # Publication figures, plots, and diagrams
│   └── referencias.bib             # BibTeX reference database
├── data/                           # (Git-ignored) Observational FITS files & datasets
├── generar_ejemplo_psf.py          # Didactic script generating diffraction figure
├── requirements.txt                # Lightweight Python package requirements
├── requirements-freeze.txt         # Full reproducibility lockfile
├── LICENSE                         # MIT License
└── README.md
```

---

## Key Achievements & Findings

1. **Algorithms Programmed from Scratch**: All classical deconvolution schemes (Richardson-Lucy with Wiener pre-filtering, apodization, damping) and the coupled forward-modeling algorithm were implemented without black-box libraries.
2. **Reconcentration of Magnetic Flux**: On synthetic solar models, the coupled method recovers up to **92% of the ground-truth peak magnetic flux** in narrow granules, compared to only **58%** in raw degraded observations.
3. **Robustness to Noise**: The coupled approach prevents high-frequency noise amplification by anchoring the magnetic inversion to the deconvolved intensity gradient $\partial I / \partial \lambda$.
4. **Real TuMag Data Application**: Applied to quiet Sun observations from SUNRISE III, sharpening intergranular magnetic concentrations and confirming the theoretical predictions.

---

## Quickstart & Usage

### 1. Installation
Clone the repository and install dependencies:
```bash
git clone https://github.com/<your-username>/solar-weak-field-deconvolution.git
cd solar-weak-field-deconvolution
pip install -r requirements.txt
```

### 2. Run Verification Tests
Verify mathematical consistency and adjointness of the spatial operators:
```bash
python3 tests/test_adjoint_operator.py
```

### 3. Generate Didactic Degradation Figure
```bash
python3 generar_ejemplo_psf.py
```

---

## Citation & Academic Reference

If you use this work or its methodology in your research, please cite:

```bibtex
@mastersthesis{Orellana2026TFG,
  author       = {Pablo Orellana Chornyak},
  title        = {Análisis del efecto de la difracción en la estimación del campo magnético solar mediante la aproximación de campo débil (WFA) a partir de datos del instrumento TuMag a bordo de Sunrise III},
  school       = {Universidad de Granada},
  year         = {2026},
  note         = {Trabajo de Fin de Grado en Física. En colaboración con el Instituto de Astrofísica de Andalucía (IAA-CSIC)}
}
```

---

## License

This project is licensed under the [MIT License](LICENSE).