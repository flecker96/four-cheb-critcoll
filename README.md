## 🚧 Project Status: Work in Progress

This project is actively under development. Core functionality is implemented, but documentation and some features are still being finalized.

This code solves for the Choptuik spacetime in D dimensions using a pseudospectral collocation method, Fourier in time and Chebyshev in the radial direction. 
The equations of motion that are solved are given by 
## Equations of motion

**E1 (constraint):**

$$2y \partial_y\Omega = -\Omega \big(D-1 - y(\Pi^2 + y \Psi ^2)\big) - y\Omega ^2 + (D-3)(\Pi^2+y\Psi^2) $$
   
**E2:**

$$ (1-y)\partial_y F = F + \frac{1}{2}\left(1+(1-y)F\right)\Omega $$

**E3:**

$$\big(\partial_\tau + 2y\partial_y\big)\Pi = \big(1+(1-y)F\big)\big((D-1 + y\Omega)\Psi + 2y\partial _y \Psi \big) -\Pi$$

**E4:**

$$\big(\partial_\tau + 2y\partial_y\big)\Psi = \big(1+(1-y)F\big) \big(2\partial_y \Pi + \Omega\Pi\big) - 2\Psi$$

**E5 (only monitored):**

$$\big(\partial_\tau + 2y\partial_y\big)\Omega = 2\big(D-3+y\Omega\big)\big(1+(1-y)F\big)\Pi\Psi-2\Omega$$

see also Eq. (76) in https://arxiv.org/pdf/2602.10185. Additionally, a change of variables $f(y)\to F(y)$ where 
$f(y)=1+(1-y)F(y)$ was performed such that $f(y=1)=1$ manifestly. 

## Features
- Transformation to Chebyshev and Fourier space by FFT (fftw package)
- Minimizing residual by Newton-Raphson algorithm (LAPACK)
- All derivatives are implemented in spectral space
- Anti-aliasing implemented by doubling modes before transformation
- Chebyshev collocation grid $z\in [-1,1]$ mapped linearly to $y\in [0,1]$
- Data stored in HDF5 format

## Planned
- Improve stability at $D>6$, investigate solution branches
- Implement Krylov solver to make higher resolutions accessible
- Add full documentation

## Usage

```bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config \
    libfftw3-dev liblapacke-dev \
    libopenmpi-dev openmpi-bin \
    libomp-dev

