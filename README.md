This code solves for the Choptuik spacetime in D dimensions using pseudospectral methods, Fourier in time and Chebyshev in the radial direction. 
The equations of motion are given by 
## Equations of motion

**E1:**

$$y\partial_y\Omega + y\Omega^2 + (d-1)\Omega - y\Omega\left(\Pi^2 + y\Psi^2\right) + (d-3)\left(\Pi^2 + y\Psi^2\right) = 0$$

**E2:**

$$-F + \frac{1}{2}\left(1+(y-1)F\right)\Omega - (y-1)\partial_y F = 0$$

**E3:**

$$\Pi + 2y\partial_y\Pi + \partial_\tau\Pi - \left(1+(y-1)F\right)\left((d-1+y\Omega)\Psi + 2y\partial_y\Psi\right) = 0$$

**E4:**

$$2\Psi + 2y\partial_y\Psi + \partial_\tau\Psi + \left(1+(y-1)F\right)\left(\Omega\Pi + 2\partial_y\Pi\right) = 0$$

**E5:**

$$2\Omega + 2y\partial_y\Omega + \partial_\tau\Omega + 2(d-3)\left(1+(y-1)F\right)\Pi\Psi + 2y\Omega\left(1+(y-1)F\right)\Pi\Psi = 0$$
see also Eq. (76) in https://arxiv.org/pdf/2602.10185 which the additional the change of variables $f(y)\to F(y)$ where 
$f(y)=1+(1-y)F(y)$ such that $f(y=1)=1$ manifestly. 

