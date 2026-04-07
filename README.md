This code solves for the Choptuik spacetime in D dimensions using pseudospectral methods, Fourier in time and Chebyshev in the radial direction. 
The equations of motion are given by 
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

