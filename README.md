# Driven Damped Oscillator
This project is organized into four main topics:
- Free Damped Oscillator
- Driven Damped Oscillator
- Q-factor analysis
- Energy analysis of the driven damped system
 
## Free Damped Oscillator
We first consider the homogeneous equation (no external force):

$m \frac{d^2 x}{dt^2} + b \frac{dx}{dt} + k x = 0$ (Newton equation)

$\frac{d^2 x}{dt^2} + 2 \gamma \frac{dx}{dt} + \omega_0^2 x = 0$ (standard form)

where:

$\gamma = \frac{b}{2m}, \omega_0 = \sqrt{\frac{k}{m}}$

1. Physical Meaning

This equation describes a damped harmonic oscillator without external driving force.

Depending on the relation between $\gamma$ and $\omega_0$, the system behavior changes.

2. Damping Regimes

The form of the solution depends on the relation between $\gamma$ and $\omega_0$.

- Underdamping $(\gamma < \omega_0)$

Oscillatory motion with exponentially decaying amplitude.

$x(t) = e^{-\gamma t} (A\cos(\omega t)+B\sin(\omega t))$

where:

$\omega = \sqrt{\omega_0^2-\gamma^2}$

- Critical damping $(\gamma = \omega_0)$

Fastest return to equilibrium without oscillation.

$x(t) = e^{-\gamma t} (A+Bt)$

- Overdamping $(\gamma > \omega_0)$

No oscillation.

$x(t) = e^{-\gamma t} (Ae^{\sqrt{\gamma^2-\omega_0^2}t}+Be^{-\sqrt{\gamma^2-\omega_0^2}t})$

3. Non-dimensionalization

Define $\tau = \omega_0 t$, then the differential equation becomes:

$\frac{d^2 x}{d\tau^2} + 2 \zeta \frac{dx}{d\tau} + x = 0$ 

where:

$\zeta = \frac{\gamma}{\omega_0}$

This transformation makes the system's behavior depend only on the value of $\zeta$.

4. Theoretical vs Numerical solution

- Analytical solution derived from characteristic equation.
  
- Numerical solution using MATLAB ODE solver 'ode45'

- Error analysis between analytical and numerical solution.

5. Envelope analysis

For underdamping case, the amplitude envelope is:

$x_{env} (\tau) = \pm A e^{- \zeta \tau}$ 

As $\zeta$ increases, the amplitude decays faster.
