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


## Driven Damped Oscillator
We consider the equation (with external driving force):

$m \frac{d^2 x}{dt^2} + b \frac{dx}{dt} + k x = F_0 cos(\omega t)$ (Newton equation)

$\frac{d^2 x}{dt^2} + 2 \gamma \frac{dx}{dt} + \omega_0^2 x =  F cos(\omega t)$ (standard form)

where: 

$\gamma = \frac{b}{2m}, \omega_0 = \sqrt{\frac{k}{m}}, F = \frac{F_0}{m}$

1. Physical Meaning

This equation describes a driven damped oscillator.

After the transient response decays, the system reaches a steady-state motion at the driving frequency.

2. Solution analysis (steady-state)
- Theoretical solution (steady-state)


$Let \ A = \frac{F}{\sqrt{(\omega_0^2-\omega^2)^2 + 4 \gamma^2 \omega^2}} \ , \phi \ = arctan(\frac{2 \gamma \omega}{\omega_0^2-\omega^2})$

$x(t) \ = A cos (\omega t - \phi)$

- Numerical solution
  
Solved using ODE solver 'ode45'

Compared with theoretical steady-state solution

Error analysis performed
3. Steady-state Extraction using FFT (Fast Fourier Transform)
- Separation of transient and steady-state

Amplitude extraction via FFT

Verification of dominant frequency

4. Parameter sweep $(\gamma, \omega)$
- Amplitude vs. driving frequency $\omega$ sweep

Amplitude has maxima when $\omega$ is near $\omega_0$

$\omega_{res} \ = \sqrt{\omega_0^2-2 \gamma^2}$

- Amplitude vs. damping coefficient $\gamma$ sweep

Amplitude decreases as $\gamma$ increases

5. Resonance analysis
- Resonance frequency

$\omega_{res} \ = \sqrt{\omega_0^2-2 \gamma^2}$


$\gamma < \frac{\omega_0}{\sqrt{2}}$

6. Phase analysis
- Phase vs. driving frequency

$\omega >> \omega_0$ , $\phi$ -> $\pi$

$\omega << \omega_0$ , $\phi$ -> $0$

$\omega \approx \omega_0$ , $\phi$ -> $\frac{\pi}{2}$
