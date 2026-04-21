
## Version
branch: lmj-highRe-ptSrc
commid id: d027e0e1

## Case Description

### Plane Wave (1D problem)
#### Boundary condition

At $x=x_{min}$ face, apply:
$$
\rho(x=x_{min}, t) = \rho_{\infty} + \rho_0cos(2\pi f_0 t)
$$
with:
$\rho_0=0.01$, $f_0=0.05$ (20 steps per period)

#### Analytical solution
Reference: [Viggen 2009](https://www.researchgate.net/publication/242162544_The_Lattice_Boltzmann_Method_with_Applications_in_Acoustics) (Chapter 5.4.2 Line of point sources)

$$
\rho^\prime(x, t) = 0.38 \rho_0 e^{-\alpha_x x} e^{\mathrm{i}(2\omega_0t-kx)}
$$

with:
$$
\alpha_x = \frac{\omega_0 \tau_{vi}} {2} k_0 = \frac{\omega_0^2\tau_{vi}} {2 c_s}
$$
$$ \tau_{vi} = \frac{1}{c_s^2} (\frac{4}{3} \nu + \nu_B) $$
$$ \nu_B = \frac{2}{3} \nu $$
$$ k = \frac{k_0}{1 + \frac{3}{8}(\omega_0 \tau_{vi})^2} $$
$$ k_0 = \frac{\omega_0}{c_s} $$

### Point Source (2D problem)
#### Boundary condition

At the center of the cubic domain ($x=L_x/2$, $y=L_y/2$, $z$), apply:
$$
\rho(x, y=L/2, t)= \rho_{\infty} + \rho_0cos(2\pi f_0 t)
$$
with:
$\rho_0=0.01$, $f_0=0.05$

#### Analytical solution
Hankel function of the second kind:
$$
\rho^\prime(r, t) = A H_0^{(2)} (\hat{k}r) e^{\mathrm{i} \omega_0 t}
$$
with:
$$ A = \frac{\rho^\prime_{ref}}{|H_0^{(2)} (\hat{k} r_{ref})|} $$

$$ \hat{k} = k - \mathrm{i} \alpha_x $$

Here, we use $r_{ref}=0.0001$ and $\rho^\prime_{ref}=\rho_0$

### Pulsating Sphere (3D problem)

#### Boundary condition
At the center of a cubic domain of $300^3$, a sphere of radius $r_s=20$ pulsates with sinusoidal velocity:
$$
u(|\bold{x}|=r_s, t) = u_0 cos(2\pi f_0 t)
$$
with:
$u_0=0.01$, $f_0=0.05$

#### Analytical solution
$$
\rho^\prime(r, t) = \frac{\mathrm{i} \omega_0 \rho_0 u_0 r_s^2}{r(1+\mathrm{i} \hat{k} r_s)} e^{\mathrm{i} (\omega_0 t - \hat{k} (r-r_s))} 
$$

where
$$ \frac{1}{1+\mathrm{i} \hat{k} r_s} $$
is the frequency-dependent correction factor.
