# Device Simulation Using C/C++ & Python Programming
Writting C/C++ and Python code in order to simulate the IV curve of a PN diode based on drift diffusion model.

## Drift Diffusion Model: Derivation
All aspects of any electronic semiconductor device relies on the dynamics of carriers in semicondutors. Many models depending on the type and scale of the device have been presented to describe carrier transport. One of the simplest and most widely used models is so called Drift-Diffusion (DD) model. According to the DD model, in a one-dimensional system, the electron current is given by: 

![first equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/762996309ccc3ec7dfe1148bbafd8205759801fd)

where e is electron charge n is electron concentration, μ<sub>e</sub> is electron mobility, E is the applied electric field and D<sub>e</sub> is electron diffusion coefficient. The hole current is given by a similar equation.

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/81b2f3c90844ea1845c110af143a86ff0d9f3d19)

The former DD current equations can be easily derived from the Boltzmann Transport Equation (BTE) by considering moments of the BTE. The heart of BTE is the evolution of probability density function of an ensemble of interacting particles which gives the probability of finding dN particles that all have r and p within d<sup>3</sup>rd<sup>3</sup>p.

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/2a24c0ab57c3bcfaa078c66b51329b9d68183c5d)

where r and p are the location and momentum of the particles. The set of all possible positions r and momenta p is called the phase space of the system. The general BTE equation then can be written based on force, diffusion and collision terms. 

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/5d6c3edd68a8d09e67bd377bd8906d144e1a911f)

where the "force" term corresponds to the forces exerted on the particles by an external influence (not by the particles themselves), the "diff" term represents the diffusion of particles, and "coll" is the collision term – accounting for the forces acting between particles in collisions. If we assume that F(r, t) is a force field acting on the particles, and m is the mass of the particles, the the final statement of BTE would be:
