# Device Simulation Using C/C++ & Python Programming
Writting C/C++ and Python code in order to simulate the IV curve of a PN diode based on drift diffusion model.

## Drift Diffusion Model: Derivation
All aspects of any electronic semiconductor device relies on the dynamics of carriers in semicondutors. Many models depending on the type and scale of the device have been presented to describe carrier transport. One of the simplest and most widely used models is so called Drift-Diffusion (DD) model. According to the DD model, in a one-dimensional system, the current is given by: 

![first equation](https://latex.codecogs.com/gif.latex?%5Cdpi%7B200%7D%20j%20%3D%20n%5Cmu%20E%20-%20D%5Cfrac%7B%5Cpartial%20n%7D%7B%5Cpartial%20x%7D)

where n is the carrier concentration, μ is the carrier mobility, E is the applied electric field and D is the diffusion coefficient. 

The former DD current equation can be easily derived from the Boltzmann Transport Equation (BTE) by considering moments of the BTE. The heart of BTE is the evolution of probability density function of an ensemble of interacting particles which gives the probability of finding dN particles that all have r and p within d3rd3p.
