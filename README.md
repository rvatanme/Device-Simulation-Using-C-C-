# Device Simulation Using C/C++ & Python Programming
Writting C/C++ and Python code in order to simulate the IV curve of a PN diode based on drift diffusion model.

## Drift Diffusion Model: Derivation
All aspects of any electronic semiconductor device relies on the dynamics of carriers in semicondutors. Many models depending on the type and scale of the device have been presented to describe carrier transport. One of the simplest and most widely used models is so called Drift-Diffusion (DD) model. According to the DD model, in a one-dimensional system, the electron current is given by: 

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/762996309ccc3ec7dfe1148bbafd8205759801fd) (1)

where e, n, μ<sub>e</sub>, E and D<sub>e</sub> are electron charge, electron concentration, electron mobility, applied electric field and electron diffusion coefficient at x, respectively. The hole current is given by a similar equation.

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/81b2f3c90844ea1845c110af143a86ff0d9f3d19) (2)

The former DD current equations can be easily derived from the Boltzmann Transport Equation (BTE) by considering moments of the BTE. The heart of BTE is the evolution of probability density function of an ensemble of interacting particles which gives the probability of finding dN particles that all have r and p within d<sup>3</sup>rd<sup>3</sup>p.

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/2a24c0ab57c3bcfaa078c66b51329b9d68183c5d) (3)

where r and p are the location and momentum of the particles. The set of all possible positions r and momenta p is called the phase space of the system. The general BTE equation then can be written based on force, diffusion and collision terms. 

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/5d6c3edd68a8d09e67bd377bd8906d144e1a911f) (4)

where the "force" term corresponds to the forces exerted on the particles by an external influence (not by the particles themselves), the "diff" term represents the diffusion of particles, and "coll" is the collision term – accounting for the forces acting between particles in collisions. If we assume that F(r, t) is a force field acting on the particles, and m is the mass of the particles, the the final statement of BTE would be:

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/fdf75926860831e97d1bd746278c82554fc10640) (5)

The collision less BTE, where individual collisions are replaced with long-range aggregated interactions, e.g. Coulomb interactions, is often called the Vlasov equation. The final statement of BTE is more useful than the principal one above, but it is still incomplete, since f cannot be solved unless the collision term in f is known. This term cannot be found as easily or generally as the others – it is a statistical term representing the particle collisions, and requires knowledge of the statistics the particles obey, like the Maxwell–Boltzmann, Fermi–Dirac or Bose–Einstein distributions. At time t = 0, we will switch of all external forces of the system. By scattering the system reaches the thermodynamic equilibrium state again, which will be described in a linear approximation by:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%28%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20t%7D%29%20%3D%20%28%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20t%7D%29_%7Bscat%7D%20%3D%20%5Cfrac%7Bf%28%5Cvec%7Br%7D%2C%5Cvec%7Bk%7D%2Ct%29-f_0%28%5Cvec%7Br%7D%2C%5Cvec%7Bk%7D%29%7D%7B%5Ctau%28%5Cvec%7Bk%7D%29%7D) (6)

The f<sub>0</sub>(r , k) is the distribution function at the equilibrium. The relaxation time τ(k) describes how fast the system reaches thermodynamic equilibrium again. With the use of a relaxation time approximation, the Boltzmann transport equation can be written as:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7BeE%7D%7Bm%5E*%7D%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20%5Cnu%7D%20&plus;%20%5Cnu%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20x%7D%20%3D%20%5Cfrac%7Bf_0-f%28%5Cnu%2Cx%29%7D%7B%5Ctau%7D) (7)

Here, m<sup>*</sup> is the effective mass of the given carrier. The integral on the right hand side represents the first 'moment' of the distribution function. After multiplying both sides by ν and integrating over ν, from the RHS of the equation we get:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7B1%7D%7B%5Ctau%7D%5B%5Cint%20%5Cnu%20f_0d%5Cnu%20-%20%5Cint%20%5Cnu%20f%28%5Cnu%2Cx%29d%5Cnu%5D%20%3D%20-%5Cfrac%7BJ%28x%29%7D%7Be%5Ctau%7D) (8)

![](https://latex.codecogs.com/svg.latex?%5Clarge%20J%28x%29%20%3D%20e%20%5Cint%20%5Cnu%20f%28%5Cnu%2Cx%29d%5Cnu)

The distribution function is symmetric at equilibrium in ν, and hence the first integral in equaution 8 is zero. By considering LHS of equation 7, we have:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20J%28x%29%20%3D%20-e%20%5Cfrac%7Be%5Ctau%7D%7Bm%5E*%7DE%5Cint%20%5Cnu%20%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20%5Cnu%7D%20d%5Cnu%20-e%5Ctau%20%5Cfrac%7Bd%7D%7Bdx%7D%20%5Cint%20%7B%5Cnu%7D%5E2%20f%28%5Cnu%2Cx%29%20d%5Cnu) (9)

By plugging the following equations in the equation 9

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cint%20%5Cnu%20%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20%5Cnu%7D%20d%5Cnu%20%3D%20%5B%5Cnu%20f%28%5Cnu%2Cx%29%5D_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7D-%5Cint%20f%28%5Cnu%2Cx%29d%5Cnu%20%3D%20-n%28x%29)

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cint%20%7B%5Cnu%7D%5E2f%28%5Cnu%2Cx%29d%5Cnu%20%3D%20n%28x%29%3C%7B%5Cnu%7D%5E2%3E)

we get the DD current equation as presented in the beginning of the discussion: 

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/762996309ccc3ec7dfe1148bbafd8205759801fd)

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/81b2f3c90844ea1845c110af143a86ff0d9f3d19)

where the carrier mobility, diffusion and average squared velocity are given by:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cmu%20%3D%20%5Cfrac%7Be%5Ctau%7D%7Bm*%7D)

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%3C%7B%5Cnu%7D%5E2%3E%20%3D%20%5Cfrac%7Bk_BT%7D%7Bm%5E*%7D)

![](https://latex.codecogs.com/svg.latex?%5Clarge%20D%20%3D%20%5Cmu%20k_BT)

Using thermal equilibrium velocity in driving the DD current equations means that the DD model is only valid for small perturbations of the equilibrium state that holds at low electric fields. The validity of the DD equations is empirically extended by introducing field-dependent mobility μ(E) and diffusion coefficient D(E), obtained from empirical models or detailed calculation to capture effects such as velocity saturation at high electric fields due to hot carrier effects. 


## Numerical Solution Based on DD Model
In a complete numerical solution algorithm based on DD model to simulate a semiconductor device, the following set of equations are employed in one-dimension:

1) DD Current Equations

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/762996309ccc3ec7dfe1148bbafd8205759801fd)

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/81b2f3c90844ea1845c110af143a86ff0d9f3d19)

2) Continuity Equations

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7B%5Cpartial%20n%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B1%7D%7Bq%7D%5Cfrac%7BdJ_n%7D%7Bdx%7D&plus;U_n)

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7B%5Cpartial%20p%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B1%7D%7Bq%7D%5Cfrac%7BdJ_p%7D%7Bdx%7D&plus;U_p)

3) Poisson Equation

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cepsilon%20%5Cfrac%7Bd%5E2V%7D%7Bdx%5E2%7D%20%3D%20-%28p-n&plus;N_%7BD%7D%5E%7B&plus;%7D-N_%7BA%7D%5E%7B-%7D%29)


where U<sub>n</sub> and U<sub>p</sub> are the net generation-recombination rates for electrons and holes, respectively and V is the voltage at x. Using difference finite element method, the former equations are discretized in order to be numerically solved.However, there are some limitations on the choice of mesh size and time step: 1. The mesh size must be smaller than the Debye Length 2. The time step must be smaller than the dielectric relaxation time.

A simple example of Debye length effect is the redistribution of carries at an interface between two regions with different doping concentrations. Carriers diffuse into the lower doped region creating excess carrier distribution which at equilibrium decays in space down to the bulk concentration with exponential behavior. The spatial decay constant is the Debye length. 

![](https://latex.codecogs.com/svg.latex?%5Clarge%20L_D%3D%5Csqrt%7B%5Cfrac%7B%5Cepsilon%20k_BT%7D%7Bq%5E2N%7D%7D)

In GaAs and Si, at room temperature the Debye length is approximately 400 Å when N ≈ 1E-16 /cm<sup>3</sup> and decreases to about only 50 Å when N ≈ 1E-18 /cm<sup>3</sup>. 

The dielectric relaxation time is the characteristic time for charge fluctuations to decay under the influence of the field that they produce. The dielectric relaxation time may be estimated using:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20t_%7Bdr%7D%3D%5Cfrac%7B%5Cepsilon%7D%7BqN%5Cmu%7D)

The two popular methods for solving the discretized equations are the Gummel's iteration method and the Newton's method. It is common practice to perform the actual calculation using normalized units to make the algorithms more efficient, and in cases to avoid numerical overflow and underflow. It is advisable to input the data in M.K.S. or practical units (the use of centimeters is for instance very common in semiconductor practice, instead of meters) and then provide a conversion block before and after the computation blocks to normalize and denormalize the variables.
