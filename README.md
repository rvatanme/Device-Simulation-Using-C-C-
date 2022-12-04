# Device Simulation Using C/C++ & Python Programming
Writting C/C++ and Python code in order to simulate the IV curve of a one-dimensional Si PN diode based on drift diffusion model.

## Drift Diffusion Model: Derivation
All aspects of any electronic semiconductor device relies on the dynamics of carriers in semicondutors. Many models depending on the type and scale of the device have been proposed to describe carrier transport in semiconductors. One of the simplest and most widely used models is so called Drift-Diffusion (DD) model. According to this model, in a one-dimensional system, the electron current is given by: 

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


## Discretization of Equations used in DD Model
In a complete numerical solution algorithm based on DD model to simulate a semiconductor device, the following set of equations should be solved in a one-dimensional system:


1) Poisson Equation

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cepsilon%20%5Cfrac%7Bd%5E2V%7D%7Bdx%5E2%7D%20%3D%20-%28p-n&plus;N_%7BD%7D%5E%7B&plus;%7D-N_%7BA%7D%5E%7B-%7D%29)

2) DD Current Equations

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/762996309ccc3ec7dfe1148bbafd8205759801fd)

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/81b2f3c90844ea1845c110af143a86ff0d9f3d19)

3) Continuity Equations

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7B%5Cpartial%20n%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B1%7D%7Bq%7D%5Cfrac%7BdJ_n%7D%7Bdx%7D&plus;U_n)

![](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7B%5Cpartial%20p%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B1%7D%7Bq%7D%5Cfrac%7BdJ_p%7D%7Bdx%7D&plus;U_p)

where U<sub>n</sub> and U<sub>p</sub> are the net generation-recombination rates for electrons and holes, respectively and V is the voltage at x. Using difference finite element method, the former equations are discretized in order to be numerically solved.However, there are some limitations on the choice of mesh size and time step: 1. The mesh size must be smaller than the Debye Length 2. The time step must be smaller than the dielectric relaxation time.

A simple example of Debye length effect is the redistribution of carries at an interface between two regions with different doping concentrations. Carriers diffuse into the lower doped region creating excess carrier distribution which at equilibrium decays in space down to the bulk concentration with exponential behavior. The spatial decay constant is the Debye length. 

![](https://latex.codecogs.com/svg.latex?%5Clarge%20L_D%3D%5Csqrt%7B%5Cfrac%7B%5Cepsilon%20k_BT%7D%7Bq%5E2N%7D%7D)

In GaAs and Si, at room temperature the Debye length is approximately 400 Å when N ≈ 1E-16 /cm<sup>3</sup> and decreases to about only 50 Å when N ≈ 1E-18 /cm<sup>3</sup>. 

The dielectric relaxation time is the characteristic time for charge fluctuations to decay under the influence of the field that they produce. The dielectric relaxation time may be estimated using:

![](https://latex.codecogs.com/svg.latex?%5Clarge%20t_%7Bdr%7D%3D%5Cfrac%7B%5Cepsilon%7D%7BqN%5Cmu%7D)

In order to numerically solve the set of mentioned Poisson, DD current and continuity equations, each of them should be discretized first. The Poisson equation can be discritized by finit difference element as following (assuming φ<sup>new</sup> = φ<sup>old</sup>+δ):

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20%5Cfrac%7B1%7D%7Bdx%5E2%7D%5Cphi_%7Bi&plus;1%7D%5E%7Bn&plus;1%7D&plus;%28%5Cfrac%7B2%7D%7Bdx%5E2%7D&plus;n_i&plus;p_i%29%5Cphi_i%5E%7Bn&plus;1%7D&plus;%5Cfrac%7B1%7D%7Bdx%5E2%7D%5Cphi_%7Bi-1%7D%5E%7Bn&plus;1%7D%3D%20%5C%5C%20%5C%5C%20-%28p_i-n_i&plus;N_%7Bdi%7D-N_%7Bai%7D%29-%28p_i&plus;n_i%29%5Cphi_i%5En)

where φ<sub>i</sub><sup>n</sup> is the potential of ith node of the discritized device in the nth iteration. The n<sub>i</sub> and p<sub>i</sub> are donor and acceptor doping concentrations in the ith node. The n<sub>i</sub> and p<sub>i</sub> are the electron and hole densities in the ith node. At equilibrium, ion densities can be calculated by boltzman equation as follows:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20n_i%20%3D%20n_0e%5E%7B-%5Cfrac%7Bq%5Cphi_i%7D%7BkT%7D%7D%2C%20%5C%3B%5C%3B%5C%3B%5C%3B%20p_i%20%3D%20p_0e%5E%7B-%5Cfrac%7Bq%5Cphi_i%7D%7BkT%7D%7D)

where n<sub>0</sub> and p<sub>0</sub> are electron and hole densities in pure bulk phase at room temperature. The -qφ<sub>i</sub> is the required work for moving an electron or hole from infinity to the ith node. By substituting Boltzman equation in the poisson equation, the so called Poisson–Boltzmann equation is obtained which should be solved in some part of any device simulation algorithm (usually called as Poisson–Boltzmann solver). The Poisson–Boltzmann equation can take many forms throughout various scientific fields. In biophysics and certain surface chemistry applications, it is known simply as the Poisson–Boltzmann equation. It is also known in electrochemistry as Gouy-Chapman theory; in solution chemistry as Debye–Huckel theory; in colloid chemistry as Derjaguin–Landau–Verwey–Overbeek (DLVO) theory. Only minor modifications are necessary to apply the Poisson–Boltzmann equation to various interfacial models, making it a highly useful tool in determining electrostatic potential at surfaces or interfaces.

The discretization of the DD current equations in a coservative form (under continuity equation condition) requires the knowledge of current densities on the mid-points of the mesh lines connecting nieghboring grid nodes. Since solution are available only on the grid node, interpolation schemes are needed to determine the solutions. There are two schemes that one can use: (1) Linearized scheme where V, n, p, mu and D vary linearly between nieghboring grid points (2) Scharfetter-Gummel scheme where electron and hole densities follow exponential variation between grid points. Since first scheme can lead to substantial errors in the regions with high electric field and high doping concentration, only discretizatioin under second scheme would be explained. One can solve the current density equation as follows:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20J_%7Bi&plus;1/2%7D%20%3D%20-en%5Cmu%20_%7Bi&plus;1/2%7D%5Cfrac%7BV_%7Bi&plus;1%7D-V_i%7D%7Ba_i%7D&plus;eD_%7Bi&plus;1/2%7D%5Cfrac%7B%5Cpartial%20n%7D%7B%5Cpartial%20x%7D%5C%5C%5C%5C%20%3D-en%5Cmu%20_%7Bi&plus;1/2%7D%5Cfrac%7BV_%7Bi&plus;1%7D-V_i%7D%7Ba_i%7D&plus;eD_%7Bi&plus;1/2%7D%5Cfrac%7B%5Cpartial%20n%7D%7B%5Cpartial%20V%7D%5Cfrac%7B%5Cpartial%20V%7D%7B%5Cpartial%20x%7D)

where V<sub>i</sub> is the potential in the ith node. Since for n(V<sub>i</sub>) subject to the boundary conditions is n<sub>i</sub>, therefore the solution of this first order differential equation leads to the following discretized equation (so called Sharfetter-Gummel Discretized Continuity Equations):

![](https://github.com/rvatanme/Device-Simulation-Using-C-C-/blob/main/Continuity_Equ.png)

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20B%28x%29%3D%5Cfrac%7Bx%7D%7Be%5Ex&plus;1%7D%20%5C%3B%5C%3B%5C%3B%5C%3BBernoulli%5C%3BFunction)

where the dependency of n and p upon the potential is buried under Bernoulli function. Finaly, after the convergence, the electron current density from node (i) to (i+1) can be calculated by the following equation:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20J_%7Bi&plus;1/2%7D%20%3D%20%5Cfrac%7BeD%20_%7Bi&plus;1/2%7D%7D%7Ba_i%7D%5Bn_%7Bi&plus;1%7DB%28%5Cfrac%7BV_%7Bi&plus;1%7D-V_i%7D%7BV_t%7D%29-n_iB%28%5Cfrac%7BV_i-V_%7Bi&plus;1%7D%7D%7BV_t%7D%29%5D)

## Numerical Sulotion of Discretized Equations
The result of the discretization using finite difference, is several coupled of linear algebraic equations Au=f. Each of these linear algebraic equations can be solved with the variety of methods including: (1) Direct methods a: Guassian elemination b: LU decomposition method (2) Iterative method a: mesh relaxation methods (Successive Over-Relaxation (SOR), or LSOR) b: Matrix methods (Jacobi, Gauss-Seidel, ...). Here in the code, the LU decomposition method is used in order to solve each linear equation. For a given tridiagonal matrix Au=f (where A and f is known and the goal is finding u), A can be decomposed to L and U as follows: 

![](https://github.com/rvatanme/Device-Simulation-Using-C-C-/blob/main/LU_M.png)

where the L and U matrix are given by:

![](https://github.com/rvatanme/Device-Simulation-Using-C-C-/blob/main/LU.png)

The solution then would be as follows:

![](https://github.com/rvatanme/Device-Simulation-Using-C-C-/blob/main/LY.png)

![](https://github.com/rvatanme/Device-Simulation-Using-C-C-/blob/main/Uu.png)

Since, in DD model algorithm, there are several coupling linear algebraic equations, an iterative method is also needed in order to solve these equations together. The two popular methods for solving the coupled discretized equations are the Gummel's iteration method<sup>1</sup> (used in the code) and the Newton's method. 

Gummel's  method  solves  the  coupled  set  of  semiconductor  equations  together  with  the  Poisson equation via a decoupled procedure. If we choose the quasi-Fermi level formulation, we solve first a nonlinear Poisson's equation. The potential obtained from this solution is substituted into  the  continuity  equations,  which  are  now  linear,  and  are  solved  directly  to  conclude  the  iteration  step.  The  result  in  terms  of  quasi-Fermi  levels  is  then  substituted  back  into  Poisson's  equation   and   the   process   repeated   until   convergence   is   reached.   In   order   to   check   for   convergence, one can calculate the residuals obtained by positioning all the terms to the left hand side  of  the  equations  and  substituting  the  variables  with  the  iteration  values.  For  the  exact  solution  the  residuals  should  be  zero.  Convergence  is  assumed  when  the  residuals  are  smaller  than a set tolerance. The rate of convergence of the Gummel method is faster when there is little coupling between the different equations. The computational cost of one Gummel iteration is one matrix solution for each carrier type plus one iterative solution for the linearization of Poisson's equation.  Note  that  in  conditions  of  equilibrium  (zero  bias)  only  the  solution  of  Poisson's  equation  is  necessary,  since  the  equilibrium  Fermi  level  is  constant  and  coincides  with  both  quasi-Fermi levels.

## Simulation Procedure in Detailes
The puprpose of this simulation is to predict the current flowing trhough a 1d pn diode under an applied bias. This is achieved by approximating the operation of the device onto a one dimensional grid, consisting of a number of grid points called nodes. By applying the set of finit element differential equations (derived above from DD, Poisson's  and continuity equations) onto this grid one can simulate the transport of carriers through the mentioned structure. The detailed numerical procedure that must be impelemented is as following:
 
1) Solve only Poisson's equation at equilibrium (no applied bias).

   a) Choose an initial guess for the potential at each grid point for example φ<sub>i</sub> = (kT/q)ln(Na/n<sub>0</sub>) for P side and
   φ<sub>i</sub> = -(kT/q)ln(Nd/n<sub>0</sub>) for N side where φ<sub>i</sub> and n<sub>0</sub> is the potential in the ith grid point and
   intrinsic carrier concentration, respectively.
   
   b) Write the potential at the next iteration step (n+1) as φ<sup>n+1</sup> = φ<sup>n</sup> + δφ. The new potential is solved using the 
   following discritized Poisson's equations derived above,  
   
   ![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20%5Cfrac%7B1%7D%7Bdx%5E2%7D%5Cphi_%7Bi&plus;1%7D%5E%7Bn&plus;1%7D&plus;%28%5Cfrac%7B2%7D%7Bdx%5E2%7D&plus;n_i&plus;p_i%29%5Cphi_i%5E%7Bn&plus;1%7D&plus;%5Cfrac%7B1%7D%7Bdx%5E2%7D%5Cphi_%7Bi-1%7D%5E%7Bn&plus;1%7D%3D%20%5C%5C%20%5C%5C%20-%28p_i-n_i&plus;N_%7Bdi%7D-N_%7Bai%7D%29-%28p_i&plus;n_i%29%5Cphi_i%5En)
   

   ![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20n_i%20%3D%20n_0e%5E%7B-%5Cfrac%7Bq%5Cphi_i%7D%7BkT%7D%7D%2C%20%5C%3B%5C%3B%5C%3B%5C%3B%20p_i%20%3D%20p_0e%5E%7B-%5Cfrac%7Bq%5Cphi_i%7D%7BkT%7D%7D)

   The two above equatons result in to several coupled of linear algebraic equations in the form of Au = f where A, u and f are coefficints
   matrix, unknow φ<sup>n+1</sup>, and known quantities obtained from φ<sup>n</sup>. Since A is a tridiagonal matrix that can be solved
   using LU decomposision explained above.


## References
1) Drift-Diffusion Model: Solution Details, Prof. Dragica Vasileska, Arizona State University (https://nanohub.org/resources/1565/download/ddmodel_solution_details_word.pdf)
