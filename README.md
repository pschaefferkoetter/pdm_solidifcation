# pdm_solidifcation

## Description
A program that employs a strong-form meshfree collocation method for a multi-phase field model. Simulations for a low concentration (0.2 at% Sn) Al–Sn binary alloy system under periodic boundary conditions to address non-equilibrium solidification. Numerical implementation takes place through spatial discretization of the governing equations with the collocation method followed by application of the Crank–Nicolson method to integrate through time. 

## File Structure

- user needs to set geometry, material properties, and disretization params in case_1_in.m file
- case_1.sh file is used to run the program on a cluster. During exectution the **case_1.sh** will call **case_1_in.m** and **main.m** files, and generate an output folder where workspace variables can be stored for subsequent analysis and follow-on simulations.

![image](https://github.com/pschaefferkoetter/pdm_solidifcation/assets/48839550/4e699f19-7434-4ef2-9b29-1770fa1e9048)


## Program 

### case_{int}_in.m
INPUTS
| NAME | DATA TYPE | DESCRIPTION| SIZE | 
|------|:----------|:-----------|:-----|
 Length            |double |            length of box edge |
 Height            |double |         height of box edge
 discrtization_type|string  |       can be set to "uniform" or"random"
 dim               |double             spatial dimension (1D, 2D, or 3D)
 epitax_on         |boolean            used for epitaxial growth study
 num_pts_1d        |double             number of points in one Dimension
 periodic_bc_on    |boolean            sets periodic boundary conditions
 PDM_dat           |structure          data used to construct diff operat
   poly_order      |double             order of interpolation polynomial
   num_neighbor    |double             number of neighbor points
   usr_spec_diff_op array[string]     requested differential operators
 pf_dat            structure          contains following phase field data
   energy          double             sets interphase energy
   dens            double             set density of mixture
   diff_coeff      double             sets diffusion coeffecient 
   molar_volume    double             sets molar volume for the mixture
   mobility        double             sets interphase mobility
   intrf_width     double             tuning parameter for phase change
   scl             double             tuning parameter for phase change
   num_solu_phase  int                number of solutes in mixture
   pf_crit_height  double             tuning parameter for phase change
 min_spacing       double             used for random discretizations
 Reqested_Diff_Ops array[string]      requested differential operators
 seed_dat          structure          contains following seed data
   loc             cell{array[double]}user specidifed seeding position 
   radius          double             initial interphase radius
   phase           array[int]         phase of each seed 
 source_dat        structure          solute type and position for seeds
   loc             cell{array[double]}user specidifed seeding position
   val             cell{array[int]}   solute type 



 NOTE: 0 and 1 values correspond to "true" and "false" respectively 
