# pdm_solidifcation

## Description
A program that employs a strong-form meshfree collocation method for a multi-phase field model. Simulations for a low concentration (0.2 at% Sn) Al–Sn binary alloy system under periodic boundary conditions to address non-equilibrium solidification. Numerical implementation takes place through spatial discretization of the governing equations with the collocation method followed by application of the Crank–Nicolson method to integrate through time. 

## File Structure

- user needs to set geometry, material properties, and disretization params in case_1_in.m file
- case_1.sh file is used to run the program on a cluster. During exectution the **case_1.sh** will call **case_1_in.m** and **main.m** files, and generate an output folder where workspace variables can be stored for subsequent analysis and follow-on simulations.

---
![image](https://github.com/pschaefferkoetter/pdm_solidifcation/assets/48839550/498c3895-2a60-4755-89b7-b4a924581a6d)

---

### Input Parameters

| NAME | DATA TYPE | DESCRIPTION|
|------|:----------|:-----------|
 Length            |double              |length of box edge                 |
 Height            |double              |height of box edge                 |
 discrtization_type|string              |can be set to "uniform" or"random" |
 dim               |double              |spatial dimension (1D, 2D, or 3D)  |
 epitax_on         |boolean             |used for epitaxial growth study    |
 num_pts_1d        |double              |number of points in one Dimension  |
 periodic_bc_on    |boolean             |sets periodic boundary conditions  |
 PDM_dat           |structure           |data used to construct diff operat |
   poly_order      |double              |order of interpolation polynomial  |
   num_neighbor    |double              |number of neighbor points          |
   usr_spec_diff_op|array[string]       |requested differential operators   |
   pf_dat          |structure           |contains following phase field data|
   energy          |double              |sets interphase energy             |
   dens            |double              |set density of mixture             |
   diff_coeff      |double              |sets diffusion coeffecient         |
   molar_volume    |double              |sets molar volume for the mixture  |
   mobility        |double              |sets interphase mobility           |
   intrf_width     |double              |tuning parameter for phase change  |
   scl             |double              |tuning parameter for phase change  |
   num_solu_phase  |int                 |number of solutes in mixture       |
   pf_crit_height  |double              |tuning parameter for phase change  |
 min_spacing       |double              |used for random discretizations    |
 Reqested_Diff_Ops |array[string]       |requested differential operators   |
 seed_dat          |structure           |contains following seed data       |
   loc             |cell{array[double]} |user specidifed seeding position   |
   radius          |double              |initial interphase radius          |
   phase           |array[int]          |phase of each seed                 |
 source_dat        |structure           |solute type and position for seeds |
   loc             |cell{array[double]} |user specidifed seeding position   |
   val             |cell{array[int]}    |solute type                        |

# Main Program Files
---
## Al_Sn_Gibbs.m

### Description:
Thermodynamic properties of pure Elements. Note: the 1st and 2nd derivates of these functions are related to the absolute entropy (S) and heat capacity of the compound at the same temperature. Data taken from COST-507 report

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
T | double | temperature | -
c | double | concentration | -
phase | int | phase id  | -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
f_gibbs | doube | gibbs energy | -
df_gibbs | double | derivative of f_gibbs | -

### Called by
Evolve_fields_scl_calc_rev2.m
---
## box_geometry_random.m
### Description:
Generates uniform descretization for a 2D domain square domain

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 min_spacing    |double           |min spacing between pts  |-
 Length         |double           |domain length            |-
 dim            |int              |spatial dimension        |-
### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord     |     array[double]  |  coordinates            |  [dim X num_pts]
 num_pts   |     int            |  number of PDM points   |  -
 min_achieved_dist | double | min dist between two points | -

### Called by
main1.m
----
---
## box_geometry_uniform_v2.m
### Description:
Generates uniform descretization for a 2D domain

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 min_spacing    |double           |min spacing between pts  |-
 Length         |double           |domain length            |-
 Height         |double           |domain height            |-
 dim            |int              |spatial dimension        |-

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord     |     array[double]  |  coordinates            |  [dim X num_pts]
 num_pts   |     int            |  number of PDM points   |  -

### Called by
 ./main<int>_epitax.m
---
## Capillary_term_scl
### Description: Calcuates cappilary term

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 inode       |     int           |  node id            |            -
 alpha       |     int           | phase id            |            -
 phi_raw     |     array[double] |order parameter      |   [num_pts x num_pts]
 Lap_op_Mat  |    array[double]  | numerical laplacian |   [num_pts x num_pts]
   eta       |     double        |  interface width    |            -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
  I_alpha    |      double   |       cappillary term   |             -

### Called by
Evolve_fields_scl_calc_rev2.m

---

## Evolve_fields_scl_calc_rev2.m
### Description:
Calculates phase and concentration field velocities for a point in the domain

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
Lap_Op_mat    | array[double]     |  laplace diff op               |[num_pts x num_pts]
D1x_mat       | array[double]     |  first derivatitive, X op      |[num_pts x num_pts]
D1y_mat       | array[double]     |  first derivatitive, X op      |[num_pts x num_pts]
conc_raw      | array[double]     |  uninterpoloated concentration |[num_pts X num_phi]
conc          | array[double]     |  interpoloated conecntration   |[num_pts X num_phi]
PF_dat        | structure         |  see case{int}_in.m            |   -
Temp          | double            |  Temperature                   |[num_pts X 1]
TOL           | double            |  Tolerance                     |    -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 phase_vel_pnt_i double  |           phase velocity at pointi   |       -
 conc_vel_pnt_i double   |           phase velocity at pointi   |

### Called by
MPF.m

---

## calc_diff_ops.m
### Description:
Calculates differential Operators

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
coord            |array[double]    |coordinates                       |[dim X num_pts]
dialation_parm   |array[int]       |dialation param                   |[1 X num_pts]
dist_list        |cell{array[int]} |distance to nbr                   |{num_pts[num_nbrs]}
nblist           |cell{array[int]} |list of nbrs                      |{num_pts[num_nbrs]}
poly_order       |double           |polynomial order                  | -
num_nbr          |int              |number of neighbors               | -
usr_spec_diff_op | array[string]   |requested differential operators  | _

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 PDM_diff_operators | cell{array}  | differnetial operators | {num_pts[num_ops X num+pts]}

### Called by
main.m
---
## Dx_mat.m
### Description:
Calculates first derviative numerical operators Dx.

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord              |array[double]  | coordinates                    | [dim X num_pts]
 PDM_diff_operators |cell{array}    | differnetial operators         | {num_pts[num_ops X num+pts]}

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
Lap_Op_Mat      |  array[double] |   numerical laplace op         |   [num_pts x num_pts]

### Called by
MPF.m

---
## Dy_mat.m
### Description:
Calculates first derviative numerical operators Dy.

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord              |array[double]  | coordinates                    | [dim X num_pts]
 PDM_diff_operators |cell{array}    | differnetial operators         | {num_pts[num_ops X num+pts]}

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
Lap_Op_Mat      |  array[double] |   numerical laplace op         |   [num_pts x num_pts]

### Called by
MPF.m

---
## dialation.m
### Description:
  Calculates the dialation parameters used to construct differential operators for each discrete point in the domain
### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord          |array[double]    |coordinates           |[dim X num_pts]
 exp_coord      |array[double]    |used for periodic BCs |size(coord)
 nblist         |cell{array[int]} |list of nbrs          | {num_pts[num_nbrs]}
 num_nbr        |int              |number of neighbors   | -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
nblist         |cell{array[int]} |list of nbrs       |{num_pts[num_nbrs]}
exp_dist_list  |cell{array[int]} |distance to nbr    |{num_pts[size(coord]}
dialation_parm |array[int]       |dialation param    [1 X num_pts]

### Called by
get_neighbors_gen.m
---
## diffus_alpha
### Description: Calculates diffusion term for phase alpha

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 inode       |     int         |    node id               |         -
 alpha       |     int         |   phase id               |         -
 phi         |   array[double] |order parameter           |  [num_pts x num_pts]
 conc        |   array[double] |concentration             | [num_pts x num_pts]
 Lap_op_Mat  |   array[double] |  numerical laplacian     | [num_pts x num_pts]
 D1x         |   array[double] |  numerical Derivative y  |  [num_pts x num_pts]
 D1y         |   array[double] |   numerical Derivative x |  [num_pts x num_pts]

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 diffus_term   |       double      |  diffusion terms    |            -
### Called by
Evolve_fields_scl_calc_rev2

---

## GBCT_Al_bct_A5
### Description:
Calculates stable state enthalpy for soid Aluminum bct phase based on COST509

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 T      |          double        |      temperature (K)     |         -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 G_Al_bct    |  double     |      Stable Gibbes Energy       |     -

### Called by
Al_Sn_Gibbs

---
## GFCC_Sn_fcc_A1
### Description:
Calculates stable state enthalpy for solid Sn fcc phase based on COST509

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 T      |          double        |      temperature (K)     |         -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 G_Sn_bct_A5    |  double     |      Stable Gibbes Energy     |       -

### Called by
Al_Sn_Gibbs

---

## GHSR_Sn_fcc_A1
### Description:
Calculates stable state enthalpy for solid Aluminum fcc phase based on COST509

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 T      |          double        |      temperature (K)     |         -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 G_Al_fcc_A1   |  double     |      Stable Gibbes Energy      |      -

### Called by
Al_Sn_Gibbs


---
## GHSR_Sn_bct_A5
### Description:
Calculaes stable state enthalpy for Sn_BCT_A phase based on COST509

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 T      |          double        |      temperature (K)     |         -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 G_Sn_bct_A5    |  double     |      Stable Gibbes Energy     |       -

### Called by
Al_Sn_Gibbs

---
## GLIQ_Sn_liq
### Description:
Calculates Stable state enthalpy for sin liquid phase based on COST509
### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 T      |          double        |      temperature (K)     |         -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 G_Sn_bct_liq    |  double     |      Stable Gibbes Energy   |        -

### Called by
Al_Sn_Gibbs

---

## GLIQ_Al_liq
### Description:
Calculates stable state enthalpy for liquid Aluminum phase based on COST509
### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 T      |          double        |      temperature (K)     |         -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
G_Al_liq     |  double     |      Stable Gibbes Energy      |      -

### Called by
Al_Sn_Gibbs

---
## intrf_dvg_force_v1
### Description: Calculates interphase driving force

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
  inode    |        int        |     node id                    |    -
  alpha    |        int        |    phase alpha id              |    -
  beta     |        int        |    phase beta id               |    -
  f_alpha  |       double      |   gibbs energy phase alpha     |    -
  f_beta   |       double      |   gibbs energy phase beta      |    -
  df_alpha |       double      |    derivative of gibbs energy  |    -
  df_beta  |       double      |    derivative of gibbs energy  |    -
  phi      |    array[double]  |   order parameter              |    [num_pts x num_pts]
  n_c      |         int       |     number of components       |    -
  eta    |        double     |     interface width            |    -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 Delta_g     |    double   |       interphase driving force   |     -

### Called by
Evolve_fields_scl_calc_rev2

---

## get_neighbors_gen.m
### Description:
Generates neighbors and modifies coordinates to handle periodic BCs.

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 Variable Name  |Data Type        |Description         | size
 Length         |double           |domain length       | -
 Height         |double           |domain height       | -
 coord          |array[double]    |coordinates         |[dim X num_pts]
 periodic_bc_on |bool             |for periodic BCs    | -
 num_nbr        |int              |number of neighbors | -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 nblist         |cell{array[int]} |list of nbrs       |{num_pts[num_nbrs]}
 dist_list      |cell{array[int]} |distance to nbr    |{num_pts[num_nbrs]}
 exp_dist_list  |cell{array[int]} |distance to nbr    |{num_pts[size(coord]}
 dialation_parm |array[int]       |dialation param    |[1 X num_pts]
 coord          |array[double]    |modified coord     |[size(coord)]

### Called by
 ./main<int>_epitax.m

### Dependencies
- periodic_coordinates_v2
- KDTreeSearcher  (Requieres Matlab ML and Statistics toolbox)
- knnsearch       (Requieres Matlab ML and Statistics toolbox)
- dialation

----
## interp_mat.m
### Description:
Calculates C0 interpoloation operator.

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord              |array[double]  | coordinates                    | [dim X num_pts]
 PDM_diff_operators |cell{array}    | differnetial operators         | {num_pts[num_ops X num+pts]}

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
Lap_Op_Mat      |  array[double] |   numerical laplace op         |   [num_pts x num_pts]

### Called by
MPF.m


---------
## laplace_mat.m
### Description:
Calculates numerical laplace operator.

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord              |array[double]  | coordinates                    | [dim X num_pts]
 PDM_diff_operators |cell{array}    | differnetial operators         | {num_pts[num_ops X num+pts]}

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
Lap_Op_Mat      |  array[double] |   numerical laplace op         |   [num_pts x num_pts]

### Called by

## Depedencies
---
## main.m
### Description:
  - discretizes geometry, constructs numerical differential operators, and exectutes simulation
  - Called by ../input_folder/case_{int}.m

### Dependencies
  - box_geometry_uniform_v2.m
  - get_neighbors_gen.m
  - calc_diff_ops.m
---------

## MPF.m
### Description:
Main Engine for solidification simulation.
### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 restart       | bool               | used for epitaxial sim           |-
 coord         |  array[double]     | coordinates                      |[dim X num_pts]
 PF_dat        | structure          |see case{int}_in.m                |-
 seed_dat      | structure          |see case{int}_in.m                |-
 nblist        | cell{array[int]}   |list of nbrs                      |{num_pts[num_nbrs]}
 n_step        | int                | number of steps in sim           | -
 radius        | double             |interphase radius                 |-
 usr_spec_diff_op | array[string]   | requested differential operators | -
 dt            | double             | time increment                   | -
 dist_list   | cell{array[double]}  | point to neighbor list           | {num_pts[num_nbrs]}
 vec_dist_list | cell{array[double]}| vector list                      | {num_pts[dim x num_nbrs]}]

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
STORE  |.mat|phase values for each field | -

### Called by
main{int}.m

## Depedencies
- MPF_Inititialization_v4
- laplace_mat
- interp_mat
- Dx_mat
- Dy_mat
- Evolve_fields_scl_calc_rev2

---

## MPF_Initialization.m
### Description:
Iniitializes order parameter array

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 Variable Name   |Data Type        |Description               |     size
 coord           |array[double]    |coordinates               |   [dim X num_pts]
 PF_dat          |structure        |see case{int}_in.m        |       -
 seed_dat        |structure        |see case{int}_in.m        |       -
 nblist          |cell{array[int]} |list of nbrs              | {num_pts[num_nbrs]}
 Temp            |double           |Temperature               |       -
 TOL             |double           |tolerance                 |       -
 init_field_type |string           |sinus, step, or exp       |       -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
phi    |array[double]   | order parameter values   |  [num_pts x num seeds]
### Called by
MPF.m

## Depedencies
Al_Sn_Gibbs
DataVis

---
## nrst_nbr_list
### Description:
Calculates nearest neighbors within the active interphase region.

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
coord         |  array[double]     | coordinates      |            [dim X num_pts]
nblist        | cell{array[int]}   | list of nbrs     |          {num_pts[num_nbrs]}
vec_dist_list | cell{array[double]}| vector list      |    {num_pts[dim x num_nbrs]}]

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
nrst_nbr_list |   cell{int}     |    nearest_nbrs    |           {size(nblist)}

### Called by
MPF.m


----
## periodic_coordinates_v2.m
### Description:
 Modifies coordinates to handle periodic BCs

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 Variable Name  |Data Type        |Description        | size
 Length         |double           |domain length      | -
 Height         |double           |domain height      | -
 coord          |array[double]    |coordinates        |[dim X num_pts]
 dim            |int              |spatial dimension  | -
 periodic_bc_on |bool             |for periodic BCs   | -
 num_nbr        |int              |number of neighbors| -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord          |array[double]    |coordinates        |[dim X num_pts]
 exp_coord      |array[double]    |accounts for BCs   |size(coord)

### Called by
get_neighbors_gen.m

### Dependencies
- periodic_coordinates_v2
- KDTreeSearcher  (Requieres Matlab ML and Statistics toolbox)
- knnsearch       (Requieres Matlab ML and Statistics toolbox)
- dialation
----

## random_seed.m
### Description:
Randomly populates seed locations up to num_seed_terms

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 coord          |array[double]    |coordinates              |[dim X num_pts]
 PF_dat         |structure        |phase field input data   |-
 num_seed_terms |double           |number of seeds          |-
 radius         |double           |interphase radius        |-

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 PF_dat        | structure       | phase field data modified  -
 seed_dat      | structure       | seed data     -

## Called by
main.m
----
## taylor_poly_vec_P.m
### Description: Calculates coeffiecents for taylor polynomial terms eg (1/!2*delX)

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 dim            | double         | spatial dimension             |   -
 Lp             | int            | number of terms polynomial    |   -
 poly_order     | int            | order of taylor polynomail    |   -

### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|
 P              | array[double]  | array of polynomial constants   | [Lp x 1]
 coeff_arry     | array[double]  | array of polynomial factorials  | [Lp x 1]
 Pnorm          | array[double]  | distance (x-x0), (y-y0),...     | [Lp x 1]

### Called by
calc_diff_ops.m

---
## File_name
### Description:

### Inputs:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|


### Return Values:
| NAME | DATA TYPE | DESCRIPTION|SIZE|
|------|:----------|:-----------|:-----|

### Called by
main.m
## Depedencies







