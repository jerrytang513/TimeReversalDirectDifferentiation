# Demonstration code: Time reversal differentiation of FDTD for photonic inverse design

Rui Jie Tang (1,†), Soon Wei Daniel Lim (2,†), Marcus Ossiander (2,3), Xinghui Yin (4), Federico Capasso (2)

(1) University of Toronto, 27 King’s College Circle, Toronto, Ontario M5S 1A1, Canada
(2) Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University, Cambridge, MA 02138, USA
(3) Institute of Experimental Physics, Graz University of Technology, 8010 Graz, Austria
(4) LIGO – Massachusetts Institute of Technology, Cambridge, MA 02139, USA
(†) Equal contribution
(*) Email: ruijie.tang@mail.utoronto.ca, lim982@g.harvard.edu 

Manuscript in ACS Photonics: https://doi.org/10.1021/acsphotonics.3c00694

## Description

This is demonstration code for the article "Time reversal differentiation of FDTD for photonic inverse design" by RJ Tang et al. 

The main algorithm for the FDTD gradient calculation is contained in the FDTD3D.cpp class as FDTD3D::reverse_mode_sim_phi_target_half_plane. The function FDTD3D::forward_sim_phi_target_half_plane performs the forward simulation to get the complex fields. Most simulation settings are specified in the header file FDTD3D.h. The main file main.cpp contains a demonstration script for validating the numerical gradients obtained against that of finite difference approximations.

### Simulation description

The simulation geometry is identical to that of the bichromatic color sorter in the article. 

There are two settings for this FDTD controlled by the bool isPeriodic in main.cpp. If isPeriodic is false, the transverse YZ boundary conditions are TFSF/PML boundaries with a recording layer. If isPeriodic is true, the transverse YZ boundary conditions are periodic with no recording layer. The X boundary conditions are always PML. 

* Extent (excluding PMLs): 184 x 50 x 70 pixels (2530 nm x 687.5 nm x 962.5 nm)
* PMLs: 10 layers, all boundaries
* TFSF: 6 pixels away from the PMLs, all boundaries.
* Materials: Glass (Up to x pixel = 42, 577.5 nm), Air (from x pixel = 43 onwards)
* Pillars: A compact 30x60 array, where each pillar has XYZ dimensions of 72 x 1 x 1 pixels (990 nm x 13.75 nm x 13.75 nm) so that the full array has an extent of 72 x 30 x 60 pixels (990 nm x 412.5 nm x 824 nm). The optimization tunable parameters are the permittivities of these 30×60=1800 pillars.  
* Simulation Duration: 700 timesteps.

### main.cpp

The main file provided (main.cpp) performs the following actions:

1. Runs direct differentiation of FDTD and computes the 30x60=1800 element objective function gradient. 
2. Writes the objective function gradient into the output file output.txt.
3. Validates the objective function gradient from the minimal memory calculation by performing explicit finite difference calculations for three positions on the 30x60 pixel array. 

## System requirements

The main code has been tested on the following platforms:

* MacOS Ventura 13.4.1 (2.4 GHz Intel Core i9)
  
  * Peak memory usage: 1.30 GB
  * Time taken: 7.5 minutes
  
* Microsoft Windows 10 Home 10.0.19044 (Intel(R) Xeon(R) CPU E3-1505M v6 @ 3.00GHz, 4 cores)
  * Peak memory usage: 1.325 GB
  * Time taken: 11.2 minutes

## Compilation and execution

### MacOS Terminal

```{bash}
g++ -std=c++11 ./main.cpp ./FDTD3D.cpp -O2 -o ./mmdd
./mmdd
```

### Windows (Visual Studio Command Prompt)

```{bash}
cl /EHsc main.cpp FDTD3D.cpp /O2 /Fe:mmdd.exe
mmdd.exe
```
