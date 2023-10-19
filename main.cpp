#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include "FDTD3D.h"
#include <cstring>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <fstream>

// Default simulation parameters from FDTD3D.h

// Define the limits of the simulation region
int pml = 10;
// define the spacing between the pml and the pillar array
int pillar_spacing_to_pml = 10;
// define the pillar (structure) array size
int pillar_y_size = 30;
int pillar_z_size = 60;
int size_y = (pillar_y_size + 2*pillar_spacing_to_pml) + 2 * pml; // 30 dof pillars + 10 pixels between pillars and pml on each side + 2 * pml
int size_z = (pillar_z_size + 2*pillar_spacing_to_pml) + 2 * pml;
int structure_start_y = pml + pillar_spacing_to_pml; // Inclusive
int structure_end_y = size_y - pml - pillar_spacing_to_pml; // Exclusive
int structure_start_z = pml + pillar_spacing_to_pml; // Inclusive
int structure_end_z = size_z - pml - pillar_spacing_to_pml; // Exclusive
int offset_y = pml + pillar_spacing_to_pml;
int offset_z = pml + pillar_spacing_to_pml;

double finite_diff_result(int y, int z, FDTD3D& network3){
    // performs a double-sided first order finite difference approximation to the gradient
    // associated with perturbing the dielectric permittivity for one pixel located at pixel (y,z).
    double perturbation = 0.000001; // size of finite difference perturbation
    network3.set_one_permittivity(y, z, network3.get_epsilon_r_val(y, z)-perturbation/2);
    double objfun1 = network3.forward_sim_phi_target_half_plane(1);
    network3.set_one_permittivity(y, z, network3.get_epsilon_r_val(y, z)+perturbation);
    double objfun2 = network3.forward_sim_phi_target_half_plane(1);
    double finite_diff = (objfun2-objfun1)/perturbation;
    // reset the dielectric permittivity to the original value
    network3.set_one_permittivity(y, z, network3.get_epsilon_r_val(y, z)-perturbation/2);
    return finite_diff;
}

void load_random_permittivity(FDTD3D& network3){
    // Set everything to vacuum first
    network3.vacuum_epsilon_r_val();
    srand (time(NULL));

    // Set a random permittivity value for the 30 x 60 array
    for (int j = 0; j < size_y; j++) {
        for (int k = 0; k < size_z; k++) {
            // Region outside of the structure will be set to vacuum
            // (vacuum corresponds to -10 in latent permittivity)
            network3.set_one_permittivity(j, k, -10.);
            // Randomize permittivity for these structures
            if (j >= structure_start_y && j < structure_end_y &&
                k >= structure_start_z && k < structure_end_z) {
                network3.set_one_permittivity(j, k, (((double) rand() / (RAND_MAX)) - 0.5) * 20);
            }
        }
    }
}

void do_single_pixel_finite_diff_comparison(int y, int z, FDTD3D& network3){

    double finite_diff_grad = finite_diff_result(y+offset_y, z+offset_z, network3);
    double mmdd_grad = network3.get_reverse_mode_sim_phi_derivative(y+offset_y, z+offset_z);

    cout << "Finite Difference gradient at pixel x=174,y=" << y+offset_y << ",z=" << z+offset_z << ": " << setprecision(20) << finite_diff_grad << endl;
    cout << "Minimal memory AD gradient at pixel x=174,y=" << y+offset_y << ",z=" << z+offset_z << ": " << setprecision(20) << mmdd_grad << endl;
    cout << "Diff: " << finite_diff_grad - mmdd_grad << " (" << 100*(finite_diff_grad - mmdd_grad)/mmdd_grad << "%)" << endl;
}

int main(int argc, char *argv[]) {
    // Perform gradient accuracy benchmark between Finite Difference approximation and
    // Minimal Memory Reverse Mode Direct Differentiation
    bool isPeriodic = false;
    FDTD3D network3 = FDTD3D(isPeriodic);

    // load a random permittivity distribution over the 30x60 degrees of freedom pixels
    load_random_permittivity(network3);

    // Run the reverse mode direct differentiation
    network3.reverse_mode_sim_phi_target_half_plane(1);

    // Output the 1800 derivative results into a text file
    std::ofstream outfile("output.txt");
    if (outfile.is_open()) {
        for (int y = 0; y < 30; y++) {
            for (int z = 0; z < 60; z++) {
                outfile << setprecision(20) << network3.get_reverse_mode_sim_phi_derivative(y + offset_y, z + offset_z)
                        << endl;
            }
        }
    }
    outfile.close();
    cout << "Minimal memory differentiation of FDTD complete; gradient saved in output.txt file. " << endl;
    cout << "Comparing gradient obtained to finite difference approximation for 3 pixels..." << endl;

    // Pick 3 coordinates to compare the difference between finite difference result and reverse mode gradient result
    int y1 = 20;
    int z1 = 20;
    int y2 = 25;
    int z2 = 25;
    int y3 = 27;
    int z3 = 25;

    cout << "--- Point 1 finite difference validation ---" << endl;
    do_single_pixel_finite_diff_comparison(y1, z1, network3);
    cout << "--- Point 2 finite difference validation ---" << endl;
    do_single_pixel_finite_diff_comparison(y2, z2, network3);
    cout << "--- Point 3 finite difference validation ---" << endl;
    do_single_pixel_finite_diff_comparison(y3, z3, network3);
}
