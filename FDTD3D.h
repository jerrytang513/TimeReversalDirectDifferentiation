#ifndef MMDD_FDTD3D_H
#define MMDD_FDTD3D_H
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <time.h>

using namespace std;

class FDTD3D {
    bool isPeriodic;
    // default wavelength (m)
    double wavelength = 488e-9;
    double c = 299792458;

    // Simulation Parameters
    double delta_x = 55e-9 / 4.0;
    double delta_t = 0.9 / sqrt(3) * delta_x / c;
    //int time_steps = int(50.e-15/delta_t);
    int time_steps = 700;

    // number of pml layers
    int pml_x = 10;
    int pml_y = 10;
    int pml_z = 10;
    // define the spacing between the pml and the pillar (structure) array
    int pillar_spacing_to_pml = 10;
    // define the pillar (structure) array size
    int pillar_y_size = 30;
    int pillar_z_size = 60;
    // define the x extent of the simulation
    int size_x = 2*pml_x + 184;
    int size_y = (pillar_y_size + 2*pillar_spacing_to_pml) + 2 * pml_y; // 30 dof pillars + 10 pixels between pillars and pml on each side + 2 * pml
    int size_z = (pillar_z_size + 2*pillar_spacing_to_pml) + 2 * pml_z;

    double freq = c / wavelength;
    double Sc = c * delta_t / delta_x;

    // Simulation Region Settings
    // Glass Region
    int glass_start_x = 0; // Inclusive
    int glass_thickness_x = 42;
    int glass_end_x = pml_x + glass_thickness_x; // Exclusive

    // Pillars
    int pillar_height_x = 72;
    int structure_start_x = glass_end_x; // Inclusive
    int structure_end_x = glass_end_x + pillar_height_x; // Exclusive
    int structure_start_y = pml_y + pillar_spacing_to_pml; // Inclusive
    int structure_end_y = size_y - pml_y - pillar_spacing_to_pml; // Exclusive
    int structure_start_z = pml_z + pillar_spacing_to_pml; // Inclusive
    int structure_end_z = size_z - pml_z - pillar_spacing_to_pml; // Exclusive

    // Frequency Plane
    int target_freq_x = pml_x + 164;
    int target_freq_y_start = structure_start_y; // Inclusive
    int target_freq_y_end = structure_end_y; // Exclusive
    int target_freq_z_start = structure_start_z; // Inclusive
    int target_freq_z_end = structure_end_z; // Exclusive

    double epsilon_0 = 8.854187e-12;
    double mu_0 = 1.256637e-6;

    // y-z plane
    vector<vector<double>> epsilon_r_val;
    // convert epsilon_r_val to actual permittivity value
    vector<vector<double>> eps_var;
    // 3D permittivity for entire region
    vector<vector<vector<double>>> eps_arr;
    // Inverse of eps_arr
    vector<vector<vector<double>>> inv_eps_arr;
    // Inverse Square of eps_arr
    vector<vector<vector<double>>> invs_eps_arr;

    double im_avg;
    double re_avg;

    // tfsf boundaries
    int firstX = pml_x + 6; // Inclusive
    int lastX = size_x - pml_x - 6; // Exclusive
    int firstY = pml_y + 6; // Inclusive
    int lastY = size_y - pml_y - 6; // Exclusive
    int firstZ = pml_z + 6; // Inclusive
    int lastZ = size_z - pml_z - 6; // Exclusive

    // TFSF source position
    int tfsf_source_start = firstX - 3;

    double tio2_permittivity = 2.404*2.404;
    double glass_permittivity = 1.44*1.44;

    vector<double> sigma_dx;
    vector<double> sigma_dy;
    vector<double> sigma_dz;
    vector<double> sigma_hx;
    vector<double> sigma_hy;
    vector<double> sigma_hz;

    vector<vector<vector<double>>> hx_arr;
    vector<vector<vector<double>>> hy_arr;
    vector<vector<vector<double>>> hz_arr;
    vector<vector<vector<double>>> dz_arr;
    vector<vector<vector<double>>> ez_arr;
    vector<vector<vector<double>>> dx_arr;
    vector<vector<vector<double>>> ex_arr;
    vector<vector<vector<double>>> dy_arr;
    vector<vector<vector<double>>> ey_arr;

    vector<vector<vector<double>>> mHx1;
    vector<vector<vector<double>>> mHx2;
    vector<vector<vector<double>>> mHx3;
    vector<vector<vector<double>>> mHx4;

    vector<vector<vector<double>>> mHy1;
    vector<vector<vector<double>>> mHy2;
    vector<vector<vector<double>>> mHy3;
    vector<vector<vector<double>>> mHy4;

    vector<vector<vector<double>>> mHz1;
    vector<vector<vector<double>>> mHz2;
    vector<vector<vector<double>>> mHz3;
    vector<vector<vector<double>>> mHz4;

    vector<vector<vector<double>>> mDx1;
    vector<vector<vector<double>>> mDx2;
    vector<vector<vector<double>>> mDx3;
    vector<vector<vector<double>>> mDx4;

    vector<vector<vector<double>>> mDy1;
    vector<vector<vector<double>>> mDy2;
    vector<vector<vector<double>>> mDy3;
    vector<vector<vector<double>>> mDy4;

    vector<vector<vector<double>>> mDz1;
    vector<vector<vector<double>>> mDz2;
    vector<vector<vector<double>>> mDz3;
    vector<vector<vector<double>>> mDz4;

    vector<vector<vector<double>>> ICEx;
    vector<vector<vector<double>>> ICEy;
    vector<vector<vector<double>>> ICEz;

    vector<vector<vector<double>>> ICHx;
    vector<vector<vector<double>>> ICHy;
    vector<vector<vector<double>>> ICHz;

    vector<vector<vector<double>>> IHx;
    vector<vector<vector<double>>> IHy;
    vector<vector<vector<double>>> IHz;

    vector<vector<vector<double>>> IDx;
    vector<vector<vector<double>>> IDy;
    vector<vector<vector<double>>> IDz;

    // Recording
    // top bot z-coordinate
    // left right y-coordinate
    // front back x-coordinate

    vector<vector<vector<double>>> hx_top;
    vector<vector<vector<double>>> hx_bot;
    vector<vector<vector<double>>> hx_left;
    vector<vector<vector<double>>> hx_right;
    vector<vector<vector<double>>> hx_front;
    vector<vector<vector<double>>> hx_back;

    vector<vector<vector<double>>> hy_top;
    vector<vector<vector<double>>> hy_bot;
    vector<vector<vector<double>>> hy_left;
    vector<vector<vector<double>>> hy_right;
    vector<vector<vector<double>>> hy_front;
    vector<vector<vector<double>>> hy_back;

    vector<vector<vector<double>>> hz_top;
    vector<vector<vector<double>>> hz_bot;
    vector<vector<vector<double>>> hz_left;
    vector<vector<vector<double>>> hz_right;
    vector<vector<vector<double>>> hz_front;
    vector<vector<vector<double>>> hz_back;

    vector<vector<vector<double>>> dx_top;
    vector<vector<vector<double>>> dx_bot;
    vector<vector<vector<double>>> dx_left;
    vector<vector<vector<double>>> dx_right;
    vector<vector<vector<double>>> dx_front;
    vector<vector<vector<double>>> dx_back;

    vector<vector<vector<double>>> dy_top;
    vector<vector<vector<double>>> dy_bot;
    vector<vector<vector<double>>> dy_left;
    vector<vector<vector<double>>> dy_right;
    vector<vector<vector<double>>> dy_front;
    vector<vector<vector<double>>> dy_back;

    vector<vector<vector<double>>> dz_top;
    vector<vector<vector<double>>> dz_bot;
    vector<vector<vector<double>>> dz_left;
    vector<vector<vector<double>>> dz_right;
    vector<vector<vector<double>>> dz_front;
    vector<vector<vector<double>>> dz_back;

    // Derivative Calculation
    vector<vector<vector<vector<double>>>> storeDx;
    vector<vector<vector<vector<double>>>> storeDy;
    vector<vector<vector<vector<double>>>> storeDz;

    vector<vector<vector<double>>> dL_dHx;
    vector<vector<vector<double>>> dL_dHy;
    vector<vector<vector<double>>> dL_dHz;
    vector<vector<vector<double>>> dL_dDx;
    vector<vector<vector<double>>> dL_dDy;
    vector<vector<vector<double>>> dL_dDz;
    vector<vector<vector<double>>> dL_dEx;
    vector<vector<vector<double>>> dL_dEy;
    vector<vector<vector<double>>> dL_dEz;

    vector<vector<vector<double>>> dL_dIHx;
    vector<vector<vector<double>>> dL_dIHy;
    vector<vector<vector<double>>> dL_dIHz;
    vector<vector<vector<double>>> dL_dIDx;
    vector<vector<vector<double>>> dL_dIDy;
    vector<vector<vector<double>>> dL_dIDz;

    vector<vector<vector<double>>> dL_dICEx;
    vector<vector<vector<double>>> dL_dICEy;
    vector<vector<vector<double>>> dL_dICEz;
    vector<vector<vector<double>>> dL_dICHx;
    vector<vector<vector<double>>> dL_dICHy;
    vector<vector<vector<double>>> dL_dICHz;

    vector<vector<double>> dL_deps;
    vector<vector<double>> dL_depsr;
    vector<vector<double>> real_fft_ex;
    vector<vector<double>> imag_fft_ex;

    vector<vector<double>> real_fft_ey;
    vector<vector<double>> imag_fft_ey;

    vector<vector<double>> real_fft_ez;
    vector<vector<double>> imag_fft_ez;


    // 1D TFSF
    vector<double> hy1D_arr;
    vector<double> ez1D_arr;
    vector<double> hy1D_coef1;
    vector<double> hy1D_coef2;
    vector<double> ez1D_coef1;
    vector<double> ez1D_coef2;
    double max_loss = 0.35;


public:
    FDTD3D(bool isPeriodic);
    // phi forward simulation, returns the phase phi
    double forward_sim_phi_target_half_plane(int plane_num);
    void backward_sim_step(int cur_timestep);
    double reverse_mode_sim_phi_target_half_plane(int plane_num);
    double get_reverse_mode_sim_phi_derivative(int y, int z);
    double get_epsilon_r_val(int y, int z);
    void init_epsilon_r_val();
    void vacuum_epsilon_r_val();
    void set_one_permittivity(int y, int z, double val);
private:
    // Gaussian Function
    double gaussian(int x, int y, int center_x, int center_y, int start_x, int start_y, int radius);
    // Generate matrices for field
    void init_field();
    // Reset all field values
    void set_field();
    // Generate matrices for derivative
    void init_derivative_field();
    // Reset all field values
    void set_derivative_field();
    void tfsf(int t);
    void init_recording_field();
    // Generate matrices for pml and field
    void init_mat();
    // Generate matrices for pml
    void init_pml();
    // Set up permittivities
    void set_permittivity();
    // Set individual permittivity
    void record_field(int cur_timestep);
    double sineFn(double time);

};


#endif //MMDD_FDTD3D_H
