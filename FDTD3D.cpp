#include "FDTD3D.h"
#include <fstream>
#include <iostream>

FDTD3D::FDTD3D(bool isPeriodic){
    this->isPeriodic = isPeriodic;
    if(isPeriodic){
        pml_y = 0;
        pml_z = 0;
    }
    init_mat();
    init_epsilon_r_val();
    init_pml();
}

double FDTD3D::forward_sim_phi_target_half_plane(int plane_num){
    // UPML formulation and initialization process refers to Ceviche and Empossible
    // Ceviche: https://github.com/fancompute/ceviche  paper: https://doi.org/10.1021/acsphotonics.9b01238
    // Empossible: https://www.youtube.com/watch?v=w_NnRZlNuAA&t=2454s
    // TFSF formulation refers to Understanding the FDTD Method by John B. Schneider https://eecs.wsu.edu/~schneidj/ufdtd/

    set_field();
    set_permittivity();
    double energy_sum = 0.;

    for(int cur_timestep = 0; cur_timestep < time_steps; cur_timestep++) {

        // update hx
        for (int x = 1; x < size_x - 1; x++) {
            for (int y = 1; y < size_y - 1; y++) {
                for (int z = 1; z < size_z - 1; z++) {
                    int xl = x;
                    int yl = y;
                    int zl = z;
                    if (y == size_y - 2 && isPeriodic) {
                        yl = 0;
                    }
                    if (z == size_z - 2 && isPeriodic) {
                        zl = 0;
                    }
                    double curl_e_x =
                            (ey_arr[x][y][zl + 1] - ey_arr[x][y][z] + ez_arr[x][y][z] - ez_arr[x][yl + 1][z]) / delta_x;
                    double curl_e_y =
                            (ez_arr[xl + 1][y][z] - ez_arr[x][y][z] + ex_arr[x][y][z] - ex_arr[x][y][zl + 1]) / delta_x;
                    double curl_e_z =
                            (ex_arr[x][yl + 1][z] - ex_arr[x][y][z] + ey_arr[x][y][z] - ey_arr[xl + 1][y][z]) / delta_x;
                    ICEx[x][y][z] += curl_e_x;
                    ICEy[x][y][z] += curl_e_y;
                    ICEz[x][y][z] += curl_e_z;
                    IHx[x][y][z] += hx_arr[x][y][z];
                    IHy[x][y][z] += hy_arr[x][y][z];
                    IHz[x][y][z] += hz_arr[x][y][z];
                    hx_arr[x][y][z] =
                            mHx1[x][y][z] * hx_arr[x][y][z] + mHx2[x][y][z] * curl_e_x + mHx3[x][y][z] * ICEx[x][y][z] +
                            mHx4[x][y][z] * IHx[x][y][z];
                    hy_arr[x][y][z] =
                            mHy1[x][y][z] * hy_arr[x][y][z] + mHy2[x][y][z] * curl_e_y + mHy3[x][y][z] * ICEy[x][y][z] +
                            mHy4[x][y][z] * IHy[x][y][z];
                    hz_arr[x][y][z] =
                            mHz1[x][y][z] * hz_arr[x][y][z] + mHz2[x][y][z] * curl_e_z + mHz3[x][y][z] * ICEz[x][y][z] +
                            mHz4[x][y][z] * IHz[x][y][z];
                }
            }
        }


        if(!isPeriodic){
            tfsf(cur_timestep);
        } else {
            // plane wave simulation
            for (int y = 1; y <= size_y-2; y++) {
                for (int z = 1; z <= size_z-2; z++) {
                    hy_arr[tfsf_source_start][y][z] += sineFn(cur_timestep);
                }
            }
        }

        // Update D field
        for (int x = 1; x < size_x - 1; x++) {
            for (int y = 1; y < size_y - 1; y++) {
                for (int z = 1; z < size_z - 1; z++) {

                    int xl = x;
                    int yl = y;
                    int zl = z;
                    if (y == 1 && isPeriodic) {
                        yl = size_y-1;
                    }
                    if (z == 1 && isPeriodic) {
                        zl = size_z-1;
                    }

                    double curl_h_x =
                            (hz_arr[x][y][z] - hz_arr[x][yl - 1][z] - hy_arr[x][y][z] + hy_arr[x][y][zl - 1]) / delta_x;
                    double curl_h_y =
                            (hx_arr[x][y][z] - hx_arr[x][y][zl - 1] - hz_arr[x][y][z] + hz_arr[xl - 1][y][z]) / delta_x;
                    double curl_h_z =
                            (hy_arr[x][y][z] - hy_arr[xl - 1][y][z] - hx_arr[x][y][z] + hx_arr[x][yl - 1][z]) / delta_x;
                    ICHx[x][y][z] += curl_h_x;
                    ICHy[x][y][z] += curl_h_y;
                    ICHz[x][y][z] += curl_h_z;
                    IDx[x][y][z] += dx_arr[x][y][z];
                    IDy[x][y][z] += dy_arr[x][y][z];
                    IDz[x][y][z] += dz_arr[x][y][z];

                    dx_arr[x][y][z] =
                            mDx1[x][y][z] * dx_arr[x][y][z] + mDx2[x][y][z] * curl_h_x + mDx3[x][y][z] * ICHx[x][y][z] +
                            mDx4[x][y][z] * IDx[x][y][z];
                    dy_arr[x][y][z] =
                            mDy1[x][y][z] * dy_arr[x][y][z] + mDy2[x][y][z] * curl_h_y + mDy3[x][y][z] * ICHy[x][y][z] +
                            mDy4[x][y][z] * IDy[x][y][z];
                    dz_arr[x][y][z] =
                            mDz1[x][y][z] * dz_arr[x][y][z] + mDz2[x][y][z] * curl_h_z + mDz3[x][y][z] * ICHz[x][y][z] +
                            mDz4[x][y][z] * IDz[x][y][z];
                }
            }
        }

        // Use the recording boundary condition
        record_field(cur_timestep);

        // Update E field
        for(int x = 1; x < size_x-1; x ++) {
            for (int y = 1; y < size_y-1; y++) {
                for (int z = 1; z < size_z-1; z++) {
                    ex_arr[x][y][z] = dx_arr[x][y][z] * inv_eps_arr[x][y][z];
                    ey_arr[x][y][z] = dy_arr[x][y][z] * inv_eps_arr[x][y][z];
                    ez_arr[x][y][z] = dz_arr[x][y][z] * inv_eps_arr[x][y][z];
                }
            }
        }
        if(!isPeriodic){
            for(int x=0; x <size_x; x++){
                for(int y=0; y <size_y; y++) {
                    for (int z = 0; z < size_z; z++) {
                        if (x == 0 || x == size_x - 1 || y == 0 || y == size_y - 1 || z == 0 || z == size_z - 1) {
                            ex_arr[x][y][z] = 0;
                            ey_arr[x][y][z] = 0;
                            ez_arr[x][y][z] = 0;
                            dx_arr[x][y][z] = 0;
                            dy_arr[x][y][z] = 0;
                            dz_arr[x][y][z] = 0;
                            hx_arr[x][y][z] = 0;
                            hy_arr[x][y][z] = 0;
                            hz_arr[x][y][z] = 0;
                        }
                    }
                }
            }
        }

        energy_sum = 0;
        for(int x=0; x <size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                for (int z = 0; z < size_z; z++) {
                    energy_sum += ez_arr[x][y][z] * ez_arr[x][y][z];
                }
            }
        }

        // Evaluate FT
        for(int y = target_freq_y_start; y < target_freq_y_end; y++){
            for(int z = target_freq_z_start; z < target_freq_z_end; z++) {
                real_fft_ez[y-target_freq_y_start][z-target_freq_z_start] += cos(2 * M_PI * freq * delta_t * cur_timestep) * ez_arr[target_freq_x][y][z];
                imag_fft_ez[y-target_freq_y_start][z-target_freq_z_start] += sin(2 * M_PI * freq * delta_t * cur_timestep) * ez_arr[target_freq_x][y][z];
            }
        }

    }
    double sum = 0.;
    if(plane_num == 0) {
        // Then the objective function is the total weighted intensity in the left half plane
        for(int y = 0; y < 30; y++){
            for(int z = 0; z < 30; z++) {
                sum += (real_fft_ez[y][z] * real_fft_ez[y][z] + imag_fft_ez[y][z] * imag_fft_ez[y][z]) * gaussian(y, z, 15, 10, 0, 0, 8);
            }
        }
    } else {
        // Then the objective function is the total weighted intensity in the right half plane
        for(int y = 0; y < 30; y++){
            for(int z = 30; z < 60; z++) {
                sum += (real_fft_ez[y][z] * real_fft_ez[y][z] + imag_fft_ez[y][z] * imag_fft_ez[y][z]) * gaussian(y, z, 15, 10, 0, 40, 8);
            }
        }
    }
    return sum;
}

void FDTD3D::backward_sim_step(int cur_timestep) {

    // Inject recording boundary components back into the simulator
    // init d boundary (top, bot, left, right, front, back)
    // top z
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y ++){
            hx_arr[x][y][structure_end_z-1] = hx_top[x-structure_start_x][y-structure_start_y][cur_timestep];
            hy_arr[x][y][structure_end_z-1] = hy_top[x-structure_start_x][y-structure_start_y][cur_timestep];
            hz_arr[x][y][structure_end_z-1] = hz_top[x-structure_start_x][y-structure_start_y][cur_timestep];
        }
    }
    // bottom z
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y ++){
            hx_arr[x][y][structure_start_z] = hx_bot[x-structure_start_x][y-structure_start_y][cur_timestep];
            hy_arr[x][y][structure_start_z] = hy_bot[x-structure_start_x][y-structure_start_y][cur_timestep];
            hz_arr[x][y][structure_start_z] = hz_bot[x-structure_start_x][y-structure_start_y][cur_timestep];
        }
    }
    // left y
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_arr[x][structure_start_y][z] = hx_left[x-structure_start_x][z-structure_start_z][cur_timestep];
            hy_arr[x][structure_start_y][z] = hy_left[x-structure_start_x][z-structure_start_z][cur_timestep];
            hz_arr[x][structure_start_y][z] = hz_left[x-structure_start_x][z-structure_start_z][cur_timestep];
        }
    }
    // right y
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_arr[x][structure_end_y-1][z] = hx_right[x-structure_start_x][z-structure_start_z][cur_timestep];
            hy_arr[x][structure_end_y-1][z] = hy_right[x-structure_start_x][z-structure_start_z][cur_timestep];
            hz_arr[x][structure_end_y-1][z] = hz_right[x-structure_start_x][z-structure_start_z][cur_timestep];
        }
    }
    // front x
    for(int y = structure_start_y; y < structure_end_y; y ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_arr[structure_start_x][y][z] = hx_front[y-structure_start_y][z-structure_start_z][cur_timestep];
            hy_arr[structure_start_x][y][z] = hy_front[y-structure_start_y][z-structure_start_z][cur_timestep];
            hz_arr[structure_start_x][y][z] = hz_front[y-structure_start_y][z-structure_start_z][cur_timestep];
        }
    }
    // back x
    for(int y = structure_start_y; y < structure_end_y; y ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_arr[structure_end_x-1][y][z] = hx_back[y-structure_start_y][z-structure_start_z][cur_timestep];
            hy_arr[structure_end_x-1][y][z] = hy_back[y-structure_start_y][z-structure_start_z][cur_timestep];
            hz_arr[structure_end_x-1][y][z] = hz_back[y-structure_start_y][z-structure_start_z][cur_timestep];
        }
    }
    // Update D field
    for(int x = structure_start_x+1; x < structure_end_x - 1; x ++){
        for(int y = structure_start_y+1; y < structure_end_y-1; y++){
            for(int z = structure_start_z+1; z < structure_end_z-1; z++){
                double curl_h_x = (hz_arr[x][y][z] - hz_arr[x][y-1][z] - hy_arr[x][y][z] + hy_arr[x][y][z-1]) / delta_x;
                double curl_h_y = (hx_arr[x][y][z] - hx_arr[x][y][z-1] - hz_arr[x][y][z] + hz_arr[x-1][y][z]) / delta_x;
                double curl_h_z = (hy_arr[x][y][z] - hy_arr[x-1][y][z] - hx_arr[x][y][z] + hx_arr[x][y-1][z]) / delta_x;

                dx_arr[x][y][z] -= mDx2[x][y][z] * curl_h_x;
                dy_arr[x][y][z] -= mDy2[x][y][z] * curl_h_y;
                dz_arr[x][y][z] -= mDz2[x][y][z] * curl_h_z;
            }
        }
    }

    // top z
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y ++){
            dx_arr[x][y][structure_end_z-1] = dx_top[x-structure_start_x][y-structure_start_y][cur_timestep-1];
            dy_arr[x][y][structure_end_z-1] = dy_top[x-structure_start_x][y-structure_start_y][cur_timestep-1];
            dz_arr[x][y][structure_end_z-1] = dz_top[x-structure_start_x][y-structure_start_y][cur_timestep-1];
        }
    }
    // bot z
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y ++){
            dx_arr[x][y][structure_start_z] = dx_bot[x-structure_start_x][y-structure_start_y][cur_timestep-1];
            dy_arr[x][y][structure_start_z] = dy_bot[x-structure_start_x][y-structure_start_y][cur_timestep-1];
            dz_arr[x][y][structure_start_z] = dz_bot[x-structure_start_x][y-structure_start_y][cur_timestep-1];
        }
    }
    // left y
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            dx_arr[x][structure_start_y][z] = dx_left[x-structure_start_x][z-structure_start_z][cur_timestep-1];
            dy_arr[x][structure_start_y][z] = dy_left[x-structure_start_x][z-structure_start_z][cur_timestep-1];
            dz_arr[x][structure_start_y][z] = dz_left[x-structure_start_x][z-structure_start_z][cur_timestep-1];
        }
    }
    // right y
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            dx_arr[x][structure_end_y-1][z] = dx_right[x-structure_start_x][z-structure_start_z][cur_timestep-1];
            dy_arr[x][structure_end_y-1][z] = dy_right[x-structure_start_x][z-structure_start_z][cur_timestep-1];
            dz_arr[x][structure_end_y-1][z] = dz_right[x-structure_start_x][z-structure_start_z][cur_timestep-1];
        }
    }
    // front x
    for(int y = structure_start_y; y < structure_end_y; y ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            dx_arr[structure_start_x][y][z] = dx_front[y-structure_start_y][z-structure_start_z][cur_timestep-1];
            dy_arr[structure_start_x][y][z] = dy_front[y-structure_start_y][z-structure_start_z][cur_timestep-1];
            dz_arr[structure_start_x][y][z] = dz_front[y-structure_start_y][z-structure_start_z][cur_timestep-1];
        }
    }
    // back x
    for(int y = structure_start_y; y < structure_end_y; y ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            dx_arr[structure_end_x-1][y][z] = dx_back[y-structure_start_y][z-structure_start_z][cur_timestep-1];
            dy_arr[structure_end_x-1][y][z] = dy_back[y-structure_start_y][z-structure_start_z][cur_timestep-1];
            dz_arr[structure_end_x-1][y][z] = dz_back[y-structure_start_y][z-structure_start_z][cur_timestep-1];
        }
    }

    // Update E field
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y++){
            for(int z = structure_start_z; z < structure_end_z; z++){
                ex_arr[x][y][z] = dx_arr[x][y][z] * inv_eps_arr[x][y][z];
                ey_arr[x][y][z] = dy_arr[x][y][z] * inv_eps_arr[x][y][z];
                ez_arr[x][y][z] = dz_arr[x][y][z] * inv_eps_arr[x][y][z];
            }
        }
    }

    // update hx
    for(int x = structure_start_x+1; x < structure_end_x - 1; x ++){
        for(int y = structure_start_y+1; y < structure_end_y-1; y++){
            for(int z = structure_start_z+1; z < structure_end_z-1; z++){
                double curl_e_x = (ey_arr[x][y][z+1] - ey_arr[x][y][z] + ez_arr[x][y][z] - ez_arr[x][y+1][z]) / delta_x;
                double curl_e_y = (ez_arr[x+1][y][z] - ez_arr[x][y][z] + ex_arr[x][y][z] - ex_arr[x][y][z+1]) / delta_x;
                double curl_e_z = (ex_arr[x][y+1][z] - ex_arr[x][y][z] + ey_arr[x][y][z] - ey_arr[x+1][y][z]) / delta_x;
                hx_arr[x][y][z] -= mHx2[x][y][z] * curl_e_x;
                hy_arr[x][y][z] -= mHy2[x][y][z] * curl_e_y;
                hz_arr[x][y][z] -= mHz2[x][y][z] * curl_e_z;
            }
        }
    }
}

double FDTD3D::reverse_mode_sim_phi_target_half_plane(int plane_num){

    // Run the forward simulation to get the target plane fields
    double intensity = forward_sim_phi_target_half_plane(plane_num);
    // Initialize derivative fields
    set_derivative_field();

    double dL_dIm = 2 * re_avg;
    double dL_dRe = 2 * im_avg;

    // Calculate derivatives backwards in time using the chain rule
    for(int cur_time_step = 0; cur_time_step < time_steps; cur_time_step++){
        int t = time_steps - 1 - cur_time_step;

        // d(objective function) / d E(t)
        for(int x = 1; x < size_x-1; x ++){
            for(int y = 1; y < size_y-1; y++){
                for(int z = 1; z < size_z-1; z++){
                    int xl = x;
                    int yl = y;
                    int zl = z;

                    if (y == 1 && isPeriodic) {
                        yl = size_y-1;
                    }
                    if (z == 1 && isPeriodic) {
                        zl = size_z-1;
                    }

                    double reverse_curlEx = - (dL_dHz[x][y][z] * mHz2[x][y][z] - dL_dHz[x][yl-1][z] * mHz2[x][yl-1][z])
                                            + (dL_dHy[x][y][z] * mHy2[x][y][z] - dL_dHy[x][y][zl-1] * mHy2[x][y][zl-1]);
                    double reverse_curlEy = - (dL_dHx[x][y][z] * mHx2[x][y][z] - dL_dHx[x][y][zl-1] * mHx2[x][y][zl-1])
                                            + (dL_dHz[x][y][z] * mHz2[x][y][z] - dL_dHz[xl-1][y][z] * mHz2[xl-1][y][z]);
                    double reverse_curlEz = - (dL_dHy[x][y][z] * mHy2[x][y][z] - dL_dHy[xl-1][y][z] * mHy2[xl-1][y][z])
                                            + (dL_dHx[x][y][z] * mHx2[x][y][z] - dL_dHx[x][yl-1][z] * mHx2[x][yl-1][z]);
                    double reverse_i_curlEx = - (dL_dICEz[x][y][z] - dL_dICEz[x][yl-1][z])
                                              + (dL_dICEy[x][y][z] - dL_dICEy[x][y][zl-1]);
                    double reverse_i_curlEy = - (dL_dICEx[x][y][z] - dL_dICEx[x][y][zl-1])
                                              + (dL_dICEz[x][y][z] - dL_dICEz[xl-1][y][z]);
                    double reverse_i_curlEz = - (dL_dICEy[x][y][z] - dL_dICEy[xl-1][y][z])
                                              + (dL_dICEx[x][y][z] - dL_dICEx[x][yl-1][z]);
                    dL_dEx[x][y][z] = (reverse_curlEx + reverse_i_curlEx) / delta_x;
                    dL_dEy[x][y][z] = (reverse_curlEy + reverse_i_curlEy) / delta_x;
                    dL_dEz[x][y][z] = (reverse_curlEz + reverse_i_curlEz) / delta_x;
                }
            }
        }

        // Derivatives from Fourier Series
        double sum = 0.;

        double real_part = 0.;
        double imag_part = 0.;
        double dphi_dEz = 0.;

        if(plane_num == 0) {
            // Then the objective function is the total weighted intensity in the left half plane
            for(int y = 0; y < 30; y++){
                for(int z = 0; z < 30; z++) {
                    real_part = 2 * real_fft_ez[y][z];
                    imag_part = 2 * imag_fft_ez[y][z];
                    dphi_dEz = real_part * cos(2 * M_PI * freq * delta_t * t) + imag_part * sin(2 * M_PI * freq * delta_t * t);
                    dL_dEz[target_freq_x][y + target_freq_y_start][z + target_freq_z_start] += dphi_dEz * gaussian(y, z, 15, 10, 0, 0, 8);
                }
            }
        } else {
            // Then the objective function is the total weighted intensity in the right half plane
            for(int y = 0; y < 30; y++){
                for(int z = 30; z < 60; z++) {
                    real_part = 2 * real_fft_ez[y][z];
                    imag_part = 2 * imag_fft_ez[y][z];
                    dphi_dEz = real_part * cos(2 * M_PI * freq * delta_t * t) + imag_part * sin(2 * M_PI * freq * delta_t * t);
                    dL_dEz[target_freq_x][y + target_freq_y_start][z + target_freq_z_start] += dphi_dEz * gaussian(y, z, 15, 10, 0, 40, 8);
                }
            }
        }

        for(int i = structure_start_x; i < structure_end_x; i ++){
            for(int j = structure_start_y; j < structure_end_y; j ++){
                for(int k = structure_start_z; k < structure_end_z; k ++){
                    // Use time reversal simulated field for derivative calculation
                    dL_deps[j][k] -= dL_dEx[i][j][k] * invs_eps_arr[i][j][k] * dx_arr[i][j][k];
                    dL_deps[j][k] -= dL_dEy[i][j][k] * invs_eps_arr[i][j][k] * dy_arr[i][j][k];
                    dL_deps[j][k] -= dL_dEz[i][j][k] * invs_eps_arr[i][j][k] * dz_arr[i][j][k];

                }
            }
        }

        // Perform the time-reversal simulation
        backward_sim_step(t);

        // Calculate chain rules from PML
        for(int x = 1; x < size_x-1; x ++){
            for(int y = 1; y < size_y-1; y++) {
                for (int z = 1; z < size_z - 1; z++) {
                    dL_dDx[x][y][z] = dL_dEx[x][y][z] * inv_eps_arr[x][y][z] + dL_dDx[x][y][z] * mDx1[x][y][z] + dL_dIDx[x][y][z];
                    dL_dDy[x][y][z] = dL_dEy[x][y][z] * inv_eps_arr[x][y][z] + dL_dDy[x][y][z] * mDy1[x][y][z] + dL_dIDy[x][y][z];
                    dL_dDz[x][y][z] = dL_dEz[x][y][z] * inv_eps_arr[x][y][z] + dL_dDz[x][y][z] * mDz1[x][y][z] + dL_dIDz[x][y][z];

                    dL_dIDx[x][y][z] = dL_dIDx[x][y][z] + dL_dDx[x][y][z] * mDx4[x][y][z];
                    dL_dIDy[x][y][z] = dL_dIDy[x][y][z] + dL_dDy[x][y][z] * mDy4[x][y][z];
                    dL_dIDz[x][y][z] = dL_dIDz[x][y][z] + dL_dDz[x][y][z] * mDz4[x][y][z];

                    dL_dICHx[x][y][z] = dL_dICHx[x][y][z] + dL_dDx[x][y][z] * mDx3[x][y][z];
                    dL_dICHy[x][y][z] = dL_dICHy[x][y][z] + dL_dDy[x][y][z] * mDy3[x][y][z];
                    dL_dICHz[x][y][z] = dL_dICHz[x][y][z] + dL_dDz[x][y][z] * mDz3[x][y][z];
                }
            }
        }

        // d(objective function) / d H(t)
        for(int x = 1; x < size_x-1; x ++){
            for(int y = 1; y < size_y-1; y++){
                for(int z = 1; z < size_z-1; z++){
                    int xl = x;
                    int yl = y;
                    int zl = z;

                    if (y == size_y - 2 && isPeriodic) {
                        yl = 0;
                    }
                    if (z == size_z - 2 && isPeriodic) {
                        zl = 0;
                    }
                    double reverse_curlHx = (dL_dDz[x][yl+1][z] * mDz2[x][yl+1][z] - dL_dDz[x][y][z] * mDz2[x][y][z])
                                            - (dL_dDy[x][y][zl+1] * mDy2[x][y][zl+1] - dL_dDy[x][y][z] * mDy2[x][y][z]);
                    double reverse_curlHy = (dL_dDx[x][y][zl+1] * mDx2[x][y][zl+1] - dL_dDx[x][y][z] * mDx2[x][y][z])
                                            - (dL_dDz[xl+1][y][z] * mDz2[xl+1][y][z] - dL_dDz[x][y][z] * mDz2[x][y][z]);
                    double reverse_curlHz = (dL_dDy[xl+1][y][z] * mDy2[xl+1][y][z] - dL_dDy[x][y][z] * mDy2[x][y][z])
                                            - (dL_dDx[x][yl+1][z] * mDx2[x][yl+1][z] - dL_dDx[x][y][z] * mDx2[x][y][z]);

                    double reverse_i_curlHx = (dL_dICHz[x][yl+1][z] - dL_dICHz[x][y][z])
                                              - (dL_dICHy[x][y][zl+1] - dL_dICHy[x][y][z]);
                    double reverse_i_curlHy = (dL_dICHx[x][y][zl+1] - dL_dICHx[x][y][z])
                                              - (dL_dICHz[xl+1][y][z] - dL_dICHz[x][y][z]);
                    double reverse_i_curlHz = (dL_dICHy[xl+1][y][z] - dL_dICHy[x][y][z])
                                              - (dL_dICHx[x][yl+1][z] - dL_dICHx[x][y][z]);
                    dL_dHx[x][y][z] = mHx1[x][y][z] * dL_dHx[x][y][z] + (reverse_curlHx + reverse_i_curlHx) / delta_x + dL_dIHx[x][y][z];
                    dL_dHy[x][y][z] = mHy1[x][y][z] * dL_dHy[x][y][z] + (reverse_curlHy + reverse_i_curlHy) / delta_x + dL_dIHy[x][y][z];
                    dL_dHz[x][y][z] = mHz1[x][y][z] * dL_dHz[x][y][z] + (reverse_curlHz + reverse_i_curlHz) / delta_x + dL_dIHz[x][y][z];
                }
            }
        }

        // Calculate chain rules from PML
        for(int x = 1; x < size_x-1; x ++) {
            for (int y = 1; y < size_y - 1; y++) {
                for (int z = 1; z < size_z - 1; z++) {
                    dL_dIHx[x][y][z] = dL_dIHx[x][y][z] + dL_dHx[x][y][z] * mHx4[x][y][z];
                    dL_dIHy[x][y][z] = dL_dIHy[x][y][z] + dL_dHy[x][y][z] * mHy4[x][y][z];
                    dL_dIHz[x][y][z] = dL_dIHz[x][y][z] + dL_dHz[x][y][z] * mHz4[x][y][z];

                    dL_dICEx[x][y][z] = dL_dICEx[x][y][z] + dL_dHx[x][y][z] * mHx3[x][y][z];
                    dL_dICEy[x][y][z] = dL_dICEy[x][y][z] + dL_dHy[x][y][z] * mHy3[x][y][z];
                    dL_dICEz[x][y][z] = dL_dICEz[x][y][z] + dL_dHz[x][y][z] * mHz3[x][y][z];
                }
            }
        }
    }

    // Convert dL_deps to dL_depsr.  epsr is scaled and bounded to the actual permittivity (eps) through tanh functions.
    for(int y = 0; y < size_y; y++){
        for(int z = 0; z < size_z; z++){
            dL_depsr[y][z] = dL_deps[y][z] * (tio2_permittivity - 1.0) / 4.0 * (1.0 - pow(tanh(0.5 * epsilon_r_val[y][z]),2));
        }
    }
    // Return the phase. User can get the gradient through the function get_reverse_mode_sim_phi_derivative
    return intensity;
//    return 0;
}


double FDTD3D::get_reverse_mode_sim_phi_derivative(int y, int z){
    return dL_depsr[y][z];
}

void FDTD3D::init_epsilon_r_val(){
    srand (time(NULL));

    for(int j = 0; j < size_y; j++){
        for(int k = 0; k < size_z; k++) {
            if (j < 5 + pml_y || j > size_y - 6 - pml_y || k < 5 + pml_z || k > size_z - 6 - pml_z) {
                epsilon_r_val[j][k] = -10.0;
            } else {
                epsilon_r_val[j][k] = (((double) rand() / (RAND_MAX)) - 0.5) * 20;
            }
        }
    }
}

void FDTD3D::vacuum_epsilon_r_val(){
    for(int i = 0; i < size_y; i++){
        for(int k = 0; k < size_z; k++) {
            epsilon_r_val[i][k] = - 10.0;
        }
    }
}

double FDTD3D::gaussian(int x, int y, int center_x, int center_y, int start_x, int start_y, int radius){
    double x_val = (x - start_x - center_x);
    double y_val = (y - start_y - center_y);
    return exp(- (x_val*x_val + y_val*y_val) / (2 * radius * radius));
}

void FDTD3D::init_field(){
    hx_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    hy_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    hz_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dx_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dy_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dz_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    ex_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    ey_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    ez_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    imag_fft_ex = vector<vector<double>>(target_freq_y_end - target_freq_y_start,vector<double>(target_freq_z_end - target_freq_z_start,0));
    real_fft_ex = vector<vector<double>>(target_freq_y_end - target_freq_y_start,vector<double>(target_freq_z_end - target_freq_z_start,0));

    imag_fft_ey = vector<vector<double>>(target_freq_y_end - target_freq_y_start,vector<double>(target_freq_z_end - target_freq_z_start,0));
    real_fft_ey = vector<vector<double>>(target_freq_y_end - target_freq_y_start,vector<double>(target_freq_z_end - target_freq_z_start,0));

    imag_fft_ez = vector<vector<double>>(target_freq_y_end - target_freq_y_start,vector<double>(target_freq_z_end - target_freq_z_start,0));
    real_fft_ez = vector<vector<double>>(target_freq_y_end - target_freq_y_start,vector<double>(target_freq_z_end - target_freq_z_start,0));

    hy1D_arr = vector<double>(size_x, 0);
    ez1D_arr = vector<double>(size_x, 0);

    ICEx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    ICEy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    ICEz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    ICHx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    ICHy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    ICHz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    IHx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    IHy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    IHz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    IDx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    IDy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    IDz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    // Record boundary
    init_recording_field();
}

void FDTD3D::init_derivative_field(){
    dL_dHx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dHy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dHz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_dDx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dDy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dDz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_dEx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dEy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dEz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_dIHx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dIHy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dIHz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_dIDx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dIDy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dIDz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_dICEx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dICEy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dICEz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_dICHx = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dICHy = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    dL_dICHz = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    dL_deps = vector<vector<double>>(size_y,vector<double>(size_z,0));
    dL_depsr = vector<vector<double>>(size_y,vector<double>(size_z,0));
}

void FDTD3D::set_field(){
    for(int x = 0; x < size_x; x++){
        for(int y = 0; y < size_y; y++){
            for(int z = 0; z < size_z; z++){
                hx_arr[x][y][z] = 0;
                hy_arr[x][y][z] = 0;
                hz_arr[x][y][z] = 0;
                dx_arr[x][y][z] = 0;
                dy_arr[x][y][z] = 0;
                dz_arr[x][y][z] = 0;
                ex_arr[x][y][z] = 0;
                ey_arr[x][y][z] = 0;
                ez_arr[x][y][z] = 0;

                ICEx[x][y][z] = 0;
                ICEy[x][y][z] = 0;
                ICEz[x][y][z] = 0;
                ICHx[x][y][z] = 0;
                ICHy[x][y][z] = 0;
                ICHz[x][y][z] = 0;

                IDx[x][y][z] = 0;
                IDy[x][y][z] = 0;
                IDz[x][y][z] = 0;
                IHx[x][y][z] = 0;
                IHy[x][y][z] = 0;
                IHz[x][y][z] = 0;
            }
        }
    }
    for(int i = 0; i < target_freq_y_end - target_freq_y_start; i++){
        for(int j = 0; j < target_freq_z_end - target_freq_z_start; j++){
            imag_fft_ex[i][j] = 0;
            real_fft_ex[i][j] = 0;
            imag_fft_ey[i][j] = 0;
            real_fft_ey[i][j] = 0;
            imag_fft_ez[i][j] = 0;
            real_fft_ez[i][j] = 0;
        }
    }
    for(int i = 0; i < size_x; i++){
        ez1D_arr[i] = 0;
        hy1D_arr[i] = 0;
    }
    im_avg = 0;
    re_avg = 0;
}

void FDTD3D::set_derivative_field(){
    for(int x = 0; x < size_x; x++){
        for(int y = 0; y < size_y; y++){
            for(int z = 0; z < size_z; z++){
                dL_dHx[x][y][z] = 0;
                dL_dHy[x][y][z] = 0;
                dL_dHz[x][y][z] = 0;
                dL_dDx[x][y][z] = 0;
                dL_dDy[x][y][z] = 0;
                dL_dDz[x][y][z] = 0;
                dL_dEx[x][y][z] = 0;
                dL_dEy[x][y][z] = 0;
                dL_dEz[x][y][z] = 0;

                dL_dIHx[x][y][z] = 0;
                dL_dIHy[x][y][z] = 0;
                dL_dIHz[x][y][z] = 0;
                dL_dIDx[x][y][z] = 0;
                dL_dIDy[x][y][z] = 0;
                dL_dIDz[x][y][z] = 0;

                dL_dICEx[x][y][z] = 0;
                dL_dICEy[x][y][z] = 0;
                dL_dICEz[x][y][z] = 0;
                dL_dICHx[x][y][z] = 0;
                dL_dICHy[x][y][z] = 0;
                dL_dICHz[x][y][z] = 0;
            }
        }
    }
    for(int y = 0; y < size_y; y++){
        for(int z = 0; z < size_z; z++){
            dL_deps[y][z] = 0;
            dL_depsr[y][z] = 0;
        }
    }
    for(int i = 0; i < size_x; i++){
        for(int j = 0; j < size_y; j++){
            for(int k = 0; k < size_z; k++){
                invs_eps_arr[i][j][k] = pow(inv_eps_arr[i][j][k], 2);
            }
        }
    }
}

void FDTD3D::tfsf(int t){
    for(int y = firstY; y < lastY+1; y++){
        for(int z = firstZ; z < lastZ; z++){
            hy_arr[firstX - 1][y][z] -= Sc * ez1D_arr[firstX];
            hy_arr[lastX][y][z] += Sc * ez1D_arr[lastX];
        }
    }
    for(int x = firstX; x < lastX+1; x++){
        for(int z = firstZ; z < lastZ; z++){
            hx_arr[x][firstY - 1][z] += Sc * ez1D_arr[x];
            hx_arr[x][lastY][z] -= Sc * ez1D_arr[x];
        }
    }
    // update H
    for(int x = 0; x < size_x - 1; x++){
        hy1D_arr[x] = hy1D_coef1[x] * hy1D_arr[x] + hy1D_coef2[x] * (ez1D_arr[x+1] - ez1D_arr[x]);
    }

    // Set sine source
    ez1D_arr[tfsf_source_start] += sineFn(t);

    // update E
    for(int x = 1; x < size_x - 1; x++){
        // In this demonstration eps_arr[x][y][z] is the same along y and z axis.
        ez1D_arr[x] = ez1D_coef1[x] * ez1D_arr[x] + ez1D_coef2[x] * (hy1D_arr[x] - hy1D_arr[x-1]) / eps_arr[x][0][0];
    }
    ez1D_arr[0] = 0;
    ez1D_arr[size_x-1] = 0;
    hy1D_arr[0] = 0;
    hy1D_arr[size_x-1] = 0;

    for(int y = firstY; y < lastY + 1; y++){
        for(int z = firstZ; z < lastZ; z++){
            dz_arr[firstX][y][z] -= Sc * hy1D_arr[firstX-1];
            dz_arr[lastX][y][z] += Sc * hy1D_arr[lastX];
        }
    }
    for(int x = firstX; x < lastX; x++){
        for(int y = firstY; y < lastY+1; y++){
            dx_arr[x][y][firstZ] += Sc * hy1D_arr[x];
            dx_arr[x][y][lastZ] -= Sc * hy1D_arr[x];
        }
    }
}

void FDTD3D::init_recording_field(){
    int pillar_x_size = structure_end_x - structure_start_x;
    int pillar_y_size = structure_end_y - structure_start_y;
    int pillar_z_size = structure_end_z - structure_start_z;

    hx_top = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    hx_bot = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    hx_left = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hx_right = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hx_front = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hx_back = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));

    hy_top = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    hy_bot = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    hy_left = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hy_right = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hy_front = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hy_back = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));

    hz_top = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    hz_bot = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    hz_left = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hz_right = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hz_front = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    hz_back = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));

    dx_top = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    dx_bot = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    dx_left = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dx_right = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dx_front = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dx_back = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));

    dy_top = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    dy_bot = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    dy_left = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dy_right = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dy_front = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dy_back = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));

    dz_top = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    dz_bot = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_y_size,vector<double>(time_steps,0)));
    dz_left = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dz_right = vector<vector<vector<double>>>(pillar_x_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dz_front = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
    dz_back = vector<vector<vector<double>>>(pillar_y_size, vector<vector<double>>(pillar_z_size,vector<double>(time_steps,0)));
}

// Initialize forward sim matrices
void FDTD3D::init_mat(){
    init_field();
    init_derivative_field();

    epsilon_r_val = vector<vector<double>>(size_y, vector<double>(size_z,0));
    eps_var = vector<vector<double>>(size_y, vector<double>(size_z,0));
    eps_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    inv_eps_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    invs_eps_arr = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    sigma_dx = vector<double>(size_x, 0);
    sigma_dy = vector<double>(size_y, 0);
    sigma_dz = vector<double>(size_z, 0);
    sigma_hx = vector<double>(size_x, 0);
    sigma_hy = vector<double>(size_y, 0);
    sigma_hz = vector<double>(size_z, 0);

    mHx1 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHx2 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHx3 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHx4 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    mHy1 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHy2 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHy3 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHy4 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    mHz1 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHz2 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHz3 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mHz4 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    mDx1 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDx2 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDx3 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDx4 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    mDy1 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDy2 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDy3 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDy4 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    mDz1 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDz2 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDz3 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));
    mDz4 = vector<vector<vector<double>>>(size_x, vector<vector<double>>(size_y,vector<double>(size_z,0)));

    hy1D_coef1 = vector<double>(size_x, 1);
    hy1D_coef2 = vector<double>(size_x, Sc);
    ez1D_coef1 = vector<double>(size_x, 1);
    ez1D_coef2 = vector<double>(size_x, Sc);
}

void FDTD3D::init_pml(){

    for(int x = 0; x < pml_x; x++){
        sigma_dx[pml_x - x] = 0.5 * epsilon_0 / delta_t * pow((double)(x+0.5)/pml_x,3);
        sigma_dx[size_x - pml_x + x] = 0.5 * epsilon_0 / delta_t * pow((double)x/pml_x,3);
        sigma_hx[pml_x - x] = 0.5 * epsilon_0 / delta_t * pow((double)(x)/pml_x,3);
        sigma_hx[size_x - pml_x + x] = 0.5 * epsilon_0 / delta_t * pow((double)(x+0.5)/pml_x,3);
    }
    for(int y = 0; y < pml_y; y++){
        sigma_dy[pml_y - y] = 0.5 * epsilon_0 / delta_t * pow((double)(y+0.5)/pml_y, 3);
        sigma_dy[size_y - pml_y + y] = 0.5 * epsilon_0 / delta_t * pow((double)y/pml_y, 3);
        sigma_hy[pml_y - y] = 0.5 * epsilon_0 / delta_t * pow((double)(y)/pml_y, 3);
        sigma_hy[size_y - pml_y + y] = 0.5 * epsilon_0 / delta_t * pow((double)(y+0.5)/pml_y, 3);
    }

    for(int z = 0; z < pml_z; z++){
        sigma_dz[pml_z - z] = 0.5 * epsilon_0 / delta_t * pow((double)(z+0.5)/pml_z, 3);
        sigma_dz[size_z - pml_z + z] = 0.5 * epsilon_0 / delta_t * pow((double)z/pml_z, 3);
        sigma_hz[pml_z - z] = 0.5 * epsilon_0 / delta_t * pow((double)z/pml_z, 3);
        sigma_hz[size_z - pml_z + z] = 0.5 * epsilon_0 / delta_t * pow((double)(z+0.5)/pml_z, 3);
    }

    // 3D PML
    for(int x = 0; x < size_x; x++){
        for(int y = 0; y < size_y; y++) {
            for (int z = 0; z < size_z; z++) {
                double hx_A = 1.0/delta_t + (sigma_hy[y]+sigma_hz[z]) / 2.0 / epsilon_0 + sigma_hy[y] * sigma_hz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double hy_A = 1.0/delta_t + (sigma_hx[x]+sigma_hz[z]) / 2.0 / epsilon_0 + sigma_hx[x] * sigma_hz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double hz_A = 1.0/delta_t + (sigma_hx[x]+sigma_hy[y]) / 2.0 / epsilon_0 + sigma_hx[x] * sigma_hy[y] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double dx_A = 1.0/delta_t + (sigma_dy[y]+sigma_dz[z]) / 2.0 / epsilon_0 + sigma_dy[y] * sigma_dz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double dy_A = 1.0/delta_t + (sigma_dx[x]+sigma_dz[z]) / 2.0 / epsilon_0 + sigma_dx[x] * sigma_dz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double dz_A = 1.0/delta_t + (sigma_dx[x]+sigma_dy[y]) / 2.0 / epsilon_0 + sigma_dx[x] * sigma_dy[y] * delta_t / 4.0 / epsilon_0 / epsilon_0;

                double hx_B = 1.0/delta_t - (sigma_hy[y]+sigma_hz[z]) / 2.0 / epsilon_0 - sigma_hy[y] * sigma_hz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double hy_B = 1.0/delta_t - (sigma_hx[x]+sigma_hz[z]) / 2.0 / epsilon_0 - sigma_hx[x] * sigma_hz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double hz_B = 1.0/delta_t - (sigma_hx[x]+sigma_hy[y]) / 2.0 / epsilon_0 - sigma_hx[x] * sigma_hy[y] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double dx_B = 1.0/delta_t - (sigma_dy[y]+sigma_dz[z]) / 2.0 / epsilon_0 - sigma_dy[y] * sigma_dz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double dy_B = 1.0/delta_t - (sigma_dx[x]+sigma_dz[z]) / 2.0 / epsilon_0 - sigma_dx[x] * sigma_dz[z] * delta_t / 4.0 / epsilon_0 / epsilon_0;
                double dz_B = 1.0/delta_t - (sigma_dx[x]+sigma_dy[y]) / 2.0 / epsilon_0 - sigma_dx[x] * sigma_dy[y] * delta_t / 4.0 / epsilon_0 / epsilon_0;

                double hx_C = sigma_hy[y] * sigma_hz[z] / epsilon_0 / epsilon_0;
                double hy_C = sigma_hx[x] * sigma_hz[z] / epsilon_0 / epsilon_0;
                double hz_C = sigma_hx[x] * sigma_hy[y] / epsilon_0 / epsilon_0;
                double dx_C = sigma_dy[y] * sigma_dz[z] / epsilon_0 / epsilon_0;
                double dy_C = sigma_dx[x] * sigma_dz[z] / epsilon_0 / epsilon_0;
                double dz_C = sigma_dx[x] * sigma_dy[y] / epsilon_0 / epsilon_0;

                double hx_D = c * sigma_hx[x] / epsilon_0;
                double hy_D = c * sigma_hy[y] / epsilon_0;
                double hz_D = c * sigma_hz[z] / epsilon_0;
                double dx_D = c * sigma_dx[x] / epsilon_0;
                double dy_D = c * sigma_dy[y] / epsilon_0;
                double dz_D = c * sigma_dz[z] / epsilon_0;

                mHx1[x][y][z] = hx_B / hx_A;
                mHy1[x][y][z] = hy_B / hy_A;
                mHz1[x][y][z] = hz_B / hz_A;
                mDx1[x][y][z] = dx_B / dx_A;
                mDy1[x][y][z] = dy_B / dy_A;
                mDz1[x][y][z] = dz_B / dz_A;

                mHx2[x][y][z] = c / hx_A;
                mHy2[x][y][z] = c / hy_A;
                mHz2[x][y][z] = c / hz_A;
                mDx2[x][y][z] = c / dx_A;
                mDy2[x][y][z] = c / dy_A;
                mDz2[x][y][z] = c / dz_A;

                mHx3[x][y][z] = hx_D * delta_t / hx_A;
                mHy3[x][y][z] = hy_D * delta_t / hy_A;
                mHz3[x][y][z] = hz_D * delta_t / hz_A;
                mDx3[x][y][z] = dx_D * delta_t / dx_A;
                mDy3[x][y][z] = dy_D * delta_t / dy_A;
                mDz3[x][y][z] = dz_D * delta_t / dz_A;
                mHx4[x][y][z] = - hx_C * delta_t / hx_A;
                mHy4[x][y][z] = - hy_C * delta_t / hy_A;
                mHz4[x][y][z] = - hz_C * delta_t / hz_A;
                mDx4[x][y][z] = - dx_C * delta_t / dx_A;
                mDy4[x][y][z] = - dy_C * delta_t / dy_A;
                mDz4[x][y][z] = - dz_C * delta_t / dz_A;

            }
        }
    }

    // 1D PML for tfsf
    for(int i = size_x - pml_x; i < size_x; i++){
        double layer = i - (size_x - 1 - pml_x) + 0.5;
        double sigma = max_loss * pow((double)layer / pml_x, 3);
        ez1D_coef1[i] = (1 - sigma) / (1 + sigma);
        ez1D_coef2[i] = Sc / ( 1 + sigma);
        ez1D_coef1[size_x-i] = (1 - sigma) / (1 + sigma);
        ez1D_coef2[size_x-i] = Sc / ( 1 + sigma);

        layer += 0.5;
        sigma = max_loss * pow((double)layer / pml_x, 3);
        hy1D_coef1[i] = (1 - sigma) / (1 + sigma);
        hy1D_coef2[i] = Sc / (1 + sigma);
        hy1D_coef1[size_x-i] = (1 - sigma) / (1 + sigma);
        hy1D_coef2[size_x-i] = Sc / (1 + sigma);
    }
}


void FDTD3D::set_one_permittivity(int y, int z, double val) {
    this->epsilon_r_val[y][z] = val;
}

double FDTD3D::get_epsilon_r_val(int y, int z){
    return epsilon_r_val[y][z];
}

void FDTD3D::set_permittivity(){
    //glass_permittivity = 1;
    for(int i = 0; i < epsilon_r_val.size(); i++) {
        for (int k = 0; k < epsilon_r_val[0].size(); k++) {
            eps_var[i][k] = (tanh(0.5 * epsilon_r_val[i][k]) / 2.0 + 0.5) * tio2_permittivity +
                            (1.0 - (tanh(0.5 * epsilon_r_val[i][k]) / 2.0 + 0.5));
        }
    }

    for(int i = 0; i < size_x; i++){
        for(int j = 0; j < size_y; j++) {
            for (int k = 0; k < size_z; k++) {
                eps_arr[i][j][k] = 1.;
                inv_eps_arr[i][j][k] = 1.;
                invs_eps_arr[i][j][k] = 1.;
                if (i >= glass_start_x && i < glass_end_x) {
                    eps_arr[i][j][k] = glass_permittivity;
                    inv_eps_arr[i][j][k] = 1. / glass_permittivity;
                    invs_eps_arr[i][j][k] = 1. / glass_permittivity / glass_permittivity;
                } else if(i >= structure_start_x && i < structure_end_x &&
                          j >= structure_start_y && j < structure_end_y &&
                          k >= structure_start_z && k < structure_end_z){
                    eps_arr[i][j][k] = eps_var[j][k];
                    inv_eps_arr[i][j][k] = 1. / eps_var[j][k];
                    invs_eps_arr[i][j][k] = 1. / eps_var[j][k] / eps_var[j][k];
                }
            }
        }
    }
}

void FDTD3D::record_field(int cur_timestep){

    // top z
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y ++){
            hx_top[x-structure_start_x][y-structure_start_y][cur_timestep] = hx_arr[x][y][structure_end_z-1];
            hy_top[x-structure_start_x][y-structure_start_y][cur_timestep] = hy_arr[x][y][structure_end_z-1];
            hz_top[x-structure_start_x][y-structure_start_y][cur_timestep] = hz_arr[x][y][structure_end_z-1];
            dx_top[x-structure_start_x][y-structure_start_y][cur_timestep] = dx_arr[x][y][structure_end_z-1];
            dy_top[x-structure_start_x][y-structure_start_y][cur_timestep] = dy_arr[x][y][structure_end_z-1];
            dz_top[x-structure_start_x][y-structure_start_y][cur_timestep] = dz_arr[x][y][structure_end_z-1];
        }
    }
    // bot z
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int y = structure_start_y; y < structure_end_y; y ++){
            hx_bot[x-structure_start_x][y-structure_start_y][cur_timestep] = hx_arr[x][y][structure_start_z];
            hy_bot[x-structure_start_x][y-structure_start_y][cur_timestep] = hy_arr[x][y][structure_start_z];
            hz_bot[x-structure_start_x][y-structure_start_y][cur_timestep] = hz_arr[x][y][structure_start_z];
            dx_bot[x-structure_start_x][y-structure_start_y][cur_timestep] = dx_arr[x][y][structure_start_z];
            dy_bot[x-structure_start_x][y-structure_start_y][cur_timestep] = dy_arr[x][y][structure_start_z];
            dz_bot[x-structure_start_x][y-structure_start_y][cur_timestep] = dz_arr[x][y][structure_start_z];
        }
    }
    // left y
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_left[x-structure_start_x][z-structure_start_z][cur_timestep] = hx_arr[x][structure_start_y][z];
            hy_left[x-structure_start_x][z-structure_start_z][cur_timestep] = hy_arr[x][structure_start_y][z];
            hz_left[x-structure_start_x][z-structure_start_z][cur_timestep] = hz_arr[x][structure_start_y][z];
            dx_left[x-structure_start_x][z-structure_start_z][cur_timestep] = dx_arr[x][structure_start_y][z];
            dy_left[x-structure_start_x][z-structure_start_z][cur_timestep] = dy_arr[x][structure_start_y][z];
            dz_left[x-structure_start_x][z-structure_start_z][cur_timestep] = dz_arr[x][structure_start_y][z];
        }
    }
    // right y
    for(int x = structure_start_x; x < structure_end_x; x ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_right[x-structure_start_x][z-structure_start_z][cur_timestep] = hx_arr[x][structure_end_y-1][z];
            hy_right[x-structure_start_x][z-structure_start_z][cur_timestep] = hy_arr[x][structure_end_y-1][z];
            hz_right[x-structure_start_x][z-structure_start_z][cur_timestep] = hz_arr[x][structure_end_y-1][z];
            dx_right[x-structure_start_x][z-structure_start_z][cur_timestep] = dx_arr[x][structure_end_y-1][z];
            dy_right[x-structure_start_x][z-structure_start_z][cur_timestep] = dy_arr[x][structure_end_y-1][z];
            dz_right[x-structure_start_x][z-structure_start_z][cur_timestep] = dz_arr[x][structure_end_y-1][z];
        }
    }
    // front x
    for(int y = structure_start_y; y < structure_end_y; y ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_front[y-structure_start_y][z-structure_start_z][cur_timestep] = hx_arr[structure_start_x][y][z];
            hy_front[y-structure_start_y][z-structure_start_z][cur_timestep] = hy_arr[structure_start_x][y][z];
            hz_front[y-structure_start_y][z-structure_start_z][cur_timestep] = hz_arr[structure_start_x][y][z];
            dx_front[y-structure_start_y][z-structure_start_z][cur_timestep] = dx_arr[structure_start_x][y][z];
            dy_front[y-structure_start_y][z-structure_start_z][cur_timestep] = dy_arr[structure_start_x][y][z];
            dz_front[y-structure_start_y][z-structure_start_z][cur_timestep] = dz_arr[structure_start_x][y][z];
        }
    }
    // back x
    for(int y = structure_start_y; y < structure_end_y; y ++){
        for(int z = structure_start_z; z < structure_end_z; z ++){
            hx_back[y-structure_start_y][z-structure_start_z][cur_timestep] = hx_arr[structure_end_x-1][y][z];
            hy_back[y-structure_start_y][z-structure_start_z][cur_timestep] = hy_arr[structure_end_x-1][y][z];
            hz_back[y-structure_start_y][z-structure_start_z][cur_timestep] = hz_arr[structure_end_x-1][y][z];
            dx_back[y-structure_start_y][z-structure_start_z][cur_timestep] = dx_arr[structure_end_x-1][y][z];
            dy_back[y-structure_start_y][z-structure_start_z][cur_timestep] = dy_arr[structure_end_x-1][y][z];
            dz_back[y-structure_start_y][z-structure_start_z][cur_timestep] = dz_arr[structure_end_x-1][y][z];
        }
    }
}

double FDTD3D::sineFn(double time){
    // Send out sinusoidal source wave, will terminate after 600 time steps
    if(time > 600){
        return 0;
    } else {
        return sin(2 * M_PI * delta_t * time * freq);
    }
}