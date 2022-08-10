/*
 * CUDAAllostericPatchySwapInteraction.cu
 *
 *  Created on: 14/may/2021
 *      Author: lorenzo
 */

#include "CUDAAllostericPatchySwapInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/transform.h>

/* BEGIN CUDA */
__constant__ int MD_N[1];
__constant__ int MD_N_patch_types[1];

__constant__ int MD_N_patches[CUDAAllostericPatchySwapInteraction::MAX_SPECIES];
__constant__ int MD_patch_types[CUDAAllostericPatchySwapInteraction::MAX_SPECIES][CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

// patch a1 values (for orientation)
__constant__ float4 MD_base_patch_a1s[CUDAAllostericPatchySwapInteraction::MAX_SPECIES][CUDAAllostericPatchySwapInteraction::MAX_PATCHES];
// patch a2 values (for orientation)
__constant__ float4 MD_base_patch_a2s[CUDAAllostericPatchySwapInteraction::MAX_SPECIES][CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

// allosteric control list
/**
 * My notation here is infamously quite complecated so let's refresh:
 * if indexed as MD_allosteric_condrols[a][b][c][d],
 * a is the species that the we want to get allosteric control for
 * b is the patch on species a that we want to get allosteric control for
 * c is the current binding state of the particle that contains the patch, expressed as a binary number
 * where each digit is true if the patch at that index is bound and false otherwise
 * d is the index of the patch that is being "flipped" (bound to unbound or vice versa)
 */
__constant__ bool MD_allosteric_controls[CUDAAllostericPatchySwapInteraction::MAX_SPECIES][CUDAAllostericPatchySwapInteraction::MAX_PATCHES][CUDAAllostericPatchySwapInteraction::MAX_STATES][CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_A_part[1], MD_B_part[1];
__constant__ float MD_spherical_attraction_strength[1], MD_spherical_E_cut[1];

/// KF-related quantities
__constant__ bool MD_is_KF[1];
__constant__ int MD_patch_power[1];
__constant__ float MD_patch_pow_delta[1];
__constant__ float MD_patch_pow_cosmax[1];
__constant__ float MD_patch_angular_cutoff[1];

texture<float, 1, cudaReadModeElementType> tex_patchy_eps;
texture<float4, 1, cudaReadModeElementType> tex_base_patches;

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) CUDA_FS_bond {
    int q;
    c_number4 force;
    c_number4 p_torque;
    c_number4 q_torque_ref_frame;
};

struct __align__(16) CUDA_FS_bond_list {
    int n_bonds;
    CUDA_FS_bond bonds[CUDAAllostericPatchySwapInteraction::MAX_NEIGHS];

    __device__
    CUDA_FS_bond_list() :
            n_bonds(0) {
    }
    __device__
    CUDA_FS_bond &new_bond() {
        n_bonds++;
        if(n_bonds > CUDAAllostericPatchySwapInteraction::MAX_NEIGHS) {
            printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
            for(int i = 0; i < n_bonds; i++) {
                printf("%d ", bonds[i].q);
            }
            printf("\n");
        }
        return bonds[n_bonds - 1];
    }
};

/**
 * Computes a two-body point interaction
 * @param ppos the position of particle p, as a quaternion. particle type can be derived from w
 * @param qpos the position of particle q, as a quaternion. particle type can be derived from w
 * @param a1 the first column of particle p rotation matrix
 * @param a2 the second column of particle p rotation matrix
 * @param a3 the third column of particle p rotation matrix
 * @param b1 the first column of particle q rotation matrix
 * @param b2 the second column of particle q rotation matrix
 * @param b3 the third column of particle q rotation matrix
 * @param Fs_from_quat(qo, b1, b2, b3);

        // get particle states & activations
        // I hope k_inde
 * @param torque
 * @param bonds
 * @param q_idx
 * @param box
 * @param p_activation_states the activation states of the patches of particle p
 * @param q_activation_states the activation states of the patches of particle q
 * @param p_binding_state the binding state of particle p, where each binary digit is a patch binding state
 * @param q_binding_state the binding state of particle q, where each binary digit is a patch binding state
 */
__device__ void _patchy_point_two_body_interaction(c_number4 &ppos,
                                                   c_number4 &qpos,
                                                   c_number4 &a1,
                                                   c_number4 &a2,
                                                   c_number4 &a3,
                                                   c_number4 &b1,
                                                   c_number4 &b2,
                                                   c_number4 &b3,
                                                   c_number4 &F,
                                                   c_number4 &torque,
                                                   CUDA_FS_bond_list *bonds,
                                                   int q_idx,
                                                   CUDABox *box,
                                                   bool* p_activation_states,
                                                   bool* q_activation_states,
                                                   short &p_binding_state,
                                                   short &q_binding_state
) {
    int ptype = get_particle_btype(ppos);
    int qtype = get_particle_btype(qpos);

    c_number4 r = box->minimum_image(ppos, qpos);
    c_number sqr_r = CUDA_DOT(r, r);
    if(sqr_r >= MD_sqr_rcut[0]) return;

    c_number force_module = 0.f;
    c_number spherical_energy = 0.f;

    // center-center
    if(sqr_r >= MD_sqr_rep_rcut[0]) {
        c_number ir2 = 1.f / sqr_r;
        c_number lj_part = ir2 * ir2 * ir2;
        force_module = -24.f * MD_spherical_attraction_strength[0] * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
        spherical_energy = 4.f * MD_spherical_attraction_strength[0] * (SQR(lj_part) - lj_part);
    }
    else {
        c_number ir2 = 1.f / sqr_r;
        c_number lj_part = ir2 * ir2 * ir2;
        force_module = -24.f * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
        spherical_energy = 4.f * (SQR(lj_part) - lj_part) + 1.f - MD_spherical_attraction_strength[0];
    }

    F.x -= r.x * force_module;
    F.y -= r.y * force_module;
    F.z -= r.z * force_module;
    F.w += spherical_energy - MD_spherical_E_cut[0];

    int p_N_patches = MD_N_patches[ptype];
    int q_N_patches = MD_N_patches[qtype];

    for(int p_patch = 0; p_patch < p_N_patches; p_patch++) {
        // if patch is not active, continue
        if (!p_activation_states[p_patch]){
            continue;
        }
        c_number4 p_base_patch = tex1Dfetch(tex_base_patches, p_patch + ptype * CUDAAllostericPatchySwapInteraction::MAX_PATCHES);
        c_number4 p_patch_pos = {
                a1.x * p_base_patch.x + a2.x * p_base_patch.y + a3.x * p_base_patch.z,
                a1.y * p_base_patch.x + a2.y * p_base_patch.y + a3.y * p_base_patch.z,
                a1.z * p_base_patch.x + a2.z * p_base_patch.y + a3.z * p_base_patch.z, 0.f
        };

        for(int q_patch = 0; q_patch < q_N_patches; q_patch++) {
            // if patch is not active, continue
            if (!q_activation_states[q_patch]){
                continue;
            }
            c_number4 q_base_patch = tex1Dfetch(tex_base_patches, q_patch + qtype * CUDAAllostericPatchySwapInteraction::MAX_PATCHES);
            c_number4 q_patch_pos = {
                    b1.x * q_base_patch.x + b2.x * q_base_patch.y + b3.x * q_base_patch.z,
                    b1.y * q_base_patch.x + b2.y * q_base_patch.y + b3.y * q_base_patch.z,
                    b1.z * q_base_patch.x + b2.z * q_base_patch.y + b3.z * q_base_patch.z, 0.f
            };

            c_number4 patch_dist = {
                    r.x + q_patch_pos.x - p_patch_pos.x,
                    r.y + q_patch_pos.y - p_patch_pos.y,
                    r.z + q_patch_pos.z - p_patch_pos.z, 0.f
            };

            c_number dist = CUDA_DOT(patch_dist, patch_dist);
            if(dist < MD_sqr_patch_rcut[0]) {
                int p_patch_type = MD_patch_types[ptype][p_patch];
                int q_patch_type = MD_patch_types[qtype][q_patch];
                c_number epsilon = tex1Dfetch(tex_patchy_eps, p_patch_type + MD_N_patch_types[0] * q_patch_type);

                if(epsilon != (c_number) 0.f) {
                    c_number r_p = sqrtf(dist);
                    if((r_p - MD_rcut_ss[0]) < 0.f) {
                        c_number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
                        c_number energy_part = epsilon * MD_A_part[0] * exp_part * (MD_B_part[0] / SQR(dist) - 1.f);

                        c_number force_mod = epsilon * MD_A_part[0] * exp_part * (4.f * MD_B_part[0] / (SQR(dist) * r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
                        c_number4 tmp_force = patch_dist * (force_mod / r_p);

                        c_number4 p_torque = _cross(p_patch_pos, tmp_force);

                        torque -= p_torque;
                        F.x -= tmp_force.x;
                        F.y -= tmp_force.y;
                        F.z -= tmp_force.z;
                        F.w += energy_part;

                        CUDA_FS_bond &my_bond = bonds[p_patch].new_bond();

                        my_bond.q = q_idx;

                        if(r_p > MD_sigma_ss[0]) {
                            my_bond.force = tmp_force;
                            my_bond.force.w = -energy_part;
                            my_bond.p_torque = p_torque;
                            my_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(b1, b2, b3, _cross(q_patch_pos, tmp_force));
                        }
                        else {
                            my_bond.force.w = epsilon;
                        }
                        // update both particles binding states
                        // check for flips
//                        for (int i = 0; i < CUDAAllostericPatchySwapInteraction::MAX_PATCHES; i++){
//                            // if particle p type has a patch i and that patch is flagged to be flipped by this transition
//                            if (i < MD_N_patches[ptype] && MD_allosteric_controls[ptype][i][p_binding_state][p_patch]){
//                                p_activation_states[i] = !p_activation_states[i]; // flip activation state
//                            }
//                            // if particle q type has a patch i and that patch is flagged to be flipped by this transition
//                            if (i < MD_N_patches[qtype] && MD_allosteric_controls[qtype][i][q_binding_state][q_patch]){
//                                q_activation_states[i] = !q_activation_states[i]; // flip activation state
//                            }
//                        }
                        // Bitwise xor with p_patch should flip binding state. TODO: double-check
                        p_binding_state = p_binding_state ^ (1 << p_patch);
                        q_binding_state = q_binding_state ^ (1 << q_patch);
                    }
                }
            }
        }
    }
}


/**
 * Computes a two-body kern-frankel interaction
 * @param ppos the position of particle p, as a quaternion. particle type can be derived from w
 * @param qpos the position of particle q, as a quaternion. particle type can be derived from w
 * @param a1 the first column of particle p rotation matrix
 * @param a2 the second column of particle p rotation matrix
 * @param a3 the third column of particle p rotation matrix
 * @param b1 the first column of particle q rotation matrix
 * @param b2 the second column of particle q rotation matrix
 * @param b3 the third column of particle q rotation matrix
 * @param F
 * @param torque
 * @param bonds
 * @param q_idx
 * @param box
 * @param p_activation_states the activation states of the patches of particle p
 * @param q_activation_states the activation states of the patches of particle q
 * @param p_binding_state the binding state of particle p, where each binary digit is a patch binding state
 * @param q_binding_state the binding state of particle q, where each binary digit is a patch binding state
 */
__device__ void _patchy_KF_two_body_interaction(c_number4 &ppos,
                                                c_number4 &qpos,
                                                c_number4 &a1,
                                                c_number4 &a2,
                                                c_number4 &a3,
                                                c_number4 &b1,
                                                c_number4 &b2,
                                                c_number4 &b3,
                                                c_number4 &F,
                                                c_number4 &torque,
                                                CUDA_FS_bond_list *bonds,
                                                int q_idx,
                                                CUDABox *box,
                                                bool* p_activation_states,
                                                bool* q_activation_states,
                                                short &p_binding_state,
                                                short &q_binding_state) {
    int ptype = get_particle_btype(ppos);
    int qtype = get_particle_btype(qpos);

    c_number4 r = box->minimum_image(ppos, qpos);
    c_number sqr_r = CUDA_DOT(r, r);
    if(sqr_r >= MD_sqr_rcut[0]) return;

    c_number force_module = 0.f;
    c_number spherical_energy = 0.f;

    // centre-centre
    if(sqr_r >= MD_sqr_rep_rcut[0]) {
        c_number ir2 = 1.f / sqr_r;
        c_number lj_part = ir2 * ir2 * ir2;
        force_module = -24.f * MD_spherical_attraction_strength[0] * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
        spherical_energy = 4.f * MD_spherical_attraction_strength[0] * (SQR(lj_part) - lj_part);
    }
    else {
        c_number ir2 = 1.f / sqr_r;
        c_number lj_part = ir2 * ir2 * ir2;
        force_module = -24.f * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
        spherical_energy = 4.f * (SQR(lj_part) - lj_part) + 1.f - MD_spherical_attraction_strength[0];
    }

    F.x -= r.x * force_module;
    F.y -= r.y * force_module;
    F.z -= r.z * force_module;
    F.w += spherical_energy - MD_spherical_E_cut[0];

    // patch-patch part
    c_number rmod = sqrtf(sqr_r);
    c_number4 r_versor = r / rmod;

    c_number dist_surf = rmod - 1.f;
    c_number dist_surf_sqr = SQR(dist_surf);
    c_number r8b10 = SQR(SQR(dist_surf_sqr)) / MD_patch_pow_delta[0];
    c_number exp_part = -1.001f * expf(-0.5f * r8b10 * dist_surf_sqr);

    int p_N_patches = MD_N_patches[ptype];
    int q_N_patches = MD_N_patches[qtype];

    for(int p_patch = 0; p_patch < p_N_patches; p_patch++) {
        if (!p_activation_states[p_patch]){
            continue;
        }
        c_number4 p_base_patch = tex1Dfetch(tex_base_patches, p_patch + ptype * CUDAAllostericPatchySwapInteraction::MAX_PATCHES);
        c_number4 p_patch_pos = {
                a1.x * p_base_patch.x + a2.x * p_base_patch.y + a3.x * p_base_patch.z,
                a1.y * p_base_patch.x + a2.y * p_base_patch.y + a3.y * p_base_patch.z,
                a1.z * p_base_patch.x + a2.z * p_base_patch.y + a3.z * p_base_patch.z, 0.f
        };
        p_patch_pos *= 2.f;

        c_number cospr = CUDA_DOT(p_patch_pos, r_versor);
        c_number cospr_minus_one = cospr - 1.f;
        if(cospr_minus_one < MD_patch_angular_cutoff[0]) {

            // what follows is a slightly faster way of doing (cospr - 1)^(MD_patch_power - 1) than a regular loop
            c_number part = SQR(cospr_minus_one);
            c_number cospr_base = cospr_minus_one;
            for(int i = 0; i < MD_patch_power[0] / 2 - 1; i++) {
                cospr_base *= part;
            }

            // we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
            c_number cospr_part = cospr_base * cospr_minus_one;
            c_number p_mod = expf(-cospr_part / (2.f * MD_patch_pow_cosmax[0]));

            for(int q_patch = 0; q_patch < q_N_patches; q_patch++) {
                if (!q_activation_states[q_patch]){
                    continue;
                }
                c_number4 q_base_patch = tex1Dfetch(tex_base_patches, q_patch + qtype * CUDAAllostericPatchySwapInteraction::MAX_PATCHES);
                c_number4 q_patch_pos = {
                        b1.x * q_base_patch.x + b2.x * q_base_patch.y + b3.x * q_base_patch.z,
                        b1.y * q_base_patch.x + b2.y * q_base_patch.y + b3.y * q_base_patch.z,
                        b1.z * q_base_patch.x + b2.z * q_base_patch.y + b3.z * q_base_patch.z, 0.f
                };
                q_patch_pos *= 2.f;

                c_number cosqr = -CUDA_DOT(q_patch_pos, r_versor);
                c_number cosqr_minus_one = cosqr - 1.f;
                if(cosqr_minus_one < MD_patch_angular_cutoff[0]) {
                    int p_patch_type = MD_patch_types[ptype][p_patch];
                    int q_patch_type = MD_patch_types[qtype][q_patch];
                    c_number epsilon = tex1Dfetch(tex_patchy_eps, p_patch_type + MD_N_patch_types[0] * q_patch_type);

                    if(epsilon != 0.f) {
                        part = SQR(cosqr_minus_one);
                        c_number cosqr_base = cosqr_minus_one;
                        for(int i = 0; i < MD_patch_power[0] / 2 - 1; i++) {
                            cosqr_base *= part;
                        }

                        c_number cosqr_part = cosqr_base * cosqr_minus_one;
                        c_number q_mod = expf(-cosqr_part / (2.f * MD_patch_pow_cosmax[0]));

                        c_number energy_part = exp_part * p_mod * q_mod;

                        // radial part
                        c_number4 radial_force = r_versor * (p_mod * q_mod * 5.f * (rmod - 1.f) * exp_part * r8b10);

                        // angular p part
                        c_number der_p = exp_part * q_mod * (0.5f * MD_patch_power[0] * p_mod * cospr_base / MD_patch_pow_cosmax[0]);
                        c_number4 p_ortho = p_patch_pos - cospr * r_versor;
                        c_number4 angular_force = p_ortho * (der_p / rmod);

                        // angular q part
                        c_number der_q = exp_part * p_mod * (-0.5f * MD_patch_power[0] * q_mod * cosqr_base / MD_patch_pow_cosmax[0]);
                        c_number4 q_ortho = q_patch_pos + cosqr * r_versor;
                        angular_force += q_ortho * (der_q / rmod);

                        c_number4 p_torque = _cross(r_versor, p_patch_pos) * der_p;
                        c_number4 q_torque = _cross(q_patch_pos, r_versor) * der_q;

                        c_number4 tot_force = radial_force + angular_force;

                        torque -= p_torque;
                        F.x -= tot_force.x;
                        F.y -= tot_force.y;
                        F.z -= tot_force.z;
                        F.w += energy_part;

                        if(energy_part < 0.f) {
                            CUDA_FS_bond &my_bond = bonds[p_patch].new_bond();

                            my_bond.force = (dist_surf < MD_sigma_ss[0]) ? angular_force : tot_force;
                            my_bond.force.w = (dist_surf < MD_sigma_ss[0]) ? epsilon * p_mod * q_mod : -energy_part;
                            my_bond.p_torque = p_torque;
                            my_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(b1, b2, b3, q_torque);
                            my_bond.q = q_idx;
                        }

//                        // update both particles binding states
//                        // check for flips
//                        for (int i = 0; i < CUDAAllostericPatchySwapInteraction::MAX_PATCHES; i++){
//                            // if particle p type has a patch i and that patch is flagged to be flipped by this transition
//                            if (i < MD_N_patches[ptype] && MD_allosteric_controls[ptype][i][p_binding_state][p_patch]){
//                                p_activation_states[i] = !p_activation_states[i]; // flip activation state
//                            }
//                            // if particle q type has a patch i and that patch is flagged to be flipped by this transition
//                            if (i < MD_N_patches[qtype] && MD_allosteric_controls[qtype][i][q_binding_state][q_patch]){
//                                q_activation_states[i] = !q_activation_states[i]; // flip activation state
//                            }
//                        }
//                        // Bitwise xor with p_patch should flip binding state. TODO: double-check
//                        p_binding_state = p_binding_state ^ (1 << p_patch);
//                        q_binding_state = q_binding_state ^ (1 << q_patch);

                    }

                }
            }
        }
    }
}

__device__ void _three_body(CUDA_FS_bond_list *bonds, c_number4 &F, c_number4 &T, c_number4 *forces, c_number4 *torques) {
    for(int pi = 0; pi < CUDAAllostericPatchySwapInteraction::MAX_PATCHES; pi++) {
        CUDA_FS_bond_list &bond_list = bonds[pi];

        for(int bi = 0; bi < bond_list.n_bonds; bi++) {
            CUDA_FS_bond &b1 = bond_list.bonds[bi];
            c_number curr_energy = b1.force.w;

            for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
                CUDA_FS_bond &b2 = bond_list.bonds[bj];
                c_number other_energy = b2.force.w;

                // the factor 2 takes into account the fact that the total pair energy is always counted twice
                F.w += 2.f * MD_lambda[0] * curr_energy * other_energy;

                // b1
                c_number factor = -MD_lambda[0] * other_energy;

                c_number4 tmp_force = b1.force * factor;
                tmp_force.w = 0.f;

                F -= tmp_force;
                LR_atomicAddXYZ(forces + b1.q, tmp_force);

                T -= factor * b1.p_torque;
                LR_atomicAddXYZ(torques + b1.q, b1.q_torque_ref_frame * factor);

                // b2
                factor = -MD_lambda[0] * curr_energy;

                tmp_force = b2.force * factor;
                tmp_force.w = 0.f;

                F -= tmp_force;
                LR_atomicAddXYZ(forces + b2.q, tmp_force);

                T -= factor * b2.p_torque;
                LR_atomicAddXYZ(torques + b2.q, b2.q_torque_ref_frame * factor);
            }
        }
    }
}

/** forces + second step without lists
 * @param poss
 * @param orientations
 * @param forces
 * @param three_body_forces
 * @param torques
 * @param three_body_torques
 * @param box
 * @param p_activations
 * @param q_activations
 * @param p_state
 * @param q_state
 */
__global__ void DPS_forces(c_number4 *poss,
                           GPU_quat *orientations,
                           c_number4 *forces,
                           c_number4 *three_body_forces,
                           c_number4 *torques,
                           c_number4 *three_body_torques,
                           CUDABox *box,
                           bool* activation_states,
                           short* particle_states
) {
    if(IND >= MD_N[0]) return;

    c_number4 F = forces[IND];
    c_number4 T = torques[IND];
    c_number4 ppos = poss[IND];
    GPU_quat po = orientations[IND];
    c_number4 a1, a2, a3, b1, b2, b3;
    get_vectors_from_quat(po, a1, a2, a3);
    bool* p_activations = &activation_states[IND * CUDAAllostericPatchySwapInteraction::MAX_PATCHES];
    short& p_state = particle_states[IND];

    CUDA_FS_bond_list bonds[CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

    for(int j = 0; j < MD_N[0]; j++) {
        if(j != IND) {
            c_number4 qpos = poss[j];

            GPU_quat qo = orientations[j];
            get_vectors_from_quat(qo, b1, b2, b3);

            // get particle states & activations
            // I hope j is the correct var
            bool* q_activations = &activation_states[j];
            short& q_state = particle_states[j];

            if(MD_is_KF[0]) {
                _patchy_KF_two_body_interaction(ppos,
                                                qpos,
                                                a1,
                                                a2,
                                                a3,
                                                b1,
                                                b2,
                                                b3,
                                                F,
                                                T,
                                                bonds,
                                                j,
                                                box,
                                                p_activations,
                                                q_activations,
                                                p_state,
                                                q_state);
            }
            else {
                _patchy_point_two_body_interaction(ppos,
                                                   qpos,
                                                   a1,
                                                   a2,
                                                   a3,
                                                   b1,
                                                   b2,
                                                   b3,
                                                   F,
                                                   T,
                                                   bonds,
                                                   j,
                                                   box,
                                                   p_activations,
                                                   q_activations,
                                                   p_state,
                                                   q_state);
            }
        }
    }

    _three_body(bonds, F, T, three_body_forces, three_body_torques);

    forces[IND] = F;
    torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

/** forces + second step with verlet lists
 *
 * @param poss
 * @param orientations
 * @param forces
 * @param three_body_forces
 * @param torques
 * @param three_body_torques
 * @param matrix_neighs
 * @param c_number_neighs
 * @param box
 * @param p_activations
 * @param q_activations
 * @param p_state
 * @param q_state
 */
__global__ void DPS_forces(c_number4 *poss,
                           GPU_quat *orientations,
                           c_number4 *forces,
                           c_number4 *three_body_forces,
                           c_number4 *torques,
                           c_number4 *three_body_torques,
                           int *matrix_neighs,
                           int *c_number_neighs,
                           CUDABox *box,
                           bool* activation_states,
                           short* particle_states) {
    if(IND >= MD_N[0]) return;

    c_number4 F = forces[IND];
    c_number4 T = torques[IND];
    c_number4 ppos = poss[IND];
    GPU_quat po = orientations[IND];
    c_number4 a1, a2, a3, b1, b2, b3;
    get_vectors_from_quat(po, a1, a2, a3);
    // is IND the correct indexer here? I hope so!
    bool* p_activations = &activation_states[IND * CUDAAllostericPatchySwapInteraction::MAX_PATCHES];
    short& p_state = particle_states[IND];

    CUDA_FS_bond_list bonds[CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

    int num_neighs = c_number_neighs[IND];
    for(int j = 0; j < num_neighs; j++) {
        int k_index = matrix_neighs[j * MD_N[0] + IND];

        c_number4 qpos = poss[k_index];

        GPU_quat qo = orientations[k_index];
        get_vectors_from_quat(qo, b1, b2, b3);

        // get particle states & activations
        // I hope k_index is the correct var
        bool* q_activations = &activation_states[k_index * CUDAAllostericPatchySwapInteraction::MAX_PATCHES];
        short& q_state = particle_states[k_index];
        if(MD_is_KF[0]) {
            _patchy_KF_two_body_interaction(ppos,
                                            qpos,
                                            a1,
                                            a2,
                                            a3,
                                            b1,
                                            b2,
                                            b3,
                                            F,
                                            T,
                                            bonds,
                                            k_index,
                                            box,
                                            p_activations,
                                            q_activations,
                                            p_state,
                                            q_state);
        }
        else {
            _patchy_point_two_body_interaction(ppos,
                                               qpos,
                                               a1,
                                               a2,
                                               a3,
                                               b1,
                                               b2,
                                               b3,
                                               F,
                                               T,
                                               bonds,
                                               k_index,
                                               box,
                                               p_activations,
                                               q_activations,
                                               p_state,
                                               q_state);
        }
    }

    _three_body(bonds, F, T, three_body_forces, three_body_torques);

    forces[IND] = F;
    torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

/* END CUDA PART */

#define HALF_ISQRT3 0.28867513459481292f

CUDAAllostericPatchySwapInteraction::CUDAAllostericPatchySwapInteraction() :
        CUDABaseInteraction(),
        AllostericPatchySwapInteraction() {
    _step = 0;
}

CUDAAllostericPatchySwapInteraction::~CUDAAllostericPatchySwapInteraction() {
    if(_d_three_body_forces != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
    }

    if(_d_three_body_torques != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_three_body_torques));
    }

    if(_d_patchy_eps != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_patchy_eps));
        cudaUnbindTexture(tex_patchy_eps);
    }

    if(_d_base_patches != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_base_patches));
        cudaUnbindTexture(tex_base_patches);
    }

    if (particle_states != nullptr) {
        CUDA_SAFE_CALL(cudaFree(particle_states));
    }

    if (activation_states != nullptr) {
        CUDA_SAFE_CALL(cudaFree(activation_states));
    }
}

void CUDAAllostericPatchySwapInteraction::get_settings(input_file &inp) {
    AllostericPatchySwapInteraction::get_settings(inp);

    int sort_every = 0;
    getInputInt(&inp, "CUDA_sort_every", &sort_every, 0);
}

void CUDAAllostericPatchySwapInteraction::cuda_init(c_number box_side, int N) {
    CUDABaseInteraction::cuda_init(box_side, N);
    AllostericPatchySwapInteraction::init();

    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_forces, N * sizeof(c_number4)));
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_torques, N * sizeof(c_number4)));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

    COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, _sqr_rcut);
    COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, _sqr_rep_rcut);
    COPY_NUMBER_TO_FLOAT(MD_sqr_patch_rcut, _sqr_patch_rcut);
    COPY_NUMBER_TO_FLOAT(MD_sigma_ss, _sigma_ss);
    COPY_NUMBER_TO_FLOAT(MD_rcut_ss, _rcut_ss);
    COPY_NUMBER_TO_FLOAT(MD_lambda, _lambda);
    COPY_NUMBER_TO_FLOAT(MD_A_part, _A_part);
    COPY_NUMBER_TO_FLOAT(MD_B_part, _B_part);
    COPY_NUMBER_TO_FLOAT(MD_spherical_E_cut, _spherical_E_cut);
    COPY_NUMBER_TO_FLOAT(MD_spherical_attraction_strength, _spherical_attraction_strength);

    // KF stuff
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_is_KF, &_is_KF, sizeof(bool)));

    if(_is_KF) {
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_power, &_patch_power, sizeof(int)));
        COPY_NUMBER_TO_FLOAT(MD_patch_pow_delta, _patch_pow_delta);
        COPY_NUMBER_TO_FLOAT(MD_patch_pow_cosmax, _patch_pow_cosmax);
        COPY_NUMBER_TO_FLOAT(MD_patch_angular_cutoff, _patch_angular_cutoff);
    }

    int N_strands;
    std::vector<BaseParticle *> particles(N);
    AllostericPatchySwapInteraction::read_topology(&N_strands, particles);

    // init particle state vars
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&particle_states, sizeof(short) * N));
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&activation_states, sizeof(bool) * N * MAX_PATCHES));
    bool* unbound_state = new bool[MAX_PATCHES];
    short* particle_states_empty = new short[N];
    bool* all_activation_states = new bool[MAX_PATCHES * N];
    std::fill(unbound_state, unbound_state + MAX_PATCHES, false);
    std::fill(particle_states_empty, particle_states_empty + N, 0);
    for (int i = 0; i < N; i++){
        AllostericPatchyParticle* pp = dynamic_cast<AllostericPatchyParticle*>(particles[i]);
        for (int x = 0; x < pp->n_patches(); x++){
            all_activation_states[i * MAX_PATCHES + x] = pp->patch_status(unbound_state, x);
        }
    }
    CUDA_SAFE_CALL(cudaMemcpy(activation_states, all_activation_states, sizeof(bool) * N * MAX_PATCHES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(particle_states, particle_states_empty, sizeof (short) * N, cudaMemcpyHostToDevice));
    delete[] unbound_state;
    delete[] particle_states_empty;
    delete[] all_activation_states;

    for(auto particle : particles) {
        delete particle;
    }

    int N_species = _N_patches.size();
    if(N_species > MAX_SPECIES) {
        throw oxDNAException("PatchySwapInteraction: cannot simulate more than %d species. You can increase this number in the PatchySwapInteraction.h file", MAX_SPECIES);
    }

    // the following quantities are initialised by read_topology and hence have to be copied over to the GPU after its call
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patch_types, &_N_patch_types, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, _N_patches.data(), sizeof(int) * N_species));

    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_patchy_eps, _patchy_eps.size() * sizeof(float)));
    std::vector<float> h_patchy_eps(_patchy_eps.begin(), _patchy_eps.end());
    CUDA_SAFE_CALL(cudaMemcpy(_d_patchy_eps, h_patchy_eps.data(), _patchy_eps.size() * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_patchy_eps, _d_patchy_eps, _patchy_eps.size() * sizeof(float)));

    int N_base_patches = MAX_PATCHES * N_species;
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_base_patches, N_base_patches * sizeof(float4)));
    std::vector<float4> h_base_patches(N_base_patches, make_float4(0., 0., 0., 0.));
    for(uint ns = 0; ns < N_species; ns++) {
        AllostericPatchyParticle& particle_type = _base_particle_types[ns];
        for(uint np = 0; np < particle_type.n_patches(); np++) {
            AllostericPatch& patch = particle_type.patches[np];
            float4 bp_f4 = make_float4(patch.position().x, patch.position().y, patch.position().z, 0.);
            h_base_patches[ns * MAX_PATCHES + np] = bp_f4;
        }
    }

    CUDA_SAFE_CALL(cudaMemcpy(_d_base_patches, h_base_patches.data(), N_base_patches * sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_base_patches, _d_base_patches, N_base_patches * sizeof(float4)));

    for(int i = 0; i < N_species; i++) {
        int n_patches = _N_patches[i];

        if(n_patches > MAX_PATCHES) {
            throw oxDNAException("CUDAAllostericPatchySwapInteraction: cannot simulate particles with more than %d patches. You can increase this number in the AllostericPatchySwapInteraction.h file", MAX_PATCHES);
        }

        int patch_types[MAX_PATCHES];
        for(int p = 0; p < n_patches; p++) {
            patch_types[p] = _base_particle_types[i].patches[p].color();
        }

        float4 base_patches[MAX_PATCHES];
        float4 patch_a1s[MAX_PATCHES];
        float4 patch_a2s[MAX_PATCHES];
        // allocate memory for patch position
        for(int p = 0; p < n_patches; p++) {
            // patch position
            LR_vector patch_position = _base_particle_types[i].patches[p].position();
            base_patches[p] = make_c_number4(patch_position.x, patch_position.y, patch_position.z, 0);
            // patch orientation
            LR_vector a1 = _base_particle_types[i].patches[p].a1();
            LR_vector a2 = _base_particle_types[i].patches[i].a2();
            patch_a1s[p] = make_c_number4(a1.x, a1.y, a1.z, 0);
            patch_a2s[p] = make_c_number4(a2.x, a2.y, a2.z, 0);

            // time to deal with allostery!
            bool patches_allosteric_flips[MAX_STATES][MAX_PATCHES];

            bool state[MAX_PATCHES];
            for (int q = 0; q < MAX_STATES; q++){
                // each unique state can be expressed as an MAX_STATES-digit binary number where
                // each digit is a patch binding state

                // first decode state
                int n = q;
                for (int x = 0; x < MAX_PATCHES; x++){
                    state[x] = n & 1;
                    n /= 2;
                }

                // encode flip value for each patch x in relation to q
                for (int x = 0; x < MAX_PATCHES; x++) {
                    // get the particle state change originating at `state` when patch `x` is flipped
                    ParticleStateChange state_change(state, MAX_PATCHES, x);
                    // get the state change result, specifically the effect on patch p
                    patches_allosteric_flips[q][x] = _base_particle_types[i].get_state_change_result(state_change, p);
                }
            }
            // I'm like 79% sure these values are right
            int allo_mem_count = sizeof(bool) * MAX_STATES * MAX_PATCHES;
            int allo_mem_offset = (i * MAX_PATCHES + p) * allo_mem_offset;
            CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_allosteric_controls, patches_allosteric_flips, allo_mem_count, allo_mem_offset));
        }

        // fourth argument is the offset
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_types, patch_types, sizeof(int) * n_patches, i * sizeof(int) * MAX_PATCHES));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patch_a1s, patch_a1s, sizeof(float4) * n_patches, i * sizeof(float4) * MAX_PATCHES));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patch_a2s, patch_a2s, sizeof(float4) * n_patches, i * sizeof(float4) * MAX_PATCHES));
    }
}

/**
 *
 * @param lists list of particles
 * @param d_poss probably an array of particle positions?
 * @param d_orientations probably an array of particle orientations?
 * @param d_forces
 * @param d_torques
 * @param d_bonds
 * @param d_box
 */
void CUDAAllostericPatchySwapInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
    int N = CUDABaseInteraction::_N;
    thrust::device_ptr < c_number4 > t_forces = thrust::device_pointer_cast(d_forces);
    thrust::device_ptr < c_number4 > t_torques = thrust::device_pointer_cast(d_torques);
    thrust::device_ptr < c_number4 > t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
    thrust::device_ptr < c_number4 > t_three_body_torques = thrust::device_pointer_cast(_d_three_body_torques);
    thrust::fill_n(t_three_body_forces, N, make_c_number4(0, 0, 0, 0));
    thrust::fill_n(t_three_body_torques, N, make_c_number4(0, 0, 0, 0));

    CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
    if(_v_lists != NULL) {
        if(_v_lists->use_edge()) throw oxDNAException("CUDAAllostericPatchySwapInteraction: use_edge is unsupported");
        else {
            DPS_forces
            <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
                    (d_poss,
                     d_orientations,
                     d_forces,
                     _d_three_body_forces,
                     d_torques,
                     _d_three_body_torques,
                     _v_lists->d_matrix_neighs,
                     _v_lists->d_number_neighs,
                     d_box,
                     this->activation_states,
                     this->particle_states);
            CUT_CHECK_ERROR("DPS_forces simple_lists error");
        }
    }
    else {
        CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
        if(_no_lists != NULL) {
            DPS_forces
            <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
                    (d_poss,
                     d_orientations,
                     d_forces,
                     _d_three_body_forces,
                     d_torques,
                     _d_three_body_torques,
                     d_box,
                     this->activation_states,
                     this->particle_states);
            CUT_CHECK_ERROR("DPS_forces no_lists error");
        }
    }

    // add the three body contributions to the two-body forces and torques
    thrust::transform(t_forces, t_forces + N, t_three_body_forces, t_forces, thrust::plus<c_number4>());
    thrust::transform(t_torques, t_torques + N, t_three_body_torques, t_torques, thrust::plus<c_number4>());
}
