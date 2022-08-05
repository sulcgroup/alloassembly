/*
 * CUDAAllostericPatchySwapInteraction.cu
 *
 *  Created on: 15/jul/2020
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
__constant__ int MD_N_particle_types[1];

__constant__ int MD_N_patches[CUDAAllostericPatchySwapInteraction::MAX_SPECIES];
// patch positions
__constant__ float4 MD_base_patches[CUDAAllostericPatchySwapInteraction::MAX_SPECIES][CUDAAllostericPatchySwapInteraction::MAX_PATCHES];
// particle epsilon values
// .... huh?
__constant__ float MD_patchy_eps[CUDAAllostericPatchySwapInteraction::MAX_SPECIES * CUDAAllostericPatchySwapInteraction::MAX_SPECIES];

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

__constant__ float MD_sqr_rcut[1]; //
__constant__ float MD_sqr_rep_rcut[1]; // something something repulsive interaction
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_A_part[1], MD_B_part[1];
__constant__ float MD_spherical_attraction_strength[1], MD_spherical_E_cut[1];


#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) CUDA_FS_bond {
    int q;
    bool r_p_less_than_sigma;
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
 * I copied this from Lorenzo's code. I have no idea what any of the parameters mean;
 * I'll add them here if I can work it out
 * @param ppos the position of particle p, represented as a number4 where w can be used to find particle species
 * @param qpos the position of particle q, represented as a number4 where w can be used to find particle species
 * @param a1 first column of the rotation matrix for particle p
 * @param a2 second column of the rotation matrix for particle p
 * @param a3 third column of the rotation matrix for particle p
 * @param b1 first column of the rotation matrix for particle q
 * @param b2 second column of the rotation matrix for particle q
 * @param b3 third column of the rotation matrix for particle q
 * @param F probably a force vector? impossible to say for sure
 * @param torque
 * @param bonds
 * @param q_idx
 * @param box
 */
__device__ void _patchy_two_body_interaction(c_number4 &ppos,
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
                                             CUDABox *box) {
    // derive particle types from position values
    int ptype = get_particle_btype(ppos);
    int qtype = get_particle_btype(qpos);

    // calculate... radius? from particle positions and whatever a CUDA box is
    c_number4 r = box->minimum_image(ppos, qpos);
    // square of the radius by taking the dot product of the radius vector with itself
    // (x,y,z) . (x,y,z) = x*x + y*y + z*z = x^2 + y^2 + z^2
    c_number sqr_r = CUDA_DOT(r, r);
    // if the radius is too big for particles to interact, skip
    if(sqr_r >= MD_sqr_rcut[0]) return;

    c_number force_module = 0.;
    c_number spherical_energy = 0.;

    // centre-centre
    // repulsive energy interaction from two particles occupying the same space
    // I don't need to think about this code, so I won't
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
        spherical_energy = 4.f * (SQR(lj_part) - lj_part) + 1 - MD_spherical_attraction_strength[0];
    }

    // subtract the repulsive force to the overall force,
    // calculating repulsive force using the radius and the force scalar thingy
    // i'm surprised there isn't a faster way to calculate this
    F.x -= r.x * force_module;
    F.y -= r.y * force_module;
    F.z -= r.z * force_module;
    F.w += spherical_energy - MD_spherical_E_cut[0];

    // retrieve number of patches on each of the two particle types
    int p_N_patches = MD_N_patches[ptype];
    int q_N_patches = MD_N_patches[qtype];

    // calcualte epsilon value for interaction between the two particle types
    // ... huh?
    c_number epsilon = MD_patchy_eps[ptype + MD_N_particle_types[0] * qtype];
    if(epsilon == (c_number) 0.f) {
        // if epsilon is 0, indicating no interaction, return
        return;
    }

    // loop patches on particle p
    for(int pi = 0; pi < p_N_patches; pi++) {
        // TODO: CHECK IF PATCH IS ACTIVE HERE
        // compute... I'm going to say patch position?
        // by taking the cross product of the base patch position and [a1 a2 a3]
        c_number4 ppatch = {
                a1.x * MD_base_patches[ptype][pi].x + a2.x * MD_base_patches[ptype][pi].y + a3.x * MD_base_patches[ptype][pi].z,
                a1.y * MD_base_patches[ptype][pi].x + a2.y * MD_base_patches[ptype][pi].y + a3.y * MD_base_patches[ptype][pi].z,
                a1.z * MD_base_patches[ptype][pi].x + a2.z * MD_base_patches[ptype][pi].y + a3.z * MD_base_patches[ptype][pi].z, 0.f
        };

        // loop patches on particle q
        for(int pj = 0; pj < q_N_patches; pj++) {
            // TODO: CHECK IF PATCH IS ACTIVE HERE
            // again, going to say patch position
            // by taking the cross product of the base patch positon and [b1 b2 b3]
            c_number4 qpatch = {
                    b1.x * MD_base_patches[qtype][pj].x + b2.x * MD_base_patches[qtype][pj].y + b3.x * MD_base_patches[qtype][pj].z,
                    b1.y * MD_base_patches[qtype][pj].x + b2.y * MD_base_patches[qtype][pj].y + b3.y * MD_base_patches[qtype][pj].z,
                    b1.z * MD_base_patches[qtype][pj].x + b2.z * MD_base_patches[qtype][pj].y + b3.z * MD_base_patches[qtype][pj].z, 0.f
            };

            // distance between the two patches
            c_number4 patch_dist = {
                    r.x + qpatch.x - ppatch.x,
                    r.y + qpatch.y - ppatch.y,
                    r.z + qpatch.z - ppatch.z, 0.f
            };

            // calculate the square of the distance between the two patches by taking
            // the dot product of the distance vector with itself
            c_number patch_dist_sqr = CUDA_DOT(patch_dist, patch_dist);
            if(patch_dist_sqr < MD_sqr_patch_rcut[0]) {
                c_number r_p = sqrtf(patch_dist_sqr);
                if((r_p - MD_rcut_ss[0]) < 0.f) {
                    c_number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
                    c_number energy_part = epsilon * MD_A_part[0] * exp_part * (MD_B_part[0] / SQR(patch_dist_sqr) - 1.f);

                    c_number force_mod = epsilon * MD_A_part[0] * exp_part * (4.f * MD_B_part[0] / (SQR(patch_dist_sqr) * r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
                    c_number4 tmp_force = patch_dist * (force_mod / r_p);

                    // form bond; add to list
                    CUDA_FS_bond_list &bond_list = bonds[pi];
                    CUDA_FS_bond &my_bond = bond_list.new_bond();
                    // TODO: set patch activities based on state change

                    my_bond.force = tmp_force;
                    my_bond.force.w = energy_part;
                    my_bond.p_torque = _cross(ppatch, tmp_force);
                    my_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(b1, b2, b3, _cross(qpatch, tmp_force));
                    my_bond.q = q_idx;
                    my_bond.r_p_less_than_sigma = r_p < MD_sigma_ss[0];

                    torque -= my_bond.p_torque;
                    F.x -= tmp_force.x;
                    F.y -= tmp_force.y;
                    F.z -= tmp_force.z;
                    F.w += energy_part;
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
            for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
                CUDA_FS_bond &b2 = bond_list.bonds[bj];

                c_number curr_energy = (b1.r_p_less_than_sigma) ? 1.f : -b1.force.w;
                c_number other_energy = (b2.r_p_less_than_sigma) ? 1.f : -b2.force.w;

                // the factor 2 takes into account the fact that the pair energy is counted twice
                F.w += 2.f * MD_lambda[0] * curr_energy * other_energy;

                if(!b1.r_p_less_than_sigma) {
                    c_number factor = -MD_lambda[0] * other_energy;

                    c_number4 tmp_force = b1.force * factor;
                    tmp_force.w = 0.f;

                    F -= tmp_force;
                    LR_atomicAddXYZ(forces + b1.q, tmp_force);

                    T -= factor * b1.p_torque;
                    LR_atomicAddXYZ(torques + b1.q, b1.q_torque_ref_frame * factor);
                }

                if(!b2.r_p_less_than_sigma) {
                    c_number factor = -MD_lambda[0] * curr_energy;

                    c_number4 tmp_force = b2.force * factor;
                    tmp_force.w = 0.f;

                    F -= tmp_force;
                    LR_atomicAddXYZ(forces + b2.q, tmp_force);

                    T -= factor * b2.p_torque;
                    LR_atomicAddXYZ(torques + b2.q, b2.q_torque_ref_frame * factor);
                }
            }
        }
    }
}

/**
 * all lorenzo had to say about this is "forces + second step without lists"
 * which is actually more than he has to say about most of his functions
 * @param poss
 * @param orientations
 * @param forces
 * @param three_body_forces
 * @param torques
 * @param three_body_torques
 * @param box
 */
__global__ void PS_forces(c_number4 *poss,
                          GPU_quat *orientations,
                          c_number4 *forces,
                          c_number4 *three_body_forces,
                          c_number4 *torques,
                          c_number4 *three_body_torques,
                          CUDABox *box) {
    if(IND >= MD_N[0]) return;

    c_number4 F = forces[IND];
    c_number4 T = torques[IND];
    c_number4 ppos = poss[IND];
    GPU_quat po = orientations[IND];
    c_number4 a1, a2, a3, b1, b2, b3; // declare cols for particle orientation rotation transform matrices
    get_vectors_from_quat(po, a1, a2, a3);  // assign transform matrix cols from quaternion

    CUDA_FS_bond_list bonds[CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

    for(int j = 0; j < MD_N[0]; j++) {
        if(j != IND) {
            c_number4 qpos = poss[j];

            GPU_quat qo = orientations[j];
            get_vectors_from_quat(qo, b1, b2, b3);
            _patchy_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, j, box);
        }
    }

    _three_body(bonds, F, T, three_body_forces, three_body_torques);

    forces[IND] = F;
    torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}
/**
 * lorenzo speaketh: "forces + second step with verlet lists"
 * @param poss
 * @param orientations
 * @param forces
 * @param three_body_forces
 * @param torques
 * @param three_body_torques
 * @param matrix_neighs
 * @param c_number_neighs
 * @param box
 */
__global__ void PS_forces(c_number4 *poss,
                          GPU_quat *orientations,
                          c_number4 *forces,
                          c_number4 *three_body_forces,
                          c_number4 *torques,
                          c_number4 *three_body_torques,
                          int *matrix_neighs,
                          int *c_number_neighs,
                          CUDABox *box) {
    if(IND >= MD_N[0]) return;

    c_number4 F = forces[IND];
    c_number4 T = torques[IND];
    c_number4 ppos = poss[IND];
    GPU_quat po = orientations[IND];
    c_number4 a1, a2, a3, b1, b2, b3; // declare cols for particle orientation rotation transform matrices
    get_vectors_from_quat(po, a1, a2, a3);  // assign transform matrix cols from quaternion

    CUDA_FS_bond_list bonds[CUDAAllostericPatchySwapInteraction::MAX_PATCHES];

    int num_neighs = c_number_neighs[IND];
    for(int j = 0; j < num_neighs; j++) {
        int k_index = matrix_neighs[j * MD_N[0] + IND];

        c_number4 qpos = poss[k_index];

        GPU_quat qo = orientations[k_index];
        get_vectors_from_quat(qo, b1, b2, b3);
        _patchy_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, k_index, box);
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
    _d_three_body_forces = _d_three_body_torques = NULL;
    _step = 0;
}

CUDAAllostericPatchySwapInteraction::~CUDAAllostericPatchySwapInteraction() {
    if(_d_three_body_forces != NULL) {
        CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
    }
    if(_d_three_body_torques != NULL) {
        CUDA_SAFE_CALL(cudaFree(_d_three_body_torques));
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

    // throw an error if there are too many species
    if(_N_particle_types > MAX_SPECIES) {
        throw oxDNAException("PatchySwapInteraction: cannot simulate more than %d species. You can increase this number in the PatchySwapInteraction.h file", MAX_SPECIES);
    }

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

    int N_strands;
    std::vector<BaseParticle *> particles(N);
    AllostericPatchySwapInteraction::read_topology(&N_strands, particles);
    for(auto particle : particles) {
        delete particle;
    }

    // the following quantities are initialised by read_topology and hence have to be copied over to the GPU after its call
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_particle_types, &_N_particle_types, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, _N_patches.data(), sizeof(int) * _N_patches.size()));
    COPY_ARRAY_TO_CONSTANT(MD_patchy_eps, _patchy_eps.data(), _patchy_eps.size());

    // for each particle type
    for(int i = 0; i < _N_particle_types; i++) {
        int n_patches = _base_particle_types[i].n_patches();

        // throw an error if we've exceeded the maximum number of patches
        if(n_patches > MAX_PATCHES) {
            throw oxDNAException("PatchySwapInteraction: cannot simulate particles with more than %d patches. You can increase this number in the PatchySwapInteraction.h file", MAX_PATCHES);
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
                    patches_allosteric_flips[q][x] = _base_particle_types[i].get_state_change_result(state_change)[p];
                }
            }
            // I'm like 79% sure these values are right
            int allo_mem_count = sizeof(bool) * MAX_STATES * MAX_PATCHES;
            int allo_mem_offset = (i * MAX_PATCHES + p) * allo_mem_offset;
            CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_allosteric_controls, patches_allosteric_flips, allo_mem_count, allo_mem_offset));
        }

        // fourth argument is the offset
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4) * n_patches, i * sizeof(float4) * MAX_PATCHES));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patch_a1s, patch_a1s, sizeof(float4) * n_patches, i * sizeof(float4) * MAX_PATCHES));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patch_a2s, patch_a2s, sizeof(float4) * n_patches, i * sizeof(float4) * MAX_PATCHES));
    }
}
/**
 * Function called from CUDA thing that computes forces
 * What do the parameters mean? Some mysteries may never be solved
 * @param lists a pointer to the head of an array of lists?
 * @param d_poss pointer to the head of an array of particle positions?
 * @param d_orientations pointer to the head of an array of particle orientations?
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

    // This is some pchem nonsense
    CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
    if(_v_lists != NULL) {
        if(_v_lists->use_edge()) throw oxDNAException("CUDAAllostericPatchySwapInteraction: use_edge is unsupported");
        else {
            PS_forces
            <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
                    (d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_box);
            CUT_CHECK_ERROR("PS_forces simple_lists error");
        }
    }
    else {
        CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
        if(_no_lists != NULL) {
            PS_forces
            <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
                    (d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, d_box);
            CUT_CHECK_ERROR("PS_forces no_lists error");
        }
    }

    // add the three body contributions to the two-body forces and torques
    thrust::transform(t_forces, t_forces + N, t_three_body_forces, t_forces, thrust::plus<c_number4>());
    thrust::transform(t_torques, t_torques + N, t_three_body_torques, t_torques, thrust::plus<c_number4>());
}
