//
// Created by josh on 7/8/22.
// Based on Lorenzo Rogivatti's CUDADetailedPatchySwapInteraction.h
//

#ifndef OXDNA_CUDAALLOSTERICPATCHYINTERACTION_H
#define OXDNA_CUDAALLOSTERICPATCHYINTERACTION_H


#include <curand_kernel.h>
#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "AllostericPatchyInteraction.h"

struct CUDA_FS_bonding_pattern;
struct swap_event;

/**
 * @brief CUDA implementation of the {@link PatchySwapInteraction}.
 */

class CUDAAllostericPatchyInteraction: public CUDABaseInteraction, public AllostericPatchyInteraction {
protected:
    // rng (for state transitions)
    curandState *_d_rand_state;

    c_number4 *_d_three_body_forces, *_d_three_body_torques;
    float *_d_patchy_eps = nullptr;
    float4 *_d_base_patches = nullptr;

    // particle internal state vars - req'd for allostery
    // the binding states are an array (length num particles) of shorts where each bit in
    // each short is a patch binding state - 0 for unbound, 1 for bound
    // it's useful to store the binding states like this so that the binding state short
    // can be used to query the binding state transition map (MD_allosteric_controls)
    unsigned int *_cu_particle_states = nullptr;

    // flattened 3d array of particle activations as a function of type and state
    // first axis is particle type, second is state index, third is patch index
    // item at _particle_activation_map[type_idx * NUM_STATES * NUM_PATCHEs + state_idx * NUM_PATCHES + patch_idx
    // will be true if particle type type_idx patch number patch_idx is active when the particle is in state state_idx
    // this is effectively a method for mapping particle.patches[p]._activation_var to CUDA
    bool* _cu_particle_activation_map = nullptr;

    // flattened 3d array that expidites state transitions
    // first axis is particle type, second axis is state, third is a table that the program rolls on to determine
    // how the particle should change state
    unsigned int* _cu_state_transition_map = nullptr;

    // temporary variable to store patch activations
    bool* _patch_activations = nullptr;

    llint _step;
public:
//    static const int MAX_PATCHES = 6;
    static const int MAX_PATCHES = 24;
    static const int MAX_SPECIES = 5;
    static const int MAX_NEIGHS = 20;
    static const int MAX_STATE_VARS = 16;
    static const int MAX_STATES = 1 << MAX_STATE_VARS; // 2 ^ MAX_PATCHES

    CUDAAllostericPatchyInteraction();

    virtual ~CUDAAllostericPatchyInteraction();

    void sync_host();
    void sync_GPU();

    void get_settings(input_file &inp);

    void cuda_init(c_number box_side, int N);

    c_number get_cuda_rcut() {
        return this->get_rcut();
    }

    void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces,
                        c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);

    virtual number
    pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    virtual void begin_energy_computation() override;

    size_t getActivationsArrayLength() const {
        return sizeof(bool) * cudaParticleMemoryCount() * MAX_PATCHES;
    }

    size_t getBindingStateArrayLength() const {
        return sizeof(int) * cudaParticleMemoryCount();
    }

    /**
     * @return the number of particles that have had space allocated for them on the GPU
     * this may be greater than the actual number of particles
     */
    int cudaParticleMemoryCount() const {
        return CUDABaseInteraction::_N;
    }

    /**
     * @return the actual number of particles in the simulation, according to the CPU-side
     * memory
     */
    int realNumParticles() const {
        return AllostericPatchyInteraction::_N;
    }

    /**
     * @return the largest particle state size in the simulation
     */
    int maxStateSize() const {
        int iMaxSize = 0;
        for (int i = 0; i < _N_particle_types; i++){
            iMaxSize = std::max(iMaxSize, _base_particle_types[i].state_size());
        }
        return iMaxSize;
    }

};
extern "C" BaseInteraction *make_CUDAAllostericPatchyInteraction() {
    return new CUDAAllostericPatchyInteraction();
}

#endif //OXDNA_CUDAALLOSTERICPATCHYINTERACTION_H
