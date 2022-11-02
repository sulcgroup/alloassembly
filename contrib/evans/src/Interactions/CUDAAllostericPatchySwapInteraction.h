//
// Created by josh on 7/8/22.
// Based on Lorenzo Rogivatti's CUDADetailedPatchySwapInteraction.h
//

#ifndef OXDNA_CUDAALLOSTERICPATCHYSWAPINTERACTION_H
#define OXDNA_CUDAALLOSTERICPATCHYSWAPINTERACTION_H


#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "AllostericPatchySwapInteraction.h"

struct CUDA_FS_bonding_pattern;
struct swap_event;

/**
 * @brief CUDA implementation of the {@link PatchySwapInteraction}.
 */

class CUDAAllostericPatchySwapInteraction: public CUDABaseInteraction, public AllostericPatchySwapInteraction {
protected:
    c_number4 *_d_three_body_forces, *_d_three_body_torques;
    float *_d_patchy_eps = nullptr;
    float4 *_d_base_patches = nullptr;

    // particle internal state vars - req'd for allostery
    // the binding states are an array (length num particles) of shorts where each bit in
    // each short is a patch binding state - 0 for unbound, 1 for bound
    // it's useful to store the binding states like this so that the binding state short
    // can be used to query the binding state transition map (MD_allosteric_controls)
    unsigned int *_particle_binding_states = nullptr;
    // we store patch activations as an array of booleans with length
    // (num particles) x (max number of patches per particle) because it's simpler
    bool *_particle_activations = nullptr;

    llint _step;
public:
    static const int MAX_PATCHES = 5;
    static const int MAX_SPECIES = 5;
    static const int MAX_NEIGHS = 20;
    static const int MAX_STATES = 1 << MAX_PATCHES; // 2 ^ MAX_PATCHES

    CUDAAllostericPatchySwapInteraction();

    virtual ~CUDAAllostericPatchySwapInteraction();

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
        return AllostericPatchySwapInteraction::_N;
    }
};
extern "C" BaseInteraction *make_CUDAAllostericPatchySwapInteraction() {
    return new CUDAAllostericPatchySwapInteraction();
}

#endif //OXDNA_CUDAALLOSTERICPATCHYSWAPINTERACTION_H
