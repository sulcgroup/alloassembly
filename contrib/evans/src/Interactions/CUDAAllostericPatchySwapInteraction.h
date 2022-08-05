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
    llint _step;
public:
    static const int MAX_PATCHES = 5;
    static const int MAX_SPECIES = 5;
    static const int MAX_NEIGHS = 20;
    static const int MAX_STATES = 1 << MAX_PATCHES; // 2 ^ MAX_PATCHES

    CUDAAllostericPatchySwapInteraction();
    virtual ~CUDAAllostericPatchySwapInteraction();

    void get_settings(input_file &inp);
    void cuda_init(c_number box_side, int N);
    c_number get_cuda_rcut() {
        return this->get_rcut();
    }

    void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAAllostericPatchySwapInteraction() {
    return new CUDAAllostericPatchySwapInteraction();
}

#endif //OXDNA_CUDAALLOSTERICPATCHYSWAPINTERACTION_H