//
// Created by josh on 7/8/22.
//

#ifndef OXDNA_CUDAALLOSTERICPATCHYINTERACTION_H
#define OXDNA_CUDAALLOSTERICPATCHYINTERACTION_H

#include "AllostericPatchyInteraction.h"
#include "CUDA/Interactions/CUDABaseInteraction.h"

template<typename number>
class CUDAAllostericPatchyInteraction : public CUDABaseInteraction, public AllostericPatchyInteraction<number>{

};

#endif //OXDNA_CUDAALLOSTERICPATCHYINTERACTION_H
