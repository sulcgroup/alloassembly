//
// Created by josh on 7/8/22.
//

#ifndef OXDNA_PARTICLESTATECHANGE_H
#define OXDNA_PARTICLESTATECHANGE_H

#include "../../../../src/Particles/BaseParticle.h"

struct ParticleStateChange {
    const bool* _state; // binding state
    const int _state_size; // size of binding state
    const int _change; // the index of the patch that binds or unbinds
    ParticleStateChange(bool* state, int state_size, int change) : _state(state), _state_size(state_size), _change(change){}
    ~ParticleStateChange(){
        delete [] this->_state;
    }
    bool operator==(const ParticleStateChange &other) const{
        // highly unlikely to occur but if it did the consequences would be... bad
        if (this->_state_size != other._state_size){
            return false;
        }
        for (int i = 0; i < this->_state_size; i++) {
            if (this->_state[i] != other._state[i]){
                return false;
            }
        }
        return true;
    }
};

template <>
struct std::hash<ParticleStateChange>{
    std::size_t operator()(ParticleStateChange const &st) const{
//		std::size_t size_hash;
//		for (int i = 0; i < st._state_size; i++)
//		{
//			size_hash ^= int(st._state[i] << i); // HOPING THIS WILL WORK
//		}
//		return size_hash ^ (st._change << 1); // HOPING THIS WORKS
        std::size_t size_hash = 0;
        int i;
        for (i = 0; i < st._state_size; i++){
            if (st._state[i]){
                size_hash += pow(2, i);
            }
        }
//		size_hash ^= st._change >> st._state_size;
        size_hash += st._change * std::size_t(pow(2, st._state_size));
        return size_hash;
    }
};
#endif //OXDNA_PARTICLESTATECHANGE_H
