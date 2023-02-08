//
// Created by josh on 7/8/22.
//

#ifndef OXDNA_ALLOSTERICPATCHYPARTICLE_H
#define OXDNA_ALLOSTERICPATCHYPARTICLE_H

#include "../../../../src/Particles/PatchyParticle.h"
#include "AllostericPatch.h"

#define STATE_TRANSITION_SUBDIV 4


// macro grabbed from https://stackoverflow.com/a/2249738
#define GET_BIT(n,k) (n & ( 1 << k )) >> k




/** the State Transition Map is a probabilistic map of how likely the particle is to transition from any given state
    to a different state. the index of the outer vector is the state, and the inner index is probabiliistic. to
    determine if a state transition happens the program generates a random number from 0 to STATE_TRANSITION_SUBDIV
    and transitions the particle state to the state at that index
 */
typedef std::vector<std::array<unsigned int, STATE_TRANSITION_SUBDIV>> StateTransitionMap;

/// the Activation Update Map is a flattened 2d array where the x and y keys are "before" and "after" states
/// and the values are the indexes of patches that will have their activations flipped when the particle
/// flips from state x to state y
typedef std::set<int>* ActivationUpdateMap;

/**
 * An allosteric patchy particle
 * As elegant as it would be to extend PatchyParticle, AllostericPatchyParticle needs to
 * treat patches as objects rather than LR_vectors so I can't
 * TODO: ^^^ I think?
 */
class AllostericPatchyParticle : public BaseParticle {
protected:
    unsigned int _state;
    int _state_size;
    StateTransitionMap* _transitionMap;
    ActivationUpdateMap* _updateMap;
public:
    std::vector<PatchyBond> bonds;
    //LR_vector *_base_patch_positions;

    // Each patch corresponds with a value in BaseParticle::int_centers
    std::vector<AllostericPatch> patches;

    void _set_base_patches();

    void init_allostery(int state_size, StateTransitionMap *transitionMap, ActivationUpdateMap *updateMap);

    AllostericPatchyParticle(int N_patches, int type);
    AllostericPatchyParticle(const AllostericPatchyParticle &b) : AllostericPatchyParticle(b.n_patches(), b.type)
    {
        this->copy_from(b);
    }

    virtual ~AllostericPatchyParticle();

    void set_positions();

    virtual void copy_from(const BaseParticle &);

    AllostericPatchyParticle& operator = (const AllostericPatchyParticle& b) {this->copy_from(b);  return *this;}
    void add_patch(AllostericPatch &patch,int position);

    int get_patch_color(int patchid)  {return this->patches[patchid].get_color();}
    virtual bool is_rigid_body() {
        return true;
    }
//
//    void _set_vertexes(void);
//    void _set_icosahedron_vertexes(void);

//    bool locked_to_particle_id(int particle_id); // {return is_locked() && particle_id == locked_to_particle;}
//    void unlock_patches(void);
//    void set_lock(int patch_idx, int particle=-1,int patch=-1,number energy=0, bool ignore_refresh=false);
//    void unlock(int patch_idx, bool ignore_refresh=false);


    void update_active_patches(int oldState);
    void set_patch_bound(int pidx);

    unsigned int get_state() const;
    void set_state(unsigned int new_state);
    int state_size() const;
    int n_patches() const;
    int n_states() const;

    int do_state_transition(int rand_idx);

};


#endif //OXDNA_ALLOSTERICPATCHYPARTICLE_H
