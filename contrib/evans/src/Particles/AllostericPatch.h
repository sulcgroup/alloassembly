//
// Created by josh on 7/8/22.
//

#ifndef OXDNA_ALLOSTERICPATCH_H
#define OXDNA_ALLOSTERICPATCH_H
#include "../../../../src/Particles/BaseParticle.h"


/// A structure describing the patch; Each particle can have multiple patches, positioned at different places; The patches are directional and each
/// patch interacts only with its specific complementary patch;
// TODO: have this extend Patch if feasable
// TODO: unclear if this even requires to inherit BaseParticle
class AllostericPatch {
protected:
    LR_vector _position; // the position of the patch with respect to the CM of the particle. initial.
    LR_vector _a1;  //vector that is to be compared against the vector connecting the patches r_pp, needs to be parallel
    LR_vector _a2; // vector which needs to be parallel with a2 on the complementary patch
    LR_vector _init_a1; // initial position of _a1
    LR_vector _init_a2; // initial position of _a2
    int id; //the id of the patch; it is used during initialization to assign patches to particles according to input file; sets the type of patch
    bool active; //is the patch on or not
//    int locked_to_particle; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
//    int locked_to_patch; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
//    number locked_energy;

    int _color; //this is the color of the patch
//    number strength;  //sets the strength of the interaction

    // this is a string representing a conditional that can activate/deactivate the patch
    // default value true means that the patch is always active
    // important to note that this var will be processed to create the particle's overall
    // control table; it is not directly in the simulation
    std::string allostery_conditional;

    bool activation_reversible; //NOT YET IMPLEMENTED

    bool _flipped;
public:
    int index ; //this is the unique index of the patch in the simulation
    bool bound;

    AllostericPatch();

    AllostericPatch(LR_vector a1_xyz, LR_vector a2_xyz, LR_vector position, int id, int color,
                    bool active, std::string allostery_conditional, bool _activation_reversible);

//    bool is_locked() const;

    int get_color() const;
    int get_id() const;
    std::string get_allosteric_conditional() const;

//    bool locked_to(int particle_id,int patch_id) const;
//    bool locked_to_particle_id(int particle_id) const;
//    void get_lock(int& particle_id, int& patch_id) const;
//    number get_lock_energy() const;

    LR_vector a1();
    LR_vector a2();

    void set_a1(LR_vector v);
    void set_a2(LR_vector v);

    LR_vector a1_initial();
    LR_vector a2_initial();

    LR_vector position() const;

    bool is_active() const;
    void set_active(bool bNewVal);
    bool toggle_active();
    bool flipped() const;
    int color() const;
};

#endif //OXDNA_ALLOSTERICPATCH_H
