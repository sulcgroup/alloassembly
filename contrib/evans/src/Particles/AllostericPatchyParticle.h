//
// Created by josh on 7/8/22.
//

#ifndef OXDNA_ALLOSTERICPATCHYPARTICLE_H
#define OXDNA_ALLOSTERICPATCHYPARTICLE_H

#include "../../../../src/Particles/PatchyParticle.h"
#include "AllostericPatch.h"
#include "ParticleStateChange.h"

/**
 * An allosteric patchy particle
 * As elegant as it would be to extend PatchyParticle, AllostericPatchyParticle needs to
 * treat patches as objects rather than LR_vectors so I can't
 * TODO: ^^^ I think?
 */
class AllostericPatchyParticle : public BaseParticle {
protected:
    std::unordered_map<ParticleStateChange, std::vector<int>>* allostery_map;
public:
    std::vector<PatchyBond> bonds;
    //LR_vector *_base_patches;

    // Each patch corresponds with a value in BaseParticle::int_centers
    std::vector<AllostericPatch> patches;

    void _set_base_patches();

    void init_allostery();

public:
    AllostericPatchyParticle(int N_patches, int type);
    AllostericPatchyParticle(const AllostericPatchyParticle &b)
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

    virtual bool patch_status(bool* binding_status, int patch_idx) const;
    bool patch_status(bool* particle_status, std::string logic) const;


        virtual void update_active_patches(int toggle_idx);

    bool* get_binding_state() const;
    int n_patches() const;

};


#endif //OXDNA_ALLOSTERICPATCHYPARTICLE_H
