//
// Created by josh on 7/8/22.
// Copied from Lorenzo Rovigatti's DetailedPatchySwapInteraction class
//

#ifndef OXDNA_ALLOSTERICPATCHYSWAPINTERACTION_H
#define OXDNA_ALLOSTERICPATCHYSWAPINTERACTION_H

// macro grabbed from https://stackoverflow.com/a/2249738
#define GET_BIT(n,k) (n & ( 1 << k )) >> k

#include "../../../../src/Interactions/BaseInteraction.h"
#include "../Particles/AllostericPatchyParticle.h"
class AllostericPatchySwapInteraction : public BaseInteraction {
protected:
    /// Number of particles of each species
    std::vector<int> _N_per_species;
    int _N = 0; // number of individual particles
    int _N_patch_types = 0;
    int _N_particle_types = 0;

    /// Patch types
    /// positions of each patch on each particle type, cached here so CUDA has access to them
    std::vector<std::vector<LR_vector>> _base_patch_positions;
    /// List of particle types
    std::vector<AllostericPatchyParticle> _base_particle_types;
    /// List of patch types (objects)
    std::vector<AllostericPatch> _base_patch_types;
//    /// Patch type for each particle species
//    std::vector<std::vector<int>> _base_patch_positions;
//    /// Base position of the patches for each particle species
//    std::vector<std::vector<LR_vector>> _base_patch_positions;

    std::string _interaction_matrix_file;

    /// Repulsive interaction energy at the cut-off
    number _rep_rcut = 0.;
    number _sqr_rep_rcut = -1.;

    /// Patchy-related quantities
    std::vector<number> _patchy_eps;
    number _patch_rcut = -1.;
    number _sqr_patch_rcut = -1.;
    number _sigma_ss = 0.4;
    number _rcut_ss = -1.;
    number _lambda = 1.0;
    number _A_part = 0.;
    number _B_part = 0.;

    /// KF-related quantities
    /// Width of the patches
    bool _is_KF = false;
    /// Exponent for the Gaussian-like potential well used for the patches
    int _patch_power = 30;
    /// delta = patch width. confusingly, Flavio calls this "alpha" (I think?)
    number _patch_delta;
    /// Angular width of the patches
    number _patch_cosmax;
    /// _patch_alpha^10
    number _patch_pow_delta;
    /// _patch_cosmax^30
    number _patch_pow_cosmax;
    /// Angular cut-off for the patchy attraction
    number _patch_angular_cutoff;

    /// Optional spherical attraction
    number _spherical_attraction_strength = 0.;
    number _spherical_rcut = 2.5;
    number _sqr_spherical_rcut = 6.25;
    number _spherical_E_cut = 0.;

    /// TODO: torsion
    bool torsion;

    void _parse_interaction_matrix();
//    std::vector<LR_vector> _parse_base_patches(std::string filename, int N_patches);

    number _patchy_two_body_point(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _patchy_two_body_KF(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

    number _three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces);

    void _load_patchy_particle_files(std::string& patchy_file, std::string& particle_file);
    AllostericPatch _process_patch_type(std::string input_string);
    AllostericPatchyParticle _process_particle_type(std::string input_string);


    inline std::vector<PatchyBond> &_particle_bonds(BaseParticle *p) {
        return static_cast<AllostericPatchyParticle *>(p)->bonds;
    }

public:
    enum {
        PATCHY = 0,
        SPHERICAL = 1
    };

    bool no_three_body = false;

    AllostericPatchySwapInteraction();
    virtual ~AllostericPatchySwapInteraction();

    virtual void get_settings(input_file &inp);
    virtual void init();

    virtual void allocate_particles(std::vector<BaseParticle *> &particles);

    void begin_energy_computation() override;

    virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
    virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
    virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);


    virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
    virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

    int get_num_base_particle_types() const {return _base_particle_types.size();}
};

extern "C" AllostericPatchySwapInteraction *make_AllostericPatchySwapInteraction();


#endif //OXDNA_ALLOSTERICPATCHYSWAPINTERACTION_H
