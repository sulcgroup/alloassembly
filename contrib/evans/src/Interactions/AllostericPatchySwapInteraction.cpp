//
// Created by josh on 7/8/22.
// Copied from Lorenzo Rovigatti's DetailedPatchySwapInteraction class
//

#include "AllostericPatchyInteraction.h"
#include "../../src/Utilities/parse_input/parse_input.h"

using namespace std;


LR_vector getVector(input_file *obs_input,const char *key)
{
    double tmpf[3];
    int tmpi;
    LR_vector v;
    string vec;

    getInputString(obs_input,key,vec,1);
    tmpi = sscanf(vec.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
    if(tmpi != 3)
        throw oxDNAException ("Could not parse vector %s in the input file. Aborting", vec.c_str());
    v.x = tmpf[0];
    v.y = tmpf[1];
    v.z = tmpf[2];
    return v;
}


AllostericPatchySwapInteraction::AllostericPatchySwapInteraction() :
        BaseInteraction() {
    ADD_INTERACTION_TO_MAP(SPHERICAL, _spherical_patchy_two_body);

}

AllostericPatchySwapInteraction::~AllostericPatchySwapInteraction() {

}

void AllostericPatchySwapInteraction::get_settings(input_file &inp) {
    BaseInteraction::get_settings(inp);

    int i; // allocate placeholder int & reuse
    // read number of patches
    getInputInt(&inp, "patch_types_N", &i, 1);
    _N_patch_types = i;
    _base_patch_types.resize(_N_patch_types);

    // read number of particles
    getInputInt(&inp, "particle_types_N", &i, 1);
    _N_per_species.resize(i, 0); // resize particle type count vector
    _base_patch_positions.resize(i);

    getInputNumber(&inp, "DPS_lambda", &_lambda, 0);
    getInputString(&inp, "DPS_interaction_matrix_file", _interaction_matrix_file, 0);

    getInputBool(&inp, "DPS_is_KF", &_is_KF, 0);
    if(_is_KF) {
        getInputInt(&inp, "DPS_patch_power", &_patch_power, 0);
        getInputNumber(&inp, "DPS_KF_delta", &_patch_delta, 1);
        getInputNumber(&inp, "DPS_KF_cosmax", &_patch_cosmax, 1);
    }
    else {
        getInputNumber(&inp, "DPS_sigma_ss", &_sigma_ss, 0);
    }

    getInputNumber(&inp, "DPS_spherical_attraction_strength", &_spherical_attraction_strength, 0.);
    if(_spherical_attraction_strength > 0.) {
        getInputNumber(&inp, "DPS_spherical_rcut", &_spherical_rcut, 1.);
    }
    torsion = false;
    getInputBool(&inp, "torsion", &torsion, 0);

    std::string patchy_file; //this file contains information about types of patches
    getInputString(&inp, "patchy_file", patchy_file, 1);
    std::string particle_file; //this file contains information about types of particles
    getInputString(&inp, "particle_file", particle_file, 1);
    _load_patchy_particle_files(patchy_file, particle_file);
}

void AllostericPatchySwapInteraction::init() {
    _rep_rcut = pow(2., 1. / 6.);
    _sqr_rep_rcut = SQR(_rep_rcut);

    if(_is_KF) {
        ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body_KF);

        _patch_rcut = 1.5 * _patch_delta;

        // the patch-patch radial attraction is centred around 0 (considering that we use the distance between the two surfaces as a variable)
        _sigma_ss = 0.;
        _patch_pow_delta = pow(_patch_delta, (number) 10.);
        _patch_pow_cosmax = pow(1. - _patch_cosmax, (number) _patch_power);
        // this makes sure that at the cutoff the angular modulation is 10^-2
        _patch_angular_cutoff = (1. - _patch_cosmax) * std::pow(4 * std::log(10), 1. / _patch_power);

        OX_LOG(Logger::LOG_INFO, "FS-KF parameters: lambda = %lf, patch_delta = %lf, patch_power = %d, patch_cosmax = %lf, patch_angular_cutoff = %lf", _lambda, _patch_delta, _patch_power, _patch_cosmax, _patch_angular_cutoff);
    }
    else {
        ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body_point);

        _rcut_ss = 1.5 * _sigma_ss;
        _patch_rcut = _rcut_ss;

        number B_ss = 1. / (1. + 4. * SQR(1. - _rcut_ss / _sigma_ss));
        _A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _rcut_ss / _sigma_ss));
        _B_part = B_ss * pow(_sigma_ss, 4.);

        OX_LOG(Logger::LOG_INFO, "FS parameters: lambda = %lf, A_part = %lf, B_part = %lf", _lambda, _A_part, _B_part);
    }

    _sqr_patch_rcut = SQR(_patch_rcut);
    _rcut = 1. + _patch_rcut;

    if(_spherical_attraction_strength > 0.) {
        _sqr_spherical_rcut = SQR(_spherical_rcut);
        _spherical_E_cut = 4. * _spherical_attraction_strength * (1. / pow(_sqr_spherical_rcut, 6) - 1. / pow(_sqr_spherical_rcut, 3));

        if(_spherical_rcut > _rcut) {
            _rcut = _spherical_rcut;
        }
    }

    _sqr_rcut = SQR(_rcut);
}

number AllostericPatchySwapInteraction::_spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    number sqr_r = _computed_r.norm();
    if(sqr_r > _sqr_rcut) {
        return (number) 0.f;
    }

    number energy = (number) 0.f;

    // centre-centre
    if(sqr_r < _sqr_rep_rcut) {
        number ir2 = 1. / sqr_r;
        number lj_part = ir2 * ir2 * ir2;
        energy = 4 * (SQR(lj_part) - lj_part) + 1.0 - _spherical_attraction_strength - _spherical_E_cut;
        if(update_forces) {
            LR_vector force = _computed_r * (-24. * (lj_part - 2 * SQR(lj_part)) / sqr_r);
            p->force -= force;
            q->force += force;
        }
    }
    else {
        if(sqr_r < _sqr_spherical_rcut && _spherical_attraction_strength > 0.) {
            number ir2 = 1. / sqr_r;
            number lj_part = ir2 * ir2 * ir2;
            energy = 4 * _spherical_attraction_strength * (SQR(lj_part) - lj_part) - _spherical_E_cut;
            if(update_forces) {
                LR_vector force = _computed_r * (-24. * _spherical_attraction_strength * (lj_part - 2 * SQR(lj_part)) / sqr_r);
                p->force -= force;
                q->force += force;
            }
        }
    }

    return energy;
}

number AllostericPatchySwapInteraction::_patchy_two_body_point(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    number sqr_r = _computed_r.norm();
    if(sqr_r > _sqr_rcut) {
        return (number) 0.f;
    }

    AllostericPatchyParticle* pp = dynamic_cast<AllostericPatchyParticle*>(p);
    AllostericPatchyParticle* qq = dynamic_cast<AllostericPatchyParticle*>(q);

    number energy = (number) 0.f;
    // loop patches on particle p
    for(uint p_patch_idx = 0; p_patch_idx < p->N_int_centers(); p_patch_idx++) {
        LR_vector p_patch_pos = p->int_centers[p_patch_idx];
        // loop patches on particle q
        for(uint q_patch_idx = 0; q_patch_idx < q->N_int_centers(); q_patch_idx++) {
            AllostericPatch* p_patch = &pp->patches[p_patch_idx];
            AllostericPatch* q_patch = &qq->patches[q_patch_idx];
            if ((p_patch->bound || p_patch->is_active()) && (q_patch->bound || q_patch->is_active())) {
                LR_vector q_patch_pos = q->int_centers[q_patch_idx];

                LR_vector patch_dist = _computed_r + q_patch_pos - p_patch_pos;
                number dist = patch_dist.norm();
                if (dist < _sqr_patch_rcut) {
                    uint p_patch_type = abs(pp->get_patch_color(p_patch_idx));
                    uint q_patch_type = abs(qq->get_patch_color(q_patch_idx));
                    //                uint p_patch_type = _base_patch_positions[p->type][p_patch_idx];
                    //                uint q_patch_type = _base_patch_positions[q->type][q_patch_idx];
                    number epsilon = _patchy_eps[p_patch_type + _N_patch_types * q_patch_type];

                    if (epsilon != 0.) {
                        number r_p = sqrt(dist);
                        number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
                        number tmp_energy = epsilon * _A_part * exp_part * (_B_part / SQR(dist) - 1.);

                        energy += tmp_energy;

                        number tb_energy = (r_p < _sigma_ss) ? epsilon : -tmp_energy;

                        PatchyBond p_bond(q, r_p, p_patch_idx, q_patch_idx, tb_energy);
                        PatchyBond q_bond(p, r_p, q_patch_idx, p_patch_idx, tb_energy);

                        if (pp->patches[p_patch_idx].bound != qq->patches[q_patch_idx].bound) {
                            OX_LOG(Logger::LOG_WARNING,
                                   "Mismatch between binding states between patch %d on particle %d and patch %d on particle %d.",
                                   p->index, p_patch_idx, q->index, q_patch_idx);
                        }
                        if (!pp->patches[p_patch_idx].bound) {
                            pp->update_active_patches(p_patch_idx);
                            qq->update_active_patches(q_patch_idx);
                        }

                        if (update_forces) {
                            number force_mod = epsilon * _A_part * exp_part * (4. * _B_part / (SQR(dist) * r_p)) +
                                               _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
                            LR_vector tmp_force = patch_dist * (force_mod / r_p);

                            LR_vector p_torque = p->orientationT * p_patch_pos.cross(tmp_force);
                            LR_vector q_torque = q->orientationT * q_patch_pos.cross(tmp_force);

                            p->force -= tmp_force;
                            q->force += tmp_force;

                            p->torque -= p_torque;
                            q->torque += q_torque;

                            if (r_p > _sigma_ss) {
                                p_bond.force = tmp_force;
                                p_bond.p_torque = p_torque;
                                p_bond.q_torque = q_torque;

                                q_bond.force = -tmp_force;
                                q_bond.p_torque = -q_torque;
                                q_bond.q_torque = -p_torque;
                            }
                        }

                        _particle_bonds(p).emplace_back(p_bond);
                        _particle_bonds(q).emplace_back(q_bond);

                        if (!no_three_body) {
                            energy += _three_body(p, p_bond, update_forces);
                            energy += _three_body(q, q_bond, update_forces);

                        }
                    }
                }
            }
            else if (p_patch->bound){
                assert(q_patch->bound);
                p_patch->bound = q_patch->bound = false;
                pp->update_active_patches(p_patch_idx);
                qq->update_active_patches(q_patch_idx);
            }
        }
    }

    return energy;
}

// here we compute x^n as (x*x)^((n-1)/2) * x since we now that n is always an odd number
inline double _lr_pow(double x, size_t n){
    double res = x;
    x *= x;

    n = (n - 1) / 2;
    while(n-- > 0){
        res *= x;
    }

    return res;
}

number AllostericPatchySwapInteraction::_patchy_two_body_KF(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    number sqr_r = _computed_r.norm();
    if(sqr_r > _sqr_rcut) {
        return (number) 0.f;
    }

    AllostericPatchyParticle* pp = dynamic_cast<AllostericPatchyParticle*>(p);
    AllostericPatchyParticle* qq = dynamic_cast<AllostericPatchyParticle*>(q);

    number rmod = sqrt(sqr_r);
    LR_vector r_versor = _computed_r / rmod;

    number dist_surf = rmod - 1.;
    number dist_surf_sqr = SQR(dist_surf);
    number r8b10 = SQR(SQR(dist_surf_sqr)) / _patch_pow_delta;
    number exp_part = -1.001 * exp(-(number) 0.5 * r8b10 * dist_surf_sqr);

    number energy = (number) 0.f;
    for(uint p_patch_idx = 0; p_patch_idx < p->N_int_centers(); p_patch_idx++) {
        LR_vector p_patch_pos = p->int_centers[p_patch_idx] * 2;

        number cospr = p_patch_pos * r_versor;
        number cospr_minus_one = cospr - 1.;
        if(cospr_minus_one < _patch_angular_cutoff) {
            number cospr_base = _lr_pow(cospr_minus_one, _patch_power - 1);
            // we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
            number cospr_part = cospr_base * cospr_minus_one;
            number p_mod = exp(-cospr_part / (2. * _patch_pow_cosmax));

            for(uint q_patch_idx = 0; q_patch_idx < q->N_int_centers(); q_patch_idx++) {
                AllostericPatch* p_patch = &pp->patches[p_patch_idx];
                AllostericPatch* q_patch = &qq->patches[q_patch_idx];
                if ((p_patch->bound || p_patch->is_active()) && (q_patch->bound || q_patch->is_active())) {
                    LR_vector q_patch_pos = q->int_centers[q_patch_idx] * 2;

                    number cosqr = -(q_patch_pos * r_versor);
                    number cosqr_minus_one = cosqr - 1.;
                    if (cosqr_minus_one < _patch_angular_cutoff) {
                        uint p_patch_type = pp->get_patch_color(p_patch_idx);
                        uint q_patch_type = qq->get_patch_color(q_patch_idx);
                        //                    uint p_patch_type = _base_patch_positions[p->type][p_patch_idx];
                        //                    uint q_patch_type = _base_patch_positions[q->type][q_patch_idx];
                        number epsilon = _patchy_eps[p_patch_type + _N_patch_types * q_patch_type];

                        // if the epsilon value (binding strength coefficient) between these two patches is nonzero
                        if (epsilon != 0.) {
                            number cosqr_base = _lr_pow(cosqr_minus_one, _patch_power - 1);
                            number cosqr_part = cosqr_base * cosqr_minus_one;
                            number q_mod = exp(-cosqr_part / (2. * _patch_pow_cosmax));

                            number tmp_energy = exp_part * p_mod * q_mod;

                            if (tmp_energy < 0.) {
                                energy += tmp_energy;

                                // when we do the swapping the radial part is the one that gets to one beyond the minimum, while the angular part doesn't change
                                number tb_energy = (dist_surf < _sigma_ss) ? epsilon * p_mod * q_mod : -tmp_energy;
                                PatchyBond p_bond(q, dist_surf, p_patch_idx, q_patch_idx, tb_energy);
                                PatchyBond q_bond(p, dist_surf, q_patch_idx, p_patch_idx, tb_energy);

                                if (update_forces) {
                                    // radial part
                                    LR_vector radial_force =
                                            r_versor * (p_mod * q_mod * 5. * (rmod - 1.) * exp_part * r8b10);

                                    // angular p part
                                    number der_p = exp_part * q_mod *
                                                   (0.5 * _patch_power * p_mod * cospr_base / _patch_pow_cosmax);
                                    LR_vector p_ortho = p_patch_pos - cospr * r_versor;
                                    LR_vector angular_force = p_ortho * (der_p / rmod);

                                    // angular q part
                                    number der_q = exp_part * p_mod *
                                                   (-0.5 * _patch_power * q_mod * cosqr_base / _patch_pow_cosmax);
                                    LR_vector q_ortho = q_patch_pos + cosqr * r_versor;
                                    angular_force += q_ortho * (der_q / rmod);

                                    LR_vector tot_force = radial_force + angular_force;

                                    LR_vector p_torque = p->orientationT * (r_versor.cross(p_patch_pos) * der_p);
                                    LR_vector q_torque = q->orientationT * (q_patch_pos.cross(r_versor) * der_q);

                                    p->force -= tot_force;
                                    q->force += tot_force;

                                    p->torque -= p_torque;
                                    q->torque += q_torque;

                                    p_bond.force = (dist_surf < _sigma_ss) ? angular_force : tot_force;
                                    p_bond.p_torque = p_torque;
                                    p_bond.q_torque = q_torque;

                                    q_bond.force = (dist_surf < _sigma_ss) ? -angular_force : -tot_force;
                                    q_bond.p_torque = -q_torque;
                                    q_bond.q_torque = -p_torque;
                                }

                                _particle_bonds(p).emplace_back(p_bond);
                                _particle_bonds(q).emplace_back(q_bond);

                                if (!no_three_body) {
                                    energy += _three_body(p, p_bond, update_forces);
                                    energy += _three_body(q, q_bond, update_forces);
                                }
                            }
                        }
                    }
                }
                else if (p_patch->bound){
                    assert(q_patch->bound);
                    p_patch->bound = q_patch->bound = false;
                    pp->update_active_patches(p_patch_idx);
                    qq->update_active_patches(q_patch_idx);
                }
            }
        }
    }

    return energy;
}

number AllostericPatchySwapInteraction::_three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces) {
    number energy = 0.;

    number curr_energy = new_bond.energy;
    const auto &p_bonds = _particle_bonds(p);
    for(auto &other_bond : p_bonds) {
        // three-body interactions happen only when the same patch is involved in more than a bond
        if(other_bond.other != new_bond.other && other_bond.p_patch == new_bond.p_patch) {
            number other_energy = other_bond.energy;

            energy += _lambda * curr_energy * other_energy;

            if(update_forces) {
                {
                    BaseParticle *other = new_bond.other;

                    number factor = -_lambda * other_energy;
                    LR_vector tmp_force = factor * new_bond.force;

                    p->force -= tmp_force;
                    other->force += tmp_force;

                    p->torque -= factor * new_bond.p_torque;
                    other->torque += factor * new_bond.q_torque;
                }

                {
                    BaseParticle *other = other_bond.other;

                    number factor = -_lambda * curr_energy;
                    LR_vector tmp_force = factor * other_bond.force;

                    p->force -= tmp_force;
                    other->force += tmp_force;

                    p->torque -= factor * other_bond.p_torque;
                    other->torque += factor * other_bond.q_torque;
                }
            }
        }
    }

    return energy;
}

void AllostericPatchySwapInteraction::begin_energy_computation() {
    BaseInteraction::begin_energy_computation();

    for(int i = 0; i < _N; i++) {
        _particle_bonds(CONFIG_INFO->particles()[i]).clear();
    }
}

number AllostericPatchySwapInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r) {
        if(q != P_VIRTUAL && p != P_VIRTUAL) {
            _computed_r = _box->min_image(p->pos, q->pos);
        }
    }

    number energy = pair_interaction_bonded(p, q, false, update_forces);
    energy += pair_interaction_nonbonded(p, q, false, update_forces);
    return energy;
}

number AllostericPatchySwapInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    number energy = 0.;

    return energy;
}

number AllostericPatchySwapInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r) {
        _computed_r = _box->min_image(p->pos, q->pos);
    }

    number energy = _spherical_patchy_two_body(p, q, false, update_forces);

    if(_is_KF) {
        energy += _patchy_two_body_KF(p, q, false, update_forces);
    }
    else {
        energy += _patchy_two_body_point(p, q, false, update_forces);
    }

    return energy;
}

void AllostericPatchySwapInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
    int N = particles.size();
    int curr_limit = _N_per_species[0];
    int curr_species = 0;
    for(int i = 0; i < N; i++) {
        // if we've hit the limit for whichever species we are on
        if(i == curr_limit) {
            // advance to the next species
            curr_species++;
            curr_limit += _N_per_species[curr_species];
        }
        particles[i] = new AllostericPatchyParticle(_base_particle_types[curr_species]);

        particles[i]->index = i;
        particles[i]->strand_id = i;
        particles[i]->type = particles[i]->btype = curr_species;
    }
}

void AllostericPatchySwapInteraction::_parse_interaction_matrix() {
    // parse the interaction matrix file
    input_file inter_matrix_file;
    inter_matrix_file.init_from_filename(_interaction_matrix_file);
    _patchy_eps.resize(_N_patch_types * _N_patch_types, 0.);

    if(inter_matrix_file.state == ERROR) {
        OX_LOG(Logger::LOG_INFO, ("No interaction matrix file " + _interaction_matrix_file + " found. Using default bindings.").c_str());
        for (int i = 0; i < _N_patch_types; i++) {
            for (int j = 0; j < _N_patch_types; j++) {
                if (_base_patch_types[i].color() == -_base_patch_types[j].color()) {
                    _patchy_eps[i + _N_patch_types * j] = _patchy_eps[j + _N_patch_types * i] = 1.0;
                }
            }
        }
//        throw oxDNAException("Caught an error while opening the interaction matrix file '%s'", _interaction_matrix_file.c_str());
    }
    else {
        OX_LOG(Logger::LOG_INFO, ("Loading interaction matrix from file " + _interaction_matrix_file + ".").c_str());

        for (int i = 0; i < _N_patch_types; i++) {
            for (int j = 0; j < _N_patch_types; j++) {
                number value;
                std::string key = Utils::sformat("patchy_eps[%d][%d]", i, j);
                if (getInputNumber(&inter_matrix_file, key.c_str(), &value, 0) == KEY_FOUND) {
                    _patchy_eps[i + _N_patch_types * j] = _patchy_eps[j + _N_patch_types * i] = value;
                }
            }
        }
    }
}

void AllostericPatchySwapInteraction::_load_patchy_particle_files(std::string &patchy_file, std::string &particle_file)
{
    //first process patches
    FILE *fpatch = fopen(patchy_file.c_str(),"r");
    if(!fpatch) // if file not found
    {
        throw oxDNAException("Could not open file %s ",patchy_file.c_str());
    }
    input_file obs_input;
    obs_input.init_from_file(fpatch); // open patch file input

    // patch number
    int patch_idx = 0;
    char patch_no[1024]; //instantiate C string buffer - 1024 chars should be enough
    snprintf(patch_no, 1020, "patch_%d", patch_idx); // write patch index
    std::string patch_string;

    while(  getInputString(&obs_input,patch_no,patch_string,0) == KEY_FOUND )
    {
        AllostericPatch patch = _process_patch_type(patch_string);
        if(patch_idx >= _N_patch_types)
        {
            throw oxDNAException("Number of patch types is larger than N_patch_types = %d ",_N_patch_types);
        }
        _base_patch_types[patch_idx] = patch;
        patch_idx++;
        snprintf(patch_no, 1020, "patch_%d", patch_idx);
    }

    fclose(fpatch);

    //now process particles
    FILE *fparticle = fopen(particle_file.c_str(), "r");
    if(!fparticle)
    {
        throw oxDNAException("Could not open file %s ",particle_file.c_str());
    }
    input_file p_input;
    p_input.init_from_file(fparticle);

    int particle_idx = 0;
    char particle_no[1024];
    snprintf(particle_no, 1020, "particle_%d", particle_idx);
    std::string particle_string;

    while (getInputString(&p_input, particle_no, particle_string, 0) == KEY_FOUND) {
        AllostericPatchyParticle particle = _process_particle_type(particle_string);
        _base_particle_types.push_back(particle);
//        if(particle_idx >= _N_particle_types)
//            throw oxDNAException ("More particle types in particle config file than specified in the input file. Aborting");

        _base_patch_positions[particle_idx] = particle.int_centers; // cache patch positions so CUDA has access to them
        particle_idx++;
        snprintf(particle_no, 1020, "particle_%d", particle_idx);
    }

    fclose(fparticle);

    OX_LOG(Logger::LOG_INFO, "Loaded %d patch types and %d particle types", patch_idx, particle_idx);
    if(particle_idx != get_num_base_particle_types() || patch_idx != _N_patch_types)
        throw oxDNAException ("More (or less) particle or patches types in particle/patchy config file than specified in the input file. Aborting");
    _parse_interaction_matrix();
}

AllostericPatch AllostericPatchySwapInteraction::_process_patch_type(std::string input_string)
{
    input_file *obs_input = Utils::get_input_file_from_string(input_string);
    // the id of the patch. used to reference patch from particles.txt
    int id;
    // the color of the patch. used to determine if binding is valid
    int color;
    // the binding strength. defaults to 1
    float strength = 1.0f;
    // patch angles.
    //a1 is the patch orientation vector (should be tangent to the surface of the particle)
    //a2 is the "up" vector, should be perpendicular to the surface of the particle
    LR_vector a1, a2;
    // the position of the patch relative to the center of the particle
    LR_vector position;
    string vec;

    // a logical string representing the conditional statement controlling whether this patch is active
    std::string allostery_conditional;

    // load id, color, strength
    getInputInt(obs_input,"id",&id,1);
    getInputInt(obs_input,"color",&color,1);
    getInputFloat(obs_input,"strength",&strength,1);

    // load patch angles and position vector
    a1 = getVector(obs_input,"a1");
    a2 = getVector(obs_input,"a2");
    position = getVector(obs_input,"position");

    // normalize patch angle vectors, since vector magnitude really should not be relevant here
    a1 = a1 / a1.norm();
    a2 = a2 / a2.norm();

    getInputString(obs_input, "allostery_conditional", allostery_conditional, 1); //1?

    // construct patch. can be a straight up object since this will be immutable
    AllostericPatch loaded_patch(a1, a2, position, id, color, true, allostery_conditional, false);

    //printf("Loaded patch %d with color %d \n",loaded_patch.id,loaded_patch.color);
    // deallocate memory for file reader
    delete obs_input;
    //return patch object
    return loaded_patch;
}

AllostericPatchyParticle AllostericPatchySwapInteraction::_process_particle_type(std::string input_string)
{
    input_file *obs_input = Utils::get_input_file_from_string(input_string);
    int type;
    getInputInt(obs_input,"type",&type,1);

    // init a vector to store patches for this particle
    std::vector<AllostericPatch > particle_patches;
    std::string patches;
    if( getInputString(obs_input,"patches",patches,1) == KEY_FOUND )
    {
        //now process a list of patches: 1,2,3,4
        std::replace( patches.begin(), patches.end(), ',', ' ');
        std::stringstream s(patches);
        int patch_id;
        while( s >> patch_id)
        {
            AllostericPatch patch(this->_base_patch_types[patch_id]);
            particle_patches.push_back(patch);
            OX_LOG(Logger::LOG_INFO,"Particle of type %d adding a patch of color %d",type,patch.color());
            //s >> patch_id;
        }
    }

    int N_vertexes = 0;
    OX_LOG(Logger::LOG_INFO,"Particle of type %d has %d vertexes",type,N_vertexes);
    AllostericPatchyParticle p(particle_patches.size(), type);
    int position = 0;
    for(typename std::vector<AllostericPatch >::iterator i = particle_patches.begin(); i != particle_patches.end(); ++i)
    {
        p.add_patch(*i,position);
        position++;
    }
    p._set_base_patches();
    p.init_allostery();

//	ParticleStateChange test_change(default_state, 2, 0);

//	std::vector<int> t = (*p.allostery_map)[test_change];
    return p;
}

//std::vector<LR_vector> AllostericPatchySwapInteraction::_parse_base_patches(std::string filename, int N_patches) {
//    std::ifstream patch_file(filename);
//    if(!patch_file.good()) {
//        throw oxDNAException("Can't read patch file '%s'. Aborting", filename.c_str());
//    }
//
//    std::vector<LR_vector> base_patches(N_patches);
//    string line;
//    for(int i = 0; i < N_patches; i++) {
//        if(!patch_file.good()) {
//            throw oxDNAException("The patch file '%s' does not seem to contain enough lines (%d found, should be %d)", filename.c_str(), i, N_patches);
//        }
//        std::getline(patch_file, line);
//        auto spl = Utils::split(line);
//        if(spl.size() != 3) {
//            throw oxDNAException("Patch file '%s': invalid line '%s'", filename.c_str(), line.c_str());
//        }
//        LR_vector v;
//        v[0] = std::stof(spl[0]);
//        v[1] = std::stof(spl[1]);
//        v[2] = std::stof(spl[2]);
//
//        base_patches[i] = v;
//    }
//
//    patch_file.close();
//
//    return base_patches;
//}

void AllostericPatchySwapInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
    _N = particles.size();
    *N_strands = _N;
    int N_types;
    std::ifstream topology(this->_topology_filename, ios::in); // open a file stream to topology file
    if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
    char line[4096]; // bad code, has caused segfault issue in the past. TODO: rewrite
    topology.getline(line, 512);
    sscanf(line, "%*d %d\n", &N_types);
    //second line specifies numbero f particles of each  type
    topology.getline(line,4090);
    //printf ("N:%d FIRST LINE:--%s--\n", N, line);

    std::stringstream ss(line); // create a string stream for line

    //int count_type;
    int total_count = 0;
    int type = 0;
    while (ss >> type)
    {
        //printf("Loaded type %d, and state is %d\n",type,ss.good());
        fflush(stdout);
        if(total_count >= _N || type >= _base_particle_types.size() || type < 0)
            throw oxDNAException("The sum of number of species is larger than number of particles, or unknown type encountered. Aborting while processing file %s (%d %d %d)",
                                 this->_topology_filename,
                                 total_count,
                                 _base_particle_types.size(),
                                 type);
        int i = total_count;
        particles[i] = new AllostericPatchyParticle(this->_base_particle_types[type]);
        _N_per_species[type]++;

        /*
        printf("!!!!!!!!! JUST INITIALIZED PARTICLE %d \n",i);
        PatchyShapeParticle *p =  dynamic_cast<PatchyShapeParticle *>( particles[i]);
        printf("N_patch %d N_vertex %d N_int %d\n",p->N_patches,p->N_vertexes,p->N_int_centers);
        for(int ii = 0; ii < 12; ii++)
           {

                  printf("%f %f %f \n",p->_vertexes[ii].x,p->_vertexes[ii].y,p->_vertexes[ii].z);
           }
         */

        particles[i]->index = i;
        particles[i]->type = type;
        particles[i]->btype = type;
        particles[i]->strand_id = i;

        total_count++;
        //ss >> type;
        //printf("at the end of while, Loaded type %d, and state is %d\n",type,ss.good());

    }
    OX_LOG(Logger::LOG_INFO, "There were %d particles, %d types, and finished allocation, and line was %s, and N_particle types was %d",_N,N_types,line,_base_particle_types.size());
    int patch_index = 0;
    for(int i = 0; i < _N; i++)
    {
        //printf("Particle %d has %d patches, which have the following colors: ",i,dynamic_cast<PatchyShapeParticle *>(particles[i])->N_patches);
        particles[i]->set_positions();
        for(int c = 0; c < dynamic_cast<AllostericPatchyParticle *>( particles[i])->n_patches(); c++)
        {
            //now assign each patch its unique id
            //printf("%d ",dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].color);
            dynamic_cast<AllostericPatchyParticle *>(particles[i])->patches[c].index = patch_index;
            patch_index++;
            //printf("%d (%f %f %f) ",dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].color, dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].a1.x,dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].a1.y,dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].a1.z);
        }
        //printf("\n");
    }
    allocate_particles(particles);

}

void AllostericPatchySwapInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

extern "C" AllostericPatchySwapInteraction* make_AllostericPatchySwapInteraction() {
    return new AllostericPatchySwapInteraction();
}
