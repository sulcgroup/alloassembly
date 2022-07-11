/*
 * BaseInteraction is a class that is meant "contain" all the interaction
 * types
 */
/**
 * This file was cloned from another version of Polycubes on Jul. 8th 2022 by josh to provide backwards
 * compatibility with his existing Allosteric Patchy Particles code
 */

#ifndef LEGACY_BASE_INTERACTION_H
#define LEGACY_BASE_INTERACTION_H

#include <map>
#include <cfloat>
#include <fstream>
#include <set>
#include <vector>

#include "defs.h"
#include "Particles/BaseParticle.h"
#include "Boxes/BaseBox.h"
#include "Lists/BaseList.h"
#include "Utilities/Utils.h"
#include "Utilities/oxDNAException.h"
#include "Interactions/Mesh.h"

#include "Lists/Cells.h"

using namespace std;

namespace Legacy {
/**
 * @brief Abstract class defining the interaction interface.
 *
 * We implement the Curiously Recurring Template Pattern (see
 * http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) and therefore
 * we have to keep this class and BaseInteraction separated.
 */
    class IBaseInteraction {
    protected:
        BaseBox *_box;

        /// This is useful for "hard" potentials
        bool _is_infinite;

        /// This are needed to create initial configurations, and they are used by the generator functions
        number _energy_threshold, _temperature;
        /// If true, the generation of the initial configuration will try to take into account the fact that particles that are bonded neighbours should be close to each other
        bool _generate_consider_bonded_interactions;
        /// This controls the maximum at which bonded neighbours should be randomly placed to speed-up generation. Used by generator functions.
        number _generate_bonded_cutoff;

        number _rcut, _sqr_rcut;

        char _topology_filename[256];
    public:
        IBaseInteraction();

        virtual ~IBaseInteraction();

        virtual void set_box(BaseBox *box) { _box = box; }

        virtual void get_settings(input_file &inp);

        /**
         * Initialization of class constants.
         */
        virtual void init() = 0;

        /**
         * @brief Handles particle allocation. Child classes must implement it.
         *
         * @param particles
         * @param N
         */
        virtual void allocate_particles(BaseParticle **particles, int N) = 0;

        /**
         * @brief Returns the maximum cut-off associated to the interaction
         *
         * @return Interaction cutoff
         */
        virtual number get_rcut() { return _rcut; }

        /**
         * @brief Check whether the initial configuration makes sense.
         *
         * Since this operation is interaction-dependent, each interaction must provide its own implementation.
         * @param particles
         * @param N
         */
        virtual void check_input_sanity(BaseParticle **particles, int N) = 0;

        /**
         * @brief Computes the total interaction between particles p and q.
         *
         * If r is not given or NULL, it is computed from scratch. It can optionally update forces and torques exerted on p and q.
         * @param p
         * @param q
         * @param r
         * @param update_forces
         * @return pair-interaction energy
         */
        virtual number
        pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) = 0;

        /**
         * @brief Computes the bonded part interaction between particles p and q.
         *
         * See {\@link pair_interaction} for a description of the parameters.
         */
        virtual number
        pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) = 0;

        /**
         * @brief Computed the non bonded part of the interaction between particles p and q.
         *
         * See {\@link pair_interaction} for a description of the parameters.
         */
        virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL,
                                                  bool update_forces = false) = 0;

        /**
         * @brief Computes the requested term of the interaction energy between p and q.
         *
         * @param name identifier of the interaction method
         * @param p
         * @param q
         * @param r
         * @param update_forces
         * @return
         */
        virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r = NULL,
                                             bool update_forces = false) = 0;

        /**
         * @brief Returns the total potential energy of the system, given a box size
         *
         * @param particles
         * @param N
         * @return
         */
        virtual number get_system_energy(BaseParticle **particles, int N, BaseList *lists);

        /**
         * @brief Like get_system_energy_term, just that the box size can be specified in this case
         *
         * @param name
         * @param particles
         * @param N
         * @return the required energy contribution
         */
        virtual number get_system_energy_term(int name, BaseParticle **particles, int N, BaseList *lists);

        /**
         * @brief Read the topology of the interactions. The defaut
         * method is empty.
         *
         * This function is set up so that the topology of the interactions
         * are read in the code. The different interactions might require
         * different formats, so it is implemented here. The default
         *
         * @param N number of particles as read previously. Used for
         * checking consistency.
         * @param N_strands Pointer to an int which will be filled with the
         * number of "strands" (unbreakable clusters of connected particles)
         * @param particles array of particles.
         */
        virtual void read_topology(int N, int *N_strands, BaseParticle **particles);

        /**
         * @brief Returns the number of particles, as written in the topology.
         *
         * @return number of particles
         */
        virtual int get_N_from_topology();

        /**
         * @brief Like get_system_energy_split, just that the box size can be
         * specified in this case.
         *
         * @return map of the energy contributions
         */
        virtual map<int, number> get_system_energy_split(BaseParticle **particles, int N, BaseList *lists) = 0;

        /**
         * @brief Returns the state of the interaction
         */
        bool get_is_infinite() { return _is_infinite; }

        /**
         * @brief Sets _is_infinite
         *
         * @param arg bool
         */
        void set_is_infinite(bool arg) { _is_infinite = arg; }

        /**
         * @brief overlap criterion. Used only in the generation of configurations. Can be overloaded.
         *
         * The default implementation does not allow two particles to have an
         * interaction strength greater than 100. More subtle implementations
         * may be needed in special cases.
         *
         * @param p
         * @param q
         * @return bool whether the two particles overlap
         */
        virtual bool generate_random_configuration_overlap(BaseParticle *p, BaseParticle *q);


        /**
         * @brief generation of initial configurations. Can be overloaded.
         *
         * The default function creates a random configuration in the most simple way possible:
         * puts particles in one at a time, and using {\@link generate_random_configuration_overlap}
         * to test if two particles overlap.
         *
         * @param particles
         * @param N
         */
        virtual void generate_random_configuration(BaseParticle **particles, int N);
    };


    IBaseInteraction::IBaseInteraction() {
        _energy_threshold = (number) 100.f;
        _is_infinite = false;
        _box = NULL;
        _generate_consider_bonded_interactions = false;
        _generate_bonded_cutoff = 2.0;
    }


    IBaseInteraction::~IBaseInteraction() {

    }


    void IBaseInteraction::get_settings(input_file &inp) {
        getInputString(&inp, "topology", this->_topology_filename, 1);
        getInputNumber(&inp, "energy_threshold", &_energy_threshold, 0);
        getInputNumber(&inp, "T", &_temperature, 1);
        getInputBool(&inp, "generate_consider_bonded_interactions", &_generate_consider_bonded_interactions, 0);
        if (_generate_consider_bonded_interactions) {
            getInputNumber(&inp, "generate_bonded_cutoff", &_generate_bonded_cutoff, 0);
            OX_LOG(Logger::LOG_INFO,
                   "The generator will try to take into account bonded interactions by choosing distances between bonded neighbours no larger than %lf",
                   _generate_bonded_cutoff);
        }
    }


    int IBaseInteraction::get_N_from_topology() {
        char line[512];
        std::ifstream topology;
        topology.open(this->_topology_filename, ios::in);
        if (!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
        topology.getline(line, 512);
        topology.close();
        int ret;
        sscanf(line, "%d %*d\n", &ret);
        return ret;
    }


    void IBaseInteraction::read_topology(int N, int *N_strands, BaseParticle **particles) {
        *N_strands = N;
        allocate_particles(particles, N);
        for (int i = 0; i < N; i++) {
            particles[i]->index = i;
            particles[i]->type = 0;
            particles[i]->strand_id = i;
        }
    }


    number IBaseInteraction::get_system_energy(BaseParticle **particles, int N, BaseList *lists) {
        double energy = 0.;

        vector<ParticlePair> pairs = lists->get_potential_interactions();
        typename std::vector<ParticlePair>::iterator it;
        for (it = pairs.begin(); it != pairs.end(); it++) {
            BaseParticle *p = (*it).first;
            BaseParticle *q = (*it).second;
            energy += (double) pair_interaction(p, q);
            if (this->get_is_infinite()) return energy;
        }

        return (number) energy;
    }


    number IBaseInteraction::get_system_energy_term(int name, BaseParticle **particles, int N, BaseList *lists) {
        number energy = (number) 0.f;

        for (int i = 0; i < N; i++) {
            BaseParticle *p = particles[i];
            vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

            for (unsigned int n = 0; n < neighs.size(); n++) {
                BaseParticle *q = neighs[n];
                if (p->index > q->index) energy += pair_interaction_term(name, p, q);
                if (this->get_is_infinite()) return energy;
            }
        }

        return energy;
    }


    bool IBaseInteraction::generate_random_configuration_overlap(BaseParticle *p, BaseParticle *q) {
        LR_vector dr = _box->min_image(p, q);

        if (dr.norm() >= this->_sqr_rcut) return false;

        // number energy = pair_interaction(p, q, &dr, false);
        number energy = pair_interaction_nonbonded(p, q, &dr, false);

        // in case we create an overlap, we reset the interaction state
        this->set_is_infinite(false);

        // if energy too large, reject
        if (energy > _energy_threshold) return true;
        else return false;
    }


    void IBaseInteraction::generate_random_configuration(BaseParticle **particles, int N) {
        std::vector<BaseParticle *> particles_vector(N);
        std::copy(particles, particles + N, particles_vector.begin());
        Cells c(particles_vector, _box);
        c.init(_rcut);

        for (int i = 0; i < N; i++) {
            BaseParticle *p = particles[i];
            p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y,
                               drand48() * _box->box_sides().z);
        }

        c.global_update();

        for (int i = 0; i < N; i++) {
            BaseParticle *p = particles[i];
            bool same_strand = (_generate_consider_bonded_interactions && i > 0 && p->is_bonded(particles[i - 1]));

            bool inserted = false;
            do {
                if (same_strand)
                    p->pos = particles[i - 1]->pos +
                             LR_vector((drand48() - 0.5), (drand48() - 0.5), (drand48() - 0.5)) *
                             _generate_bonded_cutoff;
                else
                    p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y,
                                       drand48() * _box->box_sides().z);
                // random orientation
                //p->orientation = Utils::get_random_rotation_matrix (2.*M_PI);
                p->orientation = Utils::get_random_rotation_matrix_from_angle(acos(2. * (drand48() - 0.5)));
                p->orientation.orthonormalize();
                p->orientationT = p->orientation.get_transpose();

                p->set_positions();
                p->set_ext_potential(0, this->_box);
                c.single_update(p);

                inserted = true;

                // we take into account the bonded neighbours
                for (unsigned int n = 0; n < p->affected.size(); n++) {
                    BaseParticle *p1 = p->affected[n].first;
                    BaseParticle *p2 = p->affected[n].second;
                    if (p1->index <= p->index && p2->index <= p->index) {
                        number e = pair_interaction_bonded(p1, p2);
                        if (std::isnan(e) || e > _energy_threshold) inserted = false;
                    }
                }

                // here we take into account the non-bonded interactions
                vector<BaseParticle *> neighs = c.get_complete_neigh_list(p);
                for (unsigned int n = 0; n < neighs.size(); n++) {
                    BaseParticle *q = neighs[n];
                    // particles with an index larger than p->index have not been inserted yet
                    if (q->index < p->index && generate_random_configuration_overlap(p, q)) inserted = false;
                }

                // we take into account the external potential
                number boltzmann_factor = exp(-p->ext_potential / _temperature);
                if (std::isnan(p->ext_potential) || drand48() > boltzmann_factor) {
                    inserted = false;
                }

            } while (!inserted);

            if (i > 0 && N > 10 && i % (N / 10) == 0)
                OX_LOG(Logger::LOG_INFO, "Inserted %d%% of the particles (%d/%d)", i * 100 / N, i, N);
        }
    }

/**
 * @brief Base class for managing particle-particle interactions. It is an abstract class.
 *
 * This is the abstract class from which the other interaction classes inherit.
 * The children classes must implement all the public pair_interaction_*
 * methods.  We make use of the Curiously Recurring Template Pattern (see
 * http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) to make
 * it possible to call a method of any derived class from within this class
 * without knowing any detail of that derived class. Hence, we use the {\@link
 * _int_map} map to store function pointers of type child::*.
 * This class lacks a .cpp because we do not want to have template
 * instantiation of the children classes from within here.
 */
    class BaseInteraction : public IBaseInteraction {
    private:

    protected:
        using energy_function = std::function<number(BaseParticle *, BaseParticle *, LR_vector*, bool)>;
        using interaction_map = std::map<int, energy_function*>;
        /**
         * Map of function pointers of individual terms of the particle-particle interaction. Note that
         * these pointers are of type child::* (i.e. methods defined within the child class).
         */
        interaction_map _int_map;

        /**
         * @brief Calls the right interaction-computing method, chosen according to its name.
         *
         * @param that instance of the child class that actually calls the right method
         * @param name identifier of the interaction method
         * @param p
         * @param q
         * @param r
         * @param update_forces
         * @return pair-interaction energy
         */
        virtual number
        _pair_interaction_term_wrapper(BaseInteraction *that, int name, BaseParticle *p, BaseParticle *q, LR_vector *r,
                                       bool update_forces);

        /**
         * @brief Build a mesh by using a function and its derivative.
         *
         * Only derived classes can call this function.
         * @param that pointer to the class owning f and der
         * @param f function
         * @param der derivative of f
         * @param pars pointer to a structure which contains function parameters
         * @param npoints size of the mesh
         * @param xlow the mesh is defined on a finite interval. This is the lower end
         * @param xupp Upper end of the mesh interval
         * @param m mesh to be built
         */
        virtual void _build_mesh(BaseInteraction *that, number (BaseInteraction::*f)(number, void *),
                                 number (BaseInteraction::*der)(number, void *), void *pars, int npoints, number xlow,
                                 number xupp, Mesh &m);

        inline number _query_mesh(number x, Mesh &m);

        inline number _query_meshD(number x, Mesh &m);

    public:
        /**
         * @brief Basic constructor. By default, it does not need anything.
         */
        BaseInteraction();

        virtual ~BaseInteraction();

        /**
         * @brief Computes all the contributions to the total potential energy and returns them in a map.
         *
         * @param particles
         * @param N
         * @return a map storing all the potential energy contributions.
         */
        virtual map<int, number> get_system_energy_split(BaseParticle **particles, int N, BaseList *lists);

    };

    BaseInteraction::BaseInteraction() : IBaseInteraction() {
        this->_rcut = 2.5;
        this->_sqr_rcut = SQR(this->_rcut);
    }

    BaseInteraction::~BaseInteraction() {

    }

    number
    BaseInteraction::_pair_interaction_term_wrapper(BaseInteraction *that, int name, BaseParticle *p, BaseParticle *q,
                                                    LR_vector *r, bool update_forces) {
        LR_vector computed_r;
        if (r == NULL) {
            computed_r = this->_box->min_image(p->pos, q->pos);
            r = &computed_r;
        }

        // _int_map.find(name) returns an iterator which points to the
        // requested function pointer. We dereference it, apply it to the
        // that, which is a pointer to an instance of a class which inherits
        // from BaseInteraction. We then pass in the required parameters and
        // return the interaction energy.
        typename interaction_map::iterator interaction = _int_map.find(name);
        if (interaction == _int_map.end()) {
            throw oxDNAException("%s, line %d: Interaction term '%d' not found", __FILE__, __LINE__, name);
        }
        energy_function func = interaction->second;
        // TODO: make sure this refactor doesn't break EVERYTHING
        return func(p, q, r, update_forces);
    }

    inline number BaseInteraction::_query_meshD(number x, Mesh &m) {
        if (x < m.x_low()) return m.B()[0];
        if (x >= m.xupp()) x = m.xupp() - FLT_EPSILON;
        int i = (int) ((x - m.x_low()) / m.delta());
        number dx = x - m.x_low() - m.delta() * i;
        return (m.B()[i] + (2 * dx * m.C()[i] + 3 * dx * dx * m.D()[i]) * m.inv_sqr_delta());
    }

    inline number BaseInteraction::_query_mesh(number x, Mesh &m) {
        if (x <= m.x_low()) return m.A()[0];
        if (x >= m.xupp()) x = m.xupp() - FLT_EPSILON;
        int i = (int) ((x - m.x_low()) / m.delta());
        number dx = x - m.x_low() - m.delta() * i;
        return (m.A()[i] + dx * (m.B()[i] + dx * (m.C()[i] + dx * m.D()[i]) * m.inv_sqr_delta()));
    }

    void BaseInteraction::_build_mesh(BaseInteraction *that, number (BaseInteraction::*f)(number, void *),
                                      number (BaseInteraction::*der)(number, void *), void *args, int npoints,
                                      number xlow, number xupp, Mesh &m) {
        assert (xlow < xupp);
        int i;
        number x;

        m.init(npoints);

        number dx = (xupp - xlow) / (number) npoints;

        number fx0, fx1, derx0, derx1;

        std::vector<number> A, B, C, D;

        for (i = 0; i < npoints + 1; i++) {
            x = xlow + i * dx;

            fx0 = (that->*f)(x, args);
            fx1 = (that->*f)(x + dx, args);
            derx0 = (that->*der)(x, args);
            derx1 = (that->*der)(x + dx, args);

            A[i] = fx0;
            B[i] = derx0;
            D[i] = (2 * (fx0 - fx1) + (derx0 + derx1) * dx) / dx;
            C[i] = (fx1 - fx0 + (-derx0 - D[i]) * dx);
        }
        m = Mesh(xlow, xupp, dx, 1 / SQR(dx), A, B, C, D);
    }

    map<int, number> BaseInteraction::get_system_energy_split(BaseParticle **particles, int N, BaseList *lists) {
        std::map<int, number> energy_map;

        for (typename interaction_map::iterator it = _int_map.begin(); it != _int_map.end(); it++) {
            int name = it->first;
            energy_map[name] = (number) 0.f;
        }

        for (int i = 0; i < N; i++) {
            BaseParticle *p = particles[i];
            vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

            for (unsigned int n = 0; n < neighs.size(); n++) {
                BaseParticle *q = neighs[n];
                if (p->index > q->index) {
                    for (typename interaction_map::iterator it = _int_map.begin(); it != _int_map.end(); it++) {
                        int name = it->first;
                        energy_map[name] += this->pair_interaction_term(name, p, q);
                    }
                }
            }
        }

        return energy_map;
    }
}
#endif /* LEGACY_BASE_INTERACTION_H */
