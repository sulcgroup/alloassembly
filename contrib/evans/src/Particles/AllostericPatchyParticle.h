/*
 * PatchyShapeParticle.h
 *
 *  Created on: 11/june/2022
 *      Author: josh, based off lorenzo's work
 */

#ifndef ALLOPATCHYPARTICLE_H_
#define ALLOPATCHYPARTICLE_H_

#include "../../../../src/Particles/BaseParticle.h"
#include "../../../romano/src/Particles/PatchyShapeParticle.h"
#include <unordered_map>

struct ParticleStateChange {
	const bool* _state;
	const int _state_size;
	const int _change;
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


/// A structure describing the patch; Each particle can have multiple patches, positioned at different places; The patches are directional and each
/// patch interacts only with its specific complementary patch;
// TODO: have this extend Patch if feasable
template <typename number>
struct AllostericPatch  : public BaseParticle<number> {
 LR_vector<number> position; //the position of the patch with respect to the CM of the particle
 LR_vector<number> a1;  //vector that is to be compared against the vector connecting the patches r_pp, needs to be parallel
 LR_vector<number> a2; // vector which needs to be parallel with a2 on the complementary patch
 int id; //the id of the patch; it is used during initialization to assign patches to particles according to input file; sets the type of patch
 int index ; //this is the unique index of the patch in the simulation
 bool active; //is the patch on or not
 int locked_to_particle; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
 int locked_to_patch; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
 number locked_energy;

 int color; //this is the color of the patch
 number strength;  //sets the strength of the interaction
 number a1_x, a1_y, a1_z;
 number a2_x, a2_y, a2_z;

 // this is a string representing a conditional that can activate/deactivate the patch
 // default value true means that the patch is always active
 // important to note that this var will be processed to create the particle's overall
 // control table; it is not directly in the simulation
 std::string allostery_conditional;

 bool activation_reversible; //NOT YET IMPLEMENTED

 bool flipped;

    AllostericPatch() {active = false; id = 0; color = -1; strength = 1; a1_x = a1_y = a1_z = a2_x = a2_y = a2_z = 0; set_lock(-1,-1,0);}

    AllostericPatch(LR_vector<number> _a1_xyz, LR_vector<number> _a2_xyz, LR_vector<number> _position, int _id,int _color, number _strength=1.0,  bool _active = true, std::string _allostery_conditional = "(true)", bool _activation_reversible = false) :
	 position(_position), id(_id), active(_active), color(_color), strength(_strength), allostery_conditional(_allostery_conditional), activation_reversible(_activation_reversible), flipped(false)
 {
	 a1_x = _a1_xyz.x;
	 a1_y = _a1_xyz.y;
	 a1_z = _a1_xyz.z;
	 a2_x = _a2_xyz.x;
	 a2_y = _a2_xyz.y;
	 a2_z = _a2_xyz.z;

	 set_lock(-1,-1,0);
 }
 /* deprecated. Use PatchyShapeParticle.set_lock(...) instead*/
 void set_lock(int particle=-1,int patch=-1,number energy=0)
 {
	 locked_to_particle = particle;
	 locked_to_patch = patch;
	 locked_energy = energy;
 }

 bool is_locked(void) {return locked_to_particle >= 0;}

 int get_color(void) {return color;}

 bool locked_to(int particle_id,int patch_id) {return is_locked() && (particle_id == locked_to_particle && locked_to_patch == patch_id);}
 bool locked_to_particle_id(int particle_id)  {return is_locked() && particle_id == locked_to_particle;}
 void get_lock(int& particle_id, int& patch_id) {particle_id = locked_to_particle; patch_id =  locked_to_patch;}
 number get_lock_energy(void) {return locked_energy;}

 /** deprecated. use PatchyShapeParticle.unlock() instead*/
 void unlock() {locked_to_particle = -1;}

 bool is_active() const {return this->active;}
 void set_active(bool bNewVal) {this->active = bNewVal;}
 bool toggle_active()
 {
	 if (this->activation_reversible || !this->flipped)
	 {
		 this->active = !this->active;
		 this->flipped = true;
	 }
	 return this->active;
 }
};

/**
 * @brief Incapsulates a patchy particle with 2, 3, or 4 spherical patches. Used by PatchyInteraction.
 */
template<typename number>
class AllostericPatchyParticle : public BaseParticle<number> {
public:
	//LR_vector<number> *_base_patches;
	int N_patches; //number of patches  = number of patches
	int N_vertexes; //number of vertices of the shape; 0 = sphere

    AllostericPatch<number> *patches;
    LR_vector<number> *_vertexes;

	void _set_base_patches();

    void init_allostery();

public:
    AllostericPatchyParticle(int N_patches=1 , int type = 0,int N_vertexes=0);
    AllostericPatchyParticle(const AllostericPatchyParticle<number> &b)
	{
		patches = 0;
		this->_vertexes =  0;
		this->copy_from(b);
	}

	virtual ~AllostericPatchyParticle();

	void set_positions();

	virtual void copy_from(const BaseParticle<number> &);

    AllostericPatchyParticle<number>& operator = (const AllostericPatchyParticle<number>& b) {this->copy_from(b);  return *this;}
	void add_patch(AllostericPatch<number> &patch,int position);

	int get_patch_color(int patchid)  {return this->patches[patchid].get_color();}
	virtual bool is_rigid_body() {
		return true;
	}

	void _set_vertexes(void);
	void _set_icosahedron_vertexes(void);

	bool locked_to_particle_id(int particle_id); // {return is_locked() && particle_id == locked_to_particle;}
	void unlock_patches(void);
	void set_lock(int patch_idx, int particle=-1,int patch=-1,number energy=0, bool ignore_refresh=false);
	void unlock(int patch_idx, bool ignore_refresh=false);

	virtual bool patch_status(bool* binding_status, int patch_idx) const;

	virtual void update_active_patches(int toggle_idx);

	bool* get_binding_state() const;

};

template <typename number>
class SimpleAllosteryPatchyShapeParticle : public AllostericPatchyParticle<number>{
public:
	std::unordered_map<ParticleStateChange, std::vector<int>>* allostery_map;

public:
	SimpleAllosteryPatchyShapeParticle(int N_patches=1 , int type = 0,int N_vertexes=0) : AllostericPatchyParticle<number>(N_patches, type, N_vertexes){};
	SimpleAllosteryPatchyShapeParticle(const SimpleAllosteryPatchyShapeParticle<number> &p){
		this->patches = 0;
		this->_vertexes =  0;
		this->copy_from(p);
	}

	void copy_from(const BaseParticle<number> &bb);

	bool patch_status(bool* binding_status, int patch_idx) const{
	//		return this->patch_status(bindingstatus, this->patches[patch_idx].allostery_conditional);
			bool status = this->patch_status(binding_status, this->patches[patch_idx].allostery_conditional);
			//DEBUG LOGGING - might actually want to keep some of this in
			std::string status_str = "[";
			for (int i = 0; i < this->N_patches; i++) {
				status_str += (binding_status[i] ? " T" : " F");
			}
			status_str += "]";
			if ( this->patches[patch_idx].allostery_conditional != "true") {
				OX_LOG(Logger::LOG_INFO, "Patch %d of particle type %d with binding state %s will be %s according to patch logic \"%s\".", patch_idx, this->type, status_str.c_str(), std::to_string(status).c_str(), this->patches[patch_idx].allostery_conditional.c_str());
			}
			return status;
		};
	bool patch_status(bool* particle_status, std::string logic) const;
	void update_active_patches(int toggle_idx);
    void init_allostery();
};

template <typename number>
class AdvAllosteryPatchyShapeParticle : public AllostericPatchyParticle<number>{
public:
	int _internal_state_size;
	bool* _internal_state;
public:
	AdvAllosteryPatchyShapeParticle(int N_patches=1 , int type = 0,int N_vertexes=0, int state_size=0)
			: AllostericPatchyParticle<number>(N_patches, type, N_vertexes) {
		_internal_state_size = state_size;
	};
	AdvAllosteryPatchyShapeParticle(const AdvAllosteryPatchyShapeParticle<number> &p) {
		this->patches = 0;
		this->_vertexes =  0;
		this->copy_from(p);
	};
	~AdvAllosteryPatchyShapeParticle() {
		delete[] _internal_state;
	};
	void copy_from(const BaseParticle<number> &p);
};


#endif /* PLPATCHYPARTICLE_H_ */
