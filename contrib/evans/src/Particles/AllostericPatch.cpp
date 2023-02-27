//
// Created by josh on 7/8/22.
//

#include "AllostericPatch.h"

/// A structure describing the patch; Each particle can have multiple patches, positioned at different places; The patches are directional and each
/// patch interacts only with its specific complementary patch;
// TODO: have this extend Patch if feasable
// TODO: unclear if this even requires to inherit BaseParticle

AllostericPatch::AllostericPatch() {
    active = false;
    bound = false;
    id = 0;
    _color = -1;
//    strength = 1;
//    a1_x = a1_y = a1_z = a2_x = a2_y = a2_z = 0;
//    locked_to_particle = -1;
//    locked_to_patch = -1;
//    locked_energy = 0;
}

AllostericPatch::AllostericPatch(LR_vector _a1_xyz, LR_vector _position, int id, int color, int state_var,
                                 int activation_var) :
                         _a1{_a1_xyz},
                         _state_var{state_var},
                         _activation_var{activation_var},
                        _position(_position),
                        id(id),
                        _color(color), /*strength(strength),*/
                        _flipped(false),
                        bound{false} {
//    a1_x = _a1_xyz.x;
//    a1_y = _a1_xyz.y;
//    a1_z = _a1_xyz.z;
//    a2_x = _a2_xyz.x;
//    a2_y = _a2_xyz.y;
//    a2_z = _a2_xyz.z;
//
//    locked_to_particle = -1;
//    locked_to_patch = -1;
//    locked_energy = 0;
}

//bool AllostericPatch::is_locked() const{
//    return locked_to_particle >= 0;
//}

int AllostericPatch::get_color() const {
    return _color;
}

//bool AllostericPatch::locked_to(int particle_id,int patch_id) const {
//    return is_locked() && (particle_id == locked_to_particle && locked_to_patch == patch_id);
//}
//bool AllostericPatch::locked_to_particle_id(int particle_id) const
//{
//    return is_locked() && particle_id == locked_to_particle;
//}
//void AllostericPatch::get_lock(int& particle_id, int& patch_id) const {
//    particle_id = locked_to_particle;
//    patch_id =  locked_to_patch;
//}
//
//number AllostericPatch::get_lock_energy(void) const
//{
//    return locked_energy;
//}

bool AllostericPatch::is_active() const {
    return this->active;
}
void AllostericPatch::set_active(bool bNewVal) {
    this->active = bNewVal;
}

bool AllostericPatch::toggle_active()
{
    if (!this->_flipped)
    {
        this->active = !this->active;
        this->_flipped = true;
    }
    return this->active;
}

int AllostericPatch::get_id() const {
    return id;
}

LR_vector AllostericPatch::a1() {
    return _a1;
}

void AllostericPatch::set_a1(LR_vector v) {
    _a1 = v;
}

int AllostericPatch::state_var() const {
    return _state_var;
}

int AllostericPatch::activation_var() const {
    return _activation_var;
}

bool AllostericPatch::is_allosterically_controlled() const {
    return activation_var() != 0;
}

bool AllostericPatch::is_allosteric_controller() const {
    return _state_var != 0;
}

LR_vector AllostericPatch::position() const {
    return _position;
}

bool AllostericPatch::flipped() const {
    return _flipped;
}

int AllostericPatch::color() const {
    return _color;
}