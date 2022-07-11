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

AllostericPatch::AllostericPatch(LR_vector _a1_xyz, LR_vector _a2_xyz, LR_vector _position, int _id, int _color,
                                 bool _active, std::string allostery_conditional, bool activation_reversible) :
                         _init_a1{_a1_xyz},
                         _init_a2{_a2_xyz},
                        _position(_position),
                        id(_id), active(_active),
                        _color(_color), /*strength(strength),*/
                        allostery_conditional(allostery_conditional),
                        activation_reversible(activation_reversible),
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
    if (this->activation_reversible || !this->_flipped)
    {
        this->active = !this->active;
        this->_flipped = true;
    }
    return this->active;
}

LR_vector AllostericPatch::a1() {
    return _a1;
}

LR_vector AllostericPatch::a2() {
    return _a2;
}

void AllostericPatch::set_a1(LR_vector v) {
    _a1 = v;
}

void AllostericPatch::set_a2(LR_vector v) {
    _a2 = v;
}

LR_vector AllostericPatch::position() const {
    return _position;
}

LR_vector AllostericPatch::a1_initial() {
    return _init_a1;
}

LR_vector AllostericPatch::a2_initial() {
    return _init_a2;
}

bool AllostericPatch::flipped() const {
    return _flipped;
}

int AllostericPatch::color() const {
    return _color;
}