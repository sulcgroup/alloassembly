//
// Created by josh on 7/8/22.
//

#include "AllostericPatchyParticle.h"

/**
 * elementwise comparison of two LR_vector objects. I'm honestly a bit shocked this function doesn't already exist
 * Or maybe it does, but idk where\
 * There may be an issue with floating-point precision; I adivse against relying on
 * this function too strongly
 * @param v1
 * @param v2
 * @return
 */
bool operator==(const LR_vector& v1, const LR_vector& v2) {
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

AllostericPatchyParticle::AllostericPatchyParticle(int _N_patches, int _type) : BaseParticle() {
    this->type = _type;

    this->patches.resize(_N_patches);
    this->int_centers.resize(_N_patches);
    //_set_base_patches();
}


AllostericPatchyParticle::~AllostericPatchyParticle() {
    // DO NOT DELETE MAP POINTERS!!!
    // those were allocated elsewhere and will be deallocated elsewhere
}

void AllostericPatchyParticle::copy_from(const BaseParticle &b)
{
    const AllostericPatchyParticle *bb = dynamic_cast<const AllostericPatchyParticle *>(&b);
    if (bb == 0)
    {
        throw oxDNAException("Can't convert particle to PatchyShapeParticle by dynamic cast'. Aborting");
    }

    BaseParticle::copy_from(b);

    for(int i =0 ; i < this->n_patches(); i++)
    {
        // this->int_centers[i] = bb->int_centers[i];
        this->patches[i] = bb->patches[i];
    }
    this->int_centers = b.int_centers;
    _transitionMap = bb->_transitionMap;
    _state_size = bb->state_size();
    _state = bb->get_state();
}

void AllostericPatchyParticle::add_patch(AllostericPatch &patch, int position) {
    LR_vector patch_int_center = patch.position();
    // set interaction center to patch position
    int_centers[position] = patch_int_center;
    if(position < 0 || position > int_centers.size())
    {
        throw oxDNAException ("Could not process patch id, please check that the patches of id %d are correct. Aborting",position);
    }
    patches[position] = patch;
}



void AllostericPatchyParticle::_set_base_patches() {
    for(int i = 0; i < this->N_int_centers(); i++) {
        patches[i].a1().normalize();
    }
}


void AllostericPatchyParticle::set_positions() {
    /*
     if(this->index == 0)
        printf("I am in set_positions for particle %d, N_patches=%d, N_vertices=%d, N_int_centers=%d \n",this->index,this->N_patches,this->N_vertexes,this->N_int_centers);

     */
    //printf("Setting position with %f %f %f\n",this->orientationT.v1.x,this->orientationT.v1.y,this->orientationT.v1.z);
    for(int i = 0; i < this->n_patches(); i++)
    {
        this->int_centers[i] = this->orientation * this->patches[i].position();
//        this->patches[i].set_a1(patches[i].a1_x * this->orientationT.v1) + (patches[i].a1_y * this->orientationT.v2) + (patches[i].a1_z * this->orientationT.v3); //possibly can be accelerated
//        this->patches[i].set_a2(patches[i].a2_x * this->orientationT.v1) + (patches[i].a2_y * this->orientationT.v2) + (patches[i].a2_z * this->orientationT.v3); //possibly can be accelerated

    }
}


//void AllostericPatchyParticle::unlock_patches(void) {
//    for(int i = 0; i < this->N_patches; i++)
//    {
//        this->set_lock(i); //cleans the lock
//    }
//}
//
//
//void AllostericPatchyParticle::set_lock(int patch_idx, int particle,int patch,number energy, bool ignore_refresh) {
//    bool state_change = this->patches[patch_idx].locked_to_particle != particle;
//    if (state_change && !ignore_refresh){
//        this->update_active_patches(patch_idx);
//    }
//    this->patches[patch_idx].locked_to_particle = particle;
//    this->patches[patch_idx].locked_to_patch = patch;
//    this->patches[patch_idx].locked_energy = energy;
//}
//
//
//void AllostericPatchyParticle::unlock(int patch_idx, bool ignore_refresh) {
//    this->set_lock(patch_idx, -1, -1, 0, ignore_refresh);
//}
//
//
//bool AllostericPatchyParticle::locked_to_particle_id(int particle_id)
//{
//    for(int i = 0; i < this->N_patches; i++)
//    {
//        if (this->patches[i].locked_to_particle_id(particle_id) )
//            return true;
//    }
//
//    return false;
//}

//
//void AllostericPatchyParticle::_set_vertexes()
//{
//    if(N_vertexes == 12)
//    {
//        _set_icosahedron_vertexes();
//    }
//    else {
//        throw oxDNAException("Unsupported number of vertexes: %d\n",N_vertexes);
//    }
//}

//
//void AllostericPatchyParticle::_set_icosahedron_vertexes() {
//    double t = (1. + sqrt(5.))/2;      // golden radius
//    double r = 2. * sin(2. * M_PI/5);  // radius of icosahedron of edge lenght two
//    number a = 1./ (2. * r);           // will generate patches at a distance 0.5 from center
//    number b = t / (2. * r);
//
//    if (this->N_vertexes != 12)
//    {
//        throw oxDNAException("If you are using icosahedron, you need 12 vertices, but we have %d",this->N_vertexes);
//    }
//
//    // 12 vertexes of the icosahedron; circular
//    // permutations of (0, +/-a, +/-b)
//    _vertexes[ 0] = LR_vector( 0.,  a,  b);
//    _vertexes[ 1] = LR_vector( 0.,  a, -b);
//    _vertexes[ 2] = LR_vector( 0., -a,  b);
//    _vertexes[ 3] = LR_vector( 0., -a, -b);
//
//    _vertexes[ 4] = LR_vector(  b, 0.,  a);
//    _vertexes[ 5] = LR_vector(  b, 0., -a);
//    _vertexes[ 6] = LR_vector( -b, 0.,  a);
//    _vertexes[ 7] = LR_vector( -b, 0., -a);
//
//    _vertexes[ 8] = LR_vector(  a,  b, 0.);
//    _vertexes[ 9] = LR_vector(  a, -b, 0.);
//    _vertexes[10] = LR_vector( -a,  b, 0.);
//    _vertexes[11] = LR_vector( -a, -b, 0.);
//
//    // we now need to figure out which vertexes belong to which face
//    // angle between vertexes of the same face: 1.1071487177940436, ~63.435 degrees
//    // the possible angles between all pairs are: ~63.4, ~116.6 e 180.
//    // each vertex has 5 vertexes at 63, 5 more at 116 and 1 at 180 (opposite)
//
//
//    int nface = 0;
//    number thres = 0.;  // threshold angle is 90 degrees
//    for (int i = 0; i < 12; i ++) {
//        for (int j = 0; j < i; j ++) {
//            for (int k = 0; k < j; k ++) {
//                if ((_vertexes[i]*_vertexes[j] > thres) &&
//                    (_vertexes[i]*_vertexes[k] > thres) &&
//                    (_vertexes[j]*_vertexes[k] > thres)) {
//                    //_faces[3*nface + 0] = i;
//                    //_faces[3*nface + 1] = j;
//                    //_faces[3*nface + 2] = k;
//                    /*printf ("\n\n%d %d %d @\n", i, j, k);
//                        printf ("%7.5g %7.5g %7.5g\n", _vertexes[i].x, _vertexes[i].y, _vertexes[i].z);
//                        printf ("%7.5g %7.5g %7.5g\n", _vertexes[j].x, _vertexes[j].y, _vertexes[j].z);
//                        printf ("%7.5g %7.5g %7.5g\n", _vertexes[k].x, _vertexes[k].y, _vertexes[k].z);
//                        printf ("  %g %g %g\n", 4.f*(_vertexes[i]*_vertexes[j]), 4.f*(_vertexes[i]*_vertexes[k]), 4.*(_vertexes[j]*_vertexes[k]));*/
//                    nface ++;
//                }
//            }
//        }
//    }
//
//    this->set_positions();
//}

//helper function
void parse_boolean_statement(bool &status, bool operand, char &op){ //0/10 function name
    if (op == 0){
        status = operand;
    }
    else if (op == '&'){
        status &= operand;
    }
    else if (op == '|'){
        status |= operand;
    }
    else {
        throw oxDNAException("Malformed logic.");
    }
    op = 0; //reset operation
}

void AllostericPatchyParticle::set_patch_bound(int pidx) {
    patches[pidx].bound = true;
    if (patches[pidx].is_allosteric_controller()){
        int var = patches[pidx].state_var();
        // set var to true
        _state |= 1 << var;
    }
}

void AllostericPatchyParticle::update_active_patches(int oldState){
    int newState = get_state();

    //DEBUGGING
    //	for (std::unordered_map<ParticleStateChange, std::vector<int>>::iterator it = this->allostery_map->begin(); it != this->allostery_map->end(); ++it)
    //	{
    //		if (*it == change)
    //		{
    //
    //		}
    //	}
    std::set<int> updates = (*this->_updateMap)[state_size() * oldState + newState];
    if (updates.size() == 0){
        return;
    }

    //#ifdef DEBUG
    std::string flips = "[";
    //#endif
    bool state_changed = false;
    for (std::set<int>::iterator it = updates.begin(); it != updates.end(); ++it)
    {
        bool a_before = this->patches[*it].is_active();
        bool a_after = this->patches[*it].toggle_active();
        state_changed |= a_before != a_after;
        flips += std::to_string(*it) + ":" + (a_before ? "A" : "!A") + "->" + (a_after ? "A" : "NA") + ",Flp=" + std::to_string(this->patches[*it].flipped()) + ";";
    }
    if (!state_changed)
        return;
    //#ifdef DEBUG
    flips += "]";
    std::string status_before_str = "[";
    std::string status_after_str = "[";
    std::string new_activations = "[";
    for (int i = 0; i < this->state_size(); i++) {
        status_before_str += (GET_BIT(get_state(), i) ? " T" : " F");
        status_after_str += (GET_BIT(get_state(), i) != (i == oldState) ? " T" : " F"); //right?
        new_activations += std::to_string(i) + ":" + (this->patches[i].is_active() ? "T " : "F ");
    }
    status_before_str += "]";
    status_after_str += "]";
    new_activations += "]";
    OX_LOG(Logger::LOG_INFO, "Particle with ID %d changed internal state from %s to %s, affecting patches %s. \nNew patch activations: %s.",
           this->index,
           status_before_str.c_str(),
           status_after_str.c_str(),
           flips.c_str(),
           new_activations.c_str());
    //#endif
    //don't need to delete particle_status b/c that will be done automatically when change goes out of scope
}

/**
 * Initializes this particle's allosteric control
 * @param transitionMap a pointer to an ALREADY GENERATED StateTransitionMap object which will be used to gen allsotery
 * @param updateMap a pointer to an EMPTY StateTransitionMap which this method will populate
 */
void AllostericPatchyParticle::init_allostery(int state_size, StateTransitionMap *transitionMap,
                                              ActivationUpdateMap *updateMap) {
    _state_size = state_size;
    // set map pointers
    _transitionMap = transitionMap;
    _updateMap = updateMap;

    // generate update map
    for (int beforeState = 0; beforeState < state_size; beforeState++){
        for (int afterState = 0; afterState < state_size; afterState++){
            int idx = beforeState * state_size + afterState;
            for (int pidx = 0; pidx < n_patches(); pidx++){
                if (patches[pidx].is_allosterically_controlled()){
                    // get absolute value of control idx
                    // doesn't really matter whether it's positive or negative - only whether
                    // the bit at that position changes
                    int controlidx = abs(patches[pidx].activation_var());
                    // if the bit at the control idx changes
                    if (GET_BIT(beforeState, controlidx) != GET_BIT(afterState, controlidx)){
                        (*_updateMap)[idx].emplace(pidx); // add to list of flips
                    }
                }
            }
        }
    }

    //initialize state
    _state = 0;

    for (int i = 0; i < this->n_patches(); i++){
        if (this->patches[i].is_allosterically_controlled()){
            int statevar = patches[i].activation_var();
            // if statevar is negative, patch starts active. else, patch starts inactive
            this->patches[i].set_active(statevar < 0);
        }
        else {
            this->patches[i].set_active(true);
        }
    }
}

int AllostericPatchyParticle::do_state_transition(int rand_idx) {
    unsigned int oldState = get_state();
    _state = (*_transitionMap)[oldState][rand_idx];
    if (_state != oldState){
        update_active_patches(oldState);
    }
}

int AllostericPatchyParticle::n_patches() const {
    return patches.size();
}

unsigned int AllostericPatchyParticle::get_state() const {
    return _state;
}

void AllostericPatchyParticle::set_state(unsigned int new_state) {
    _state =  new_state; // set state
    // update activations
    for (int p = 0; p < n_patches(); p++){
        int activator = patches[p].activation_var();
        bool active;
        if (activator == 0){
            active = true;
        }
        else if (activator > 0){
            active = GET_BIT(new_state, activator);
        }
        else {
            active = !GET_BIT(new_state, -activator);
        }
        // update patch activation state
        patches[p].set_active(active);
        // if the patch has a binding state var, update binding state to reflect
        if (patches[p].state_var() != 0) {
            patches[p].bound = GET_BIT(new_state, patches[p].state_var());
        }
        // nothing to do otherwise
    }
}

/**
 * @return the size of the particle's state, in # bits
 */
int AllostericPatchyParticle::state_size() const {
    return _state_size;
}

/**
 * @return the number of possible states, which is 2^state_size
 */
int AllostericPatchyParticle::n_states() const {
    return 1 << state_size();
}

