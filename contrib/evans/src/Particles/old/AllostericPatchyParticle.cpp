/*
 * PatchyShapeParticle.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "AllostericPatchyParticle.h"
#include "Utilities/oxDNAException.h"
#include "../../../../romano/src/Particles/PatchyShapeParticle.h"

#define HALF_ISQRT3 0.28867513459481292f


AllostericPatchyParticle::AllostericPatchyParticle(int _N_patches, int _type, int _N_vertexes) :  BaseParticle() {
	this->type = _type;
	N_patches =  _N_patches;
	N_vertexes = _N_vertexes;

	if(N_patches + N_vertexes> 0)
	{
		this->int_centers = std::vector<LR_vector>(N_patches+N_vertexes);
		this->patches = new AllostericPatch[N_patches];
		this->_vertexes = new LR_vector[N_vertexes];

	}
	else
	{
		this->patches = 0;
		this->_vertexes = 0;
	}
	//_set_base_patches();
}


AllostericPatchyParticle::~AllostericPatchyParticle() {
	delete[] patches;
	delete [] _vertexes;
}

bool operator==(const LR_vector& v1, const LR_vector& v2){
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

void AllostericPatchyParticle::copy_from(const BaseParticle &b)
{
	const AllostericPatchyParticle *bb = dynamic_cast<const AllostericPatchyParticle *>(&b);
	if (bb == 0)
	{
		throw oxDNAException("Can't convert particle to PatchyShapeParticle by dynamic cast'. Aborting");
	}

	if( ! (this->int_centers == b.int_centers && this->N_patches == bb->N_patches && this->N_vertexes == bb->N_vertexes)      )
	{
		delete [] this->patches;
		delete [] this->_vertexes;

		this->N_patches = bb->N_patches;
		this->N_vertexes = bb->N_vertexes;

		this->int_centers = std::vector<LR_vector>(bb->N_int_centers());
		patches = new AllostericPatch[bb->N_patches];
		_vertexes = new LR_vector[bb->N_vertexes];
	}

	BaseParticle::copy_from(b);

	for(int i =0 ; i < this->N_patches; i++)
	{
		// this->int_centers[i] = bb->int_centers[i];
		this->patches[i] = bb->patches[i];
	}
	for(int i =0 ; i < this->N_vertexes; i++)
	{
		// this->int_centers[i] = bb->int_centers[i];
		this->_vertexes[i] = bb->_vertexes[i];
	}
}

 void
AllostericPatchyParticle::add_patch(AllostericPatch &patch,int position) {

	if(position < 0 || position >= this->N_int_centers())
	{
		throw oxDNAException ("Could process patch id, please check that the patches of id %d are correct. Aborting",position);
	}
	patches[position] = patch;
}



void AllostericPatchyParticle::_set_base_patches() {


	for(int i = 0; i < this->N_int_centers(); i++) {
		patches[i].a1.normalize();
		patches[i].a2.normalize();
		//patches[i] *= 0.5;
	}
}


void AllostericPatchyParticle::set_positions() {
	/*
	 if(this->index == 0)
        printf("I am in set_positions for particle %d, N_patches=%d, N_vertices=%d, N_int_centers=%d \n",this->index,this->N_patches,this->N_vertexes,this->N_int_centers);

	 */
	//printf("Setting position with %f %f %f\n",this->orientationT.v1.x,this->orientationT.v1.y,this->orientationT.v1.z);
	for(int i = 0; i < this->N_patches; i++)
	{
		this->int_centers[i] = this->orientation * this->patches[i].position;
		this->patches[i].a1 = (patches[i].a1_x * this->orientationT.v1) + (patches[i].a1_y * this->orientationT.v2) + (patches[i].a1_z * this->orientationT.v3); //possibly can be accelerated
		this->patches[i].a2 = (patches[i].a2_x * this->orientationT.v1) + (patches[i].a2_y * this->orientationT.v2) + (patches[i].a2_z * this->orientationT.v3); //possibly can be accelerated

	}
	for(int i = 0; i < this->N_vertexes; i++)
	{
		this->int_centers[this->N_patches+i] = this->orientation * this->_vertexes[i];
	}

	/*
	if(this->index == 0)
	{
		for(int i = 0; i < this->N_int_centers; i++)
			{
			printf("%d center: %f %f %f vertex: %f %f %f \n",i,this->int_centers[i].x,this->int_centers[i].y,this->int_centers[i].z,this->_vertexes[i].x,this->_vertexes[i].y,this->_vertexes[i].z);
			}
	}
	 */

}


void AllostericPatchyParticle::unlock_patches(void) {
	for(int i = 0; i < this->N_patches; i++)
	{
		this->set_lock(i); //cleans the lock
	}
}


void AllostericPatchyParticle::set_lock(int patch_idx, int particle,int patch,number energy, bool ignore_refresh){
	bool state_change = this->patches[patch_idx].locked_to_particle != particle;
	if (state_change && !ignore_refresh){
		this->update_active_patches(patch_idx);
	}
	this->patches[patch_idx].locked_to_particle = particle;
	this->patches[patch_idx].locked_to_patch = patch;
	this->patches[patch_idx].locked_energy = energy;
}


void AllostericPatchyParticle::unlock(int patch_idx, bool ignore_refresh) {
	this->set_lock(patch_idx, -1, -1, 0, ignore_refresh);
}


bool AllostericPatchyParticle::locked_to_particle_id(int particle_id)
{
	for(int i = 0; i < this->N_patches; i++)
	{
		if (this->patches[i].locked_to_particle_id(particle_id) )
			return true;
	}

	return false;
}


void AllostericPatchyParticle::_set_vertexes()
{
	if(N_vertexes == 12)
	{
		_set_icosahedron_vertexes();
	}
	else {
		throw oxDNAException("Unsupported number of vertexes: %d\n",N_vertexes);
	}
}


void AllostericPatchyParticle::_set_icosahedron_vertexes() {
	double t = (1. + sqrt(5.))/2;      // golden radius
	double r = 2. * sin(2. * M_PI/5);  // radius of icosahedron of edge lenght two
	number a = 1./ (2. * r);           // will generate patches at a distance 0.5 from center
	number b = t / (2. * r);

	if (this->N_vertexes != 12)
	{
		throw oxDNAException("If you are using icosahedron, you need 12 vertices, but we have %d",this->N_vertexes);
	}

	// 12 vertexes of the icosahedron; circular
	// permutations of (0, +/-a, +/-b)
	_vertexes[ 0] = LR_vector( 0.,  a,  b);
	_vertexes[ 1] = LR_vector( 0.,  a, -b);
	_vertexes[ 2] = LR_vector( 0., -a,  b);
	_vertexes[ 3] = LR_vector( 0., -a, -b);

	_vertexes[ 4] = LR_vector(  b, 0.,  a);
	_vertexes[ 5] = LR_vector(  b, 0., -a);
	_vertexes[ 6] = LR_vector( -b, 0.,  a);
	_vertexes[ 7] = LR_vector( -b, 0., -a);

	_vertexes[ 8] = LR_vector(  a,  b, 0.);
	_vertexes[ 9] = LR_vector(  a, -b, 0.);
	_vertexes[10] = LR_vector( -a,  b, 0.);
	_vertexes[11] = LR_vector( -a, -b, 0.);

	// we now need to figure out which vertexes belong to which face
	// angle between vertexes of the same face: 1.1071487177940436, ~63.435 degrees
	// the possible angles between all pairs are: ~63.4, ~116.6 e 180.
	// each vertex has 5 vertexes at 63, 5 more at 116 and 1 at 180 (opposite)


	int nface = 0;
	number thres = 0.;  // threshold angle is 90 degrees
	for (int i = 0; i < 12; i ++) {
		for (int j = 0; j < i; j ++) {
			for (int k = 0; k < j; k ++) {
				if ((_vertexes[i]*_vertexes[j] > thres) &&
						(_vertexes[i]*_vertexes[k] > thres) &&
						(_vertexes[j]*_vertexes[k] > thres)) {
					//_faces[3*nface + 0] = i;
					//_faces[3*nface + 1] = j;
					//_faces[3*nface + 2] = k;
					/*printf ("\n\n%d %d %d @\n", i, j, k);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[i].x, _vertexes[i].y, _vertexes[i].z);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[j].x, _vertexes[j].y, _vertexes[j].z);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[k].x, _vertexes[k].y, _vertexes[k].z);
						printf ("  %g %g %g\n", 4.f*(_vertexes[i]*_vertexes[j]), 4.f*(_vertexes[i]*_vertexes[k]), 4.*(_vertexes[j]*_vertexes[k]));*/
					nface ++;
				}
			}
		}
	}





	this->set_positions();
}

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


bool SimpleAllosteryPatchyShapeParticle::patch_status(bool* particle_status, std::string logic) const{
	if (logic == "true"){
		return true;
	}
	else if (logic == "false") {
		return false;
	}
	//	int paren_count = 0;
	int paren_count = 0;
	std::string::iterator paren_start;
	std::string numstr;
	bool prefix = true;
	bool negate_flag = true; // default to no negation
	char op = 0;
	for (std::string::iterator it = logic.begin(); it != logic.end(); ++it) {
		if (*it == '('){
			if (paren_count == 0)
			{
				paren_start = it; //should invoke copy constructor
			}
			paren_count++;
		}
		else if (*it == ')') {
			paren_count--;
			if (paren_count == 0)
			{
				std::string subexpr(paren_start + 1, it);
				bool paren_statement_val = this->patch_status(particle_status, subexpr) == negate_flag;
				parse_boolean_statement(prefix, paren_statement_val, op);
			}
		}
		else if (paren_count == 0) { //if program is mid-parentheses, continue until it finds a close-paren
			if (*it == '!'){
				negate_flag = false;
				continue;

			}
			if (int('0') <= *it && *it <= int('9')) {
				numstr += *it;
				continue;
			}
			if (*it == '&' || *it == '|') {

				op = *it;
				//deliberately no "else" here; evaluating numstr can coexist w/ operators
			}
			if (numstr != ""){ // either whitespace or an operator terminate a logical statement
				int patch = stoi(numstr);
				if (patch >= this->N_patches){
					throw oxDNAException("Patch index " + numstr + " out of bounds.");
				}
				bool val = particle_status[patch] == negate_flag;

				parse_boolean_statement(prefix, val, op);

				// clear numstr and negate_flag
				negate_flag = true;
				numstr = "";
			}
		}
	}

	//parse suffix
	if (numstr != ""){ // either whitespace or an operator terminate a logical statement
		int patch = stoi(numstr);
		if (patch >= this->N_patches){
			throw oxDNAException("Patch index " + numstr + " out of bounds.");
		}
		bool val = particle_status[patch] == negate_flag;

		parse_boolean_statement(prefix, val, op);
	}

	return prefix;
}


// WARNING: the array returned by this method allocates memory, which must be deallocated!
bool* AllostericPatchyParticle::get_binding_state() const {
	bool* particle_status = new bool[this->N_patches];
	for (int i = 0; i < this->N_patches; i++)
	{
		particle_status[i] = this->patches[i].is_locked();
	}
	return particle_status;
}


void SimpleAllosteryPatchyShapeParticle::update_active_patches(int toggle_idx){
	bool* particle_status = this->get_binding_state();
	ParticleStateChange change(particle_status, this->N_patches, toggle_idx);

	//DEBUGGING
	//	for (std::unordered_map<ParticleStateChange, std::vector<int>>::iterator it = this->allostery_map->begin(); it != this->allostery_map->end(); ++it)
	//	{
	//		if (*it == change)
	//		{
	//
	//		}
	//	}
	std::vector<int> updates = (*this->allostery_map)[change];

	//#ifdef DEBUG
	std::string flips = "[";
	//#endif
	for (std::vector<int>::iterator it = updates.begin(); it != updates.end(); ++it)
	{
		bool a_before = this->patches[*it].active;
		bool a_after = this->patches[*it].toggle_active();
		flips += std::to_string(*it) + ":" + (a_before ? "A" : "!A") + "->" + (a_after ? "A" : "NA") + ",Flp=" + std::to_string(this->patches[*it].flipped);


	}
	//#ifdef DEBUG
	flips += "]";
	std::string status_before_str = "[";
	std::string status_after_str = "[";
	std::string new_activations = "[";
	for (int i = 0; i < this->N_patches; i++) {
		status_before_str += (particle_status[i] ? " T" : " F");
		status_after_str += (particle_status[i] != (i == toggle_idx) ? " T" : " F"); //right?
		new_activations += std::to_string(i) + ":" + (this->patches[i].active ? "T " : "F ");
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

void SimpleAllosteryPatchyShapeParticle::copy_from(const BaseParticle &b){
	// pass by reference since copies of particles will inherit the base particle type's conditionals
	// should save memory, especially in large simulations with particles with many patches
	const SimpleAllosteryPatchyShapeParticle *bb = dynamic_cast<const SimpleAllosteryPatchyShapeParticle *>(&b);
	if (bb == 0)
	{
		throw oxDNAException("Can't convert particle to PatchyShapeParticle by dynamic cast'. Aborting");
	}
    AllostericPatchyParticle::copy_from(b);
	allostery_map = bb->allostery_map;
}

void SimpleAllosteryPatchyShapeParticle::init_allostery() {
    std::unordered_map<ParticleStateChange, std::vector<int>>* allosteric_control = new std::unordered_map<ParticleStateChange, std::vector<int>>();
    // construct allosteric control
    for (int i = 0; i < pow(2, this->N_patches); i++){
        for (int iSwap = 0; iSwap < this->N_patches; iSwap++)
        {
            bool* key = new bool[this->N_patches];
            for (int iPatch = 0; iPatch < this->N_patches; iPatch++) {
                key[iPatch] = i % int(pow(2, this->N_patches - iPatch)) < pow(2, this->N_patches - iPatch - 1);
            }
            ParticleStateChange* change = new ParticleStateChange(key, this->N_patches, iSwap);
            std::vector<int> affected_patches;
            bool key_after[this->N_patches];
            // REALLY feels like there's a better way of doing this
            for (int j = 0; j < this->N_patches; j++) {
                key_after[j] = key[j];
            }
            key_after[iSwap] = !key_after[iSwap];

            // loop through the patches on this particle
            for (int iPatch = 0; iPatch < this->N_patches; iPatch++){
                // if the patch
                bool patch_status_before = this->patch_status(key, iPatch);
                bool patch_status_after = this->patch_status(key_after, iPatch);
                if (patch_status_before != patch_status_after){
                    affected_patches.push_back(iPatch);
                }
            }
            if (affected_patches.size() > 0)
            {
                (*allosteric_control)[*change] = affected_patches;
            }
            else {
                delete change;
            }
        }
    }

    //initialize state
    bool default_state[this->N_patches];
    for (int i = 0; i < this->N_patches; i++){
        default_state[i] = false; // default state = no bonding
    }
    for (int i = 0; i < this->N_patches; i++){
        this->patches[i].active = this->patch_status(default_state, i);
    }
    this->allostery_map = allosteric_control;
}

void AdvAllosteryPatchyShapeParticle::copy_from(const BaseParticle &b) {
	const AdvAllosteryPatchyShapeParticle *bb = dynamic_cast<const AdvAllosteryPatchyShapeParticle *>(&b);
	if (bb == 0)
	{
		throw oxDNAException("Can't convert particle to PatchyShapeParticle by dynamic cast'. Aborting");
	}
    AllostericPatchyParticle::copy_from(b);
	this->_internal_state_size = bb->_internal_state_size;
	delete [] this->_internal_state;
	_internal_state = new bool[_internal_state_size];
	for (int i = 0; i < _internal_state_size; i++){
		_internal_state[i] = bb->_internal_state[i];
	}
}
