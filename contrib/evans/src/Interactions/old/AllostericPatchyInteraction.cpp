#include "AllostericPatchyInteraction.h"
#include "Utilities/Utils.h"
#include "Interactions/InteractionUtils.h"
#include <sstream>
#include "../../../../romano/src/Interactions/PatchyShapeInteraction.h"
#include "Utilities/parse_input/parse_input.h"


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

// implementation of inline functions:
bool AllostericPatchyInteraction::_patches_compatible(AllostericPatchyParticle  *p, AllostericPatchyParticle  *q, int pi, int pj ) {
	if(abs(p->patches[pi].color)  < 10 &&  abs(q->patches[pj].color ) < 10  )  //we are in the self-complementary regime
	{
		if(p->patches[pi].color == q->patches[pj].color ) { //patches are the same, hence complementary
			return true;
		}
		else {
			return false;
		}
	}
	else if (p->patches[pi].color + q->patches[pj].color == 0) { //the case where only complementary patches (like -10 and +10) can bind
		return true;

	}

	//if we got here, no compatibility
	return false;

}


// implementation of inline functions:
inline bool AllostericPatchyInteraction::_bonding_allowed(AllostericPatchyParticle  *p, AllostericPatchyParticle  *q, int pi, int pj )
{
	// printf("@@@@\n");
	bool allowed = false;
	//if(p->patches[pi].color == q->patches[pj].color && p->patches[pi].active && q->patches[pj].active ) //patches are complementary
	if(_patches_compatible(p,q,pi,pj) && p->patches[pi].active && q->patches[pj].active ) //patches are complementary
	{
		if(this->_same_type_bonding == false && p->type == q->type)
		{
			allowed =  false;
		}
		else allowed =  true;

		if(this->_no_multipatch) //one patch binds only one other patch
		{
			if ( p->patches[pi].locked_to(q->index,pj) )
			{
				if (  q->patches[pj].locked_to(p->index,pi)  )
				{
					return true;
				}
				else
				{
					throw oxDNAException("Assymetric lock detected for particles %d %d",p->index,q->index);
				}
			}
			else if (  ! p->patches[pi].is_locked() && !q->patches[pj].is_locked())
			{
				return true;
			}
			else {
				return false;
			}
		}
	}
	else
	{
		allowed = false;
	}


	return allowed;


	//int i1 = p->patches[pi].id;
	//int i2 = q->patches[pj].id;

	/*
	if(allowed != (bool)this->_interaction_patch_types[i1 * this->_N_patch_types + i2] )
	{
		throw oxDNAException("_bonding_allowed: incompatible values");
	}
	 */

	// return ( (bool)this->_interaction_patch_types[i1 * this->_N_patch_types + i2] && p->patches[pi].active &&  q->patches[pj].active   );
}



number AllostericPatchyInteraction:: _V_mod(int type, number t)
{
	number val = (number) 0;
	t -= PLPATCHY_THETA_T0[type];
	if (t < 0)
		t *= -1;

	if (t < PLPATCHY_THETA_TC[type]) {
		if (t > PLPATCHY_THETA_TS[type]) {
			// smoothing
			val = PLPATCHY_THETA_B[type] * SQR(PLPATCHY_THETA_TC[type] - t);
		} else
			val = (number) 1.f - PLPATCHY_THETA_A[type] * SQR(t);
	}

	return val;
}



number AllostericPatchyInteraction:: _V_modD(int type, number t)
{
	number val = (number) 0;
	number m = (number) 1;
	t -= PLPATCHY_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(t < 0) {
		t *= -1;
		m = (number) -1;
	}

	if(t < PLPATCHY_THETA_TC[type]) {
		if(t > PLPATCHY_THETA_TS[type]) {
			// smoothing
			val = m * 2 * PLPATCHY_THETA_B[type] * (t - PLPATCHY_THETA_TC[type]);
		}
		else val = -m * 2 * PLPATCHY_THETA_A[type] * t;
	}

	return val;
}



number AllostericPatchyInteraction:: _V_modDsin(int type, number t)
{
	number val = (number) 0;
	number m = (number) 1;
	number tt0 = t - PLPATCHY_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(tt0 < 0) {
		tt0 *= -1;
		m = (number) -1;
	}

	if(tt0 < PLPATCHY_THETA_TC[type]) {
		number sint = sin(t);
		if(tt0 > PLPATCHY_THETA_TS[type]) {
			// smoothing
			val = m * 2 * PLPATCHY_THETA_B[type] * (tt0 - PLPATCHY_THETA_TC[type]) / sint;
		}
		else {
			if(SQR(sint) > 1e-8) val = -m * 2 * PLPATCHY_THETA_A[type] * tt0 / sint;
			else val = -m * 2 * PLPATCHY_THETA_A[type];
		}
	}

	return val;
}




number AllostericPatchyInteraction::_repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	sigma *= (2*this->_sphere_radius);
	rstar *= (2*this->_sphere_radius);
	b /= SQR(2*this->_sphere_radius);
	rc *= (2*this->_sphere_radius);

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = PLEXCL_EPS * b * SQR(rrc);
			if(update_forces) force = -r * (2 * PLEXCL_EPS * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * PLEXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * PLEXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}



number AllostericPatchyInteraction::_repulsive_lj_n(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, int n, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = PLEXCL_EPS * b * SQR(rrc);
			if(update_forces) force = -r * (2 * PLEXCL_EPS * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = std::pow(tmp,n/2); //tmp * tmp * tmp;
			energy = 4 * PLEXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (4 * n * PLEXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}



number AllostericPatchyInteraction::_exc_LJ_vol_interaction_sphere(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {

	LR_vector force(0,0,0);


	number energy = _repulsive_lj(*r, force, PLEXCL_S, PLEXCL_R, PLEXCL_B, PLEXCL_RC, update_forces);
	/*
	for(double x = 0.7; x < 1.001; x += 0.001)
	{
		LR_vector gg (x,0,0);
		energy =  _repulsive_lj(gg, force, PLEXCL_S, PLEXCL_R, PLEXCL_B, PLEXCL_RC, false);
		printf("%f %f @\n",x,energy);
	}
	exit(-1);
	 */

	if(update_forces && energy != 0) {
		//torquep = -p->int_centers[DNANucleotide::BASE].cross(force);
		///torqueq = q->int_centers[DNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	//printf("Calculated EXC energy between %d and %d and obtained %f, distance is %f, positions are (%f,%f,%f) : (%f,%f,%f) \n",p->index,q->index ,energy,sqrt(r->norm()), p->pos.x,p->pos.y,p->pos.z, q->pos.x,q->pos.y,q->pos.z );

	return energy;
}



number AllostericPatchyInteraction::_exc_vol_hard_icosahedron(BaseParticle *ap, BaseParticle *aq, LR_vector *r,bool update_forces) {
	if (update_forces)
		throw oxDNAException("No forces, figlio di ndrocchia");

	PatchyShapeParticle *p = NULL, *q = NULL;
	LR_vector my_r;
	if (ap->index < aq->index) {
		p = static_cast<PatchyShapeParticle *>(ap);
		q = static_cast<PatchyShapeParticle *>(aq);
		my_r = *r;
	} else {
		p = static_cast<PatchyShapeParticle *>(aq);
		q = static_cast<PatchyShapeParticle *>(ap);
		my_r = this->_box->min_image(p, q);
	}

	number rnorm = my_r.norm();
	if (rnorm > this->_sqr_rcut)
		return (number) 0.f;

	if (rnorm > (number) 1.)
		return (number) 0.f;

	// radius of inscribed sphere: sqrt((1./12.) + (1./(6.*sqrt(5.))))
	if (rnorm < _tworinscribed * _tworinscribed) {
		this->set_is_infinite(true);
		return (number) 1.0e12;
	}

	// we look for the three vertexes that are the closest to
	// the distance on each particle
	number max1 = (number) -2.0f;
	number max2 = (number) -2.0f;
	number max11 = (number) -3.0f;
	number max12 = (number) -4.0f;
	number max21 = (number) -3.0f;
	number max22 = (number) -4.0f;
	int v1 = -1, v11 = -1, v21 = -1;
	int v2 = -1, v12 = -1, v22 = -1;

	// finds the three closest vertexes.
	int OFF = p->N_patches;
	for (int k = 0; k < 6; k ++) {
		number a =  (p->int_centers[OFF + k] * (my_r));
		number b = -(q->int_centers[OFF + k] * (my_r));
		if (a > max1) {
			max12 = max11; max11 = max1; max1 = a;
			v12 = v11; v11 = v1; v1 = k;
		} else if (a > max11) {
			max12 = max11; max11 = a;
			v12 = v11; v11 = k;
		} else if (a > max12) {
			max12 = a;
			v12 = k;
		}
		if (b > max2) {
			max22 = max21; max21 = max2; max2 = b;
			v22 = v21; v21 = v2; v2 = k;
		} else if (b > max21) {
			max22 = max21; max21 = b;
			v22 = v21; v21 = k;
		} else if (b > max22) {
			max22 = b;
			v22 = k;
		}
		int kk = k + 6;
		a = -a;
		b = -b;
		if (a > max1) {
			max12 = max11; max11 = max1; max1 = a;
			v12 = v11; v11 = v1; v1 = kk;
		} else if (a > max11) {
			max12 = max11; max11 = a;
			v12 = v11; v11 = kk;
		} else if (a > max12) {
			max12 = a;
			v12 = kk;
		}
		if (b > max2) {
			max22 = max21; max21 = max2; max2 = b;
			v22 = v21; v21 = v2; v2 = kk;
		} else if (b > max21) {
			max22 = max21; max21 = b;
			v22 = v21; v21 = kk;
		} else if (b > max22) {
			max22 = b;
			v22 = kk;
		}
	}

	// Early exit; this is the Separating Axis Theorem, using
	// the distance as the candidate axis. If we find light, it means
	// there cannot be an ovelap
	// number dd = rmod - max2/(0.5*rmod) - max1/(0.5*rmod);
	number dd = rnorm - max2 / 0.5 - max1 / 0.5;
	if (dd > (number) 0.)
		return 0.;

	// we try the two closest faces to see if there is separation
	LR_vector pf = ((number)(1./3.))*(p->int_centers[OFF + v1] + p->int_centers[OFF + v11] + p->int_centers[OFF + v12]);

	bool would_exit = false;
	number rin = 0.5 * _tworinscribed;

	LR_vector d = pf / rin; // now normalized
	number qpf = (my_r + q->int_centers[OFF + v2] - pf) * d;
	number tmp = (my_r + q->int_centers[OFF + v21] - pf) * d;
	if (tmp < qpf) qpf = tmp;
	tmp = (my_r + q->int_centers[OFF + v22] - pf) * d;
	if (tmp < qpf) qpf = tmp;
	// we use a safe threshold of 0.015 to account for non-orthonormal rotation matrixes
	if (qpf > (number) 0.015) would_exit = true;
	if (would_exit) return (number) 0.;

	LR_vector qf = ((number)(1./3.))*(q->int_centers[OFF + v2] + q->int_centers[OFF + v21] + q->int_centers[OFF + v22]);
	d = qf / rin;
	number ppf = (p->int_centers[OFF + v1] - my_r - qf) * d;
	tmp = (p->int_centers[OFF + v11] - my_r - qf) * d;
	if (tmp < ppf) ppf = tmp;
	tmp = (p->int_centers[OFF + v12] - my_r - qf) * d;
	if (tmp < ppf) ppf = tmp;
	// we use a safe threshold of 0.015 to account for non-orthonormal rotation matrixes
	if (ppf > (number) 0.015) would_exit = true;
	if (would_exit) return (number) 0.;

	// check edges from p with faces from q
	LR_vector s1, s2;
	s1 = p->int_centers[v1+OFF] - my_r;
	for (int i = 0; i < 5; i++) { // for each edge startging from p->int_centers[v1]
		s2 = p->int_centers[OFF+_close_vertexes[5 * v1 + i]] - my_r;
		for (int j = 0; j < 5; j++) { // for each face of q
			bool check = InteractionUtils::edge_triangle_intersection((s1), // starting point
					(s2),													// end point
					(q->int_centers[OFF+v2]),					           	// first vertex of triangle
					(q->int_centers[OFF+_close_vertexes[5 * v2 + j]]), 		// second vertex of triangle
					(q->int_centers[OFF+_close_vertexes[5 * v2 + ((j + 1) % 5)]]) // third vertex of triangle
			);
			if (check) {
				this->set_is_infinite(true);
				return (number) 1.0e12;
			}
		}
	}

	s1 = q->int_centers[OFF+v2] + my_r;
	for (int i = 0; i < 5; i++) { // for each edge startging from q->int_centers[v2]
		s2 = q->int_centers[OFF+_close_vertexes[5 * v2 + i]] + my_r;
		for (int j = 0; j < 5; j++) { // for each face of p
			bool check = InteractionUtils::edge_triangle_intersection((s1), // starting point
					(s2),                                           // end point
					(p->int_centers[OFF+v1]),	         // first vertex of triangle
					(p->int_centers[OFF+_close_vertexes[5 * v1 + j]]), // second vertex of triangle
					(p->int_centers[OFF+_close_vertexes[5 * v1 + ((j + 1) % 5)]]) // third vertex of triangle
			);
			if (check) {
				this->set_is_infinite(true);
				return (number) 1.0e12;
			}
		}
	}

	return (number) 0.;
}



number AllostericPatchyInteraction::_exc_vol_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {

	if(this->_shape == SPHERE_SHAPE) {
		return this->_exc_LJ_vol_interaction_sphere(p,q,r,update_forces);
	}
	else if (this->_shape == ICOSAHEDRON_SHAPE) {
		return this->_exc_vol_hard_icosahedron(p,q,r,update_forces);
	} else {
		throw oxDNAException("Selected particle shape type not supported");
	}
}

//USING THIS ONE!
number AllostericPatchyInteraction::_patchy_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	number rnorm = r->norm();
	if(rnorm > this->_sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

	/* repulsion energy, now moved to exc vol
	number part = 1.0f / powf(rnorm, PATCHY_POWER * 0.5f);
	energy = part - _E_cut;

	if(update_forces) {
		LR_vector force = *r * (PATCHY_POWER * part / rnorm);
		p->force -= force;
		q->force += force;
	}

	 */
	//printf("Particles %d and %d: distance %f  repulsion ene: %f\n",p->index,q->index,sqrt(rnorm),energy);

	AllostericPatchyParticle *pp = static_cast<AllostericPatchyParticle *>(p);
    AllostericPatchyParticle *qq = static_cast<AllostericPatchyParticle *>(q);

	int c = 0;
	LR_vector tmptorquep(0, 0, 0);
	LR_vector tmptorqueq(0, 0, 0);
	for(int pi = 0; pi < pp->N_patches; pi++) {
		LR_vector ppatch = p->int_centers[pi];

		for(int pj = 0; pj < qq->N_patches; pj++) {
			//printf("Patches %d and %d , colors %d %d \n",pi,pj,pp->patches[pi].color,qq->patches[pj].color);

			//if(pp->patches[pi].color == qq->patches[pj].color && pp->patches[pi].active && qq->patches[pj].active) //patches are complementary
			if(this->_bonding_allowed(pp,qq,pi,pj)  )
			{

				number K = pp->patches[pi].strength;
				LR_vector qpatch = q->int_centers[pj];

				LR_vector patch_dist = *r + qpatch - ppatch;
				number dist = patch_dist.norm();
				LR_vector patch_dist_dir = patch_dist / sqrt(dist);
				number rdist = sqrt(rnorm);
				LR_vector r_dist_dir = *r / rdist;

				//printf("Patches %d and %d distance %f  cutoff is: %f,\n",pp->patches[pi].id,qq->patches[pj].id,dist,SQR(PATCHY_CUTOFF));

				if(dist < SQR(PATCHY_CUTOFF)) {
					//printf("CRITICAL CALCULATING FORCE BETWEEN %d %d",q->index,p->index);
					c++;
					number energy_ij = 0;
					//distance part of attractive interaction
					number r8b10 = dist*dist*dist*dist / _patch_pow_alpha;
					number exp_part = -1.001f * exp(-(number)0.5f * r8b10 * dist);

					//energy += exp_part - _patch_E_cut;

					//angular part of interaction

					number cosa1 = pp->patches[pi].a1 * r_dist_dir;
					number cosb1 = -qq->patches[pj].a1 * r_dist_dir;
					number cosa2b2 = pp->patches[pi].a2 * qq->patches[pj].a2;

					number  ta1 = LRACOS(cosa1);
					number  tb1 = LRACOS(cosb1);
					number  ta2b2 = LRACOS(cosa2b2);

					number  fa1 =  _V_mod(PLPATCH_VM1,ta1);
					number  fb1 =  _V_mod(PLPATCH_VM1,tb1) ;

					number  fa2b2 =   _V_mod(PLPATCH_VM3,ta2b2);

					number f1 =  K * (exp_part - _patch_E_cut);
					number angular_part =  fa1 * fb1 * fa2b2;

					energy_ij = f1 * angular_part;
					energy += energy_ij;

					//PRO LUKASE:
					// if (energy_ij < -1.)
					// {
					//	printf("@@@@: particle: %d , patch %d binds to particle %d, patch %d \n",p->index,pi,q->index,pj);
					// }

					//patchy locking enabled for MD
					if(update_forces && this->_no_multipatch)
					{
						if (energy_ij < this->_lock_cutoff )
						{
							qq->set_lock(pj, p->index,pi,energy_ij);
							pp->set_lock(pi, q->index,pj,energy_ij);
						}
						else
						{
							qq->unlock(pj);
							pp->unlock(pi);

						}

					}

					//printf("Patches %d and %d distance %f , K:%f, attraction ene: %f, exp_part: %f, E_cut: %f, angular ene: %f, cos: %f %f %f\n",pp->patches[pi].id,qq->patches[pj].id,dist,K,(exp_part - _patch_E_cut),exp_part,_patch_E_cut,angular_part,pp->patches[pi].a1.x,pp->patches[pi].a1.y,pp->patches[pi].a1.z);

					//printf("Patches %d and %d distance %f , K:%f, attraction ene: %f, exp_part: %f, E_cut: %f, angular ene: %f\n",pp->patches[pi].id,qq->patches[pj].id,dist,K,(exp_part - _patch_E_cut),exp_part,_patch_E_cut,angular_part);
					if(update_forces ) {
						number f1D =  (5 * exp_part * r8b10);
						LR_vector tmp_force = patch_dist * (f1D * angular_part);
						//printf("CRITICAL 1 Adding %f %f %f \n",tmp_force.x,tmp_force.y,tmp_force.z);

						number fa1Dsin =  _V_modDsin(PLPATCH_VM1,ta1);
						number fb1Dsin =  _V_modDsin(PLPATCH_VM1,tb1);
						/*
					tmp_force += (pp->patches[pi].a1 -  r_dist_dir * cosa1) * (f1 * fa1Dsin *  fb1* fa2b2 / rdist);
					tmp_force += -(qq->patches[pj].a1 +  r_dist_dir * cosb1) * (f1 * fa1 *  fb1Dsin * fa2b2 / rdist);
						 */
						//printf("CRITICAL 2 Adding %f %f %f \n",tmp_force.x,tmp_force.y,tmp_force.z);
						//torque VM3
						number fa2b2Dsin =  _V_modDsin(PLPATCH_VM3,ta2b2);
						LR_vector dir = -pp->patches[pi].a2.cross(qq->patches[pj].a2) *  (f1 * fa1 * fb1 * fa2b2Dsin );
						LR_vector torqueq = dir;
						LR_vector torquep = dir;


						//torque VM1
						dir = r_dist_dir.cross(pp->patches[pi].a1);
						torquep += dir * (f1 * fa1Dsin * fb1 );

						//torque VM2
						dir = r_dist_dir.cross(qq->patches[pj].a1);
						torqueq += dir * (f1 * fa1 * fb1Dsin );


						torquep += ppatch.cross(tmp_force);
						torqueq += qpatch.cross(tmp_force);

						tmp_force += (pp->patches[pi].a1 -  r_dist_dir * cosa1) * (f1 * fa1Dsin *  fb1* fa2b2 / rdist);
						tmp_force += -(qq->patches[pj].a1 +  r_dist_dir * cosb1) * (f1 * fa1 *  fb1Dsin * fa2b2 / rdist);
						/*
					tmp_force += (pp->patches[pi].a1 -  r_dist_dir * cosa1) * (f1 * fa1Dsin *  fb1 / rdist);
					tmp_force += -(qq->patches[pj].a1 +  r_dist_dir * cosb1) * (f1 * fa1 *  fb1Dsin  / rdist);
						 */

						p->torque -= p->orientationT * torquep;
						q->torque += q->orientationT * torqueq;

						p->force -= tmp_force;
						q->force += tmp_force;

						// we are in the single bond per patch condition and hence we can safely return here
						//return energy;
					}
				}
			}
		}
	}
	return energy;
}

AllostericPatchyInteraction::AllostericPatchyInteraction() : BaseInteraction()  {
	//this->_int_map[PATCHY] = &AllostericPatchyInteraction::_patchy_interaction;
	_close_vertexes = NULL;

	_no_multipatch = 0;

	this->_int_map[EXCVOL] = &AllostericPatchyInteraction::_exc_vol_interaction;
	this->_int_map[PATCHY] = &AllostericPatchyInteraction::_patchy_interaction;

	//this->_int_map[PATCHY] = &AllostericPatchyInteraction::_patchy_LJ_interaction;
	//this->_int_map[PATCHY] = &AllostericPatchyInteraction::_patchy_CUT_LJ4896_noEXC_interaction;

	//default option, very wide!
	PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.;
	PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.7;
	PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 3.10559;
	PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 0.46;
	PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 0.133855;

	//all these variables will be initialized later from input file:
	_N_patch_types = 0;
	_N_particle_types = 0;
	N_patches = 0;
	_E_cut = 0;
	_patch_E_cut = 0;
	_patch_alpha = 0;
	_patch_pow_alpha = 0;

	_patch_types    = 0;
	_particle_types = 0;
	_narrow_type = 0;
	_use_torsion = true;
	_same_type_bonding = true;

	_interaction_tensor = false;

	_interaction_patch_types = 0; //a 2d matrix of patches (each patch has a unique id) that checks if patches can interact; is only used if interaction tensor is true
	_interaction_table_types = 0;

}

// deconstructor. deallocates memory arrays

AllostericPatchyInteraction::~AllostericPatchyInteraction() {

	delete [] _patch_types;
	delete [] _particle_types;
	delete [] _interaction_table_types;
	delete [] _interaction_patch_types;
	delete [] _close_vertexes;
}

//takes an input string from a file (patches.txt or equivelant) and returns a patch object

AllostericPatch AllostericPatchyInteraction::_process_patch_type(std::string input_string)
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
    AllostericPatch loaded_patch(a1,a2,position,id,color,strength, true, allostery_conditional);

	//printf("Loaded patch %d with color %d \n",loaded_patch.id,loaded_patch.color);
	// deallocate memory for file reader
	delete obs_input;
	//return patch object
	return loaded_patch;
}

// loads particle type data from a file (particles.txt or equivelant) and returns a PatchyShapeParticle object

AllostericPatchyParticle AllostericPatchyInteraction::_process_particle_type(std::string input_string)
{
	input_file *obs_input = Utils::get_input_file_from_string(input_string);
	int type;
	getInputInt(obs_input,"type",&type,1);

	std::vector<AllostericPatch > all_patches;
	int _N_patches;
	std::string patches;
	if( getInputString(obs_input,"patches",patches,1) == KEY_FOUND )
	{
		//now process a list of patches: 1,2,3,4
		std::replace( patches.begin(), patches.end(), ',', ' ');
		std::stringstream s(patches);
		int patch_id;
		while( s >> patch_id)
		{
            AllostericPatch patch(this->_patch_types[patch_id]);
			all_patches.push_back(patch);
			OX_LOG(Logger::LOG_INFO,"Particle of type %d adding a patch of color %d",type,patch.color);
			//s >> patch_id;
		}
	}

	_N_patches = all_patches.size();
	int N_vertexes = 0;
	if (this->_shape == ICOSAHEDRON_SHAPE)
	{
		N_vertexes = 12;
	}
	OX_LOG(Logger::LOG_INFO,"Particle of type %d has %d vertexes",type,N_vertexes);
    AllostericPatchyParticle p;
	switch (this->_allostery_logic_type){
	case ALLOSTERY_LOGIC_NONE:
		p = AllostericPatchyParticle(_N_patches,type,N_vertexes);
		break;
	case ALLOSTERY_LOGIC_SIMPLE:
		p = SimpleAllosteryPatchyShapeParticle(_N_patches,type,N_vertexes);
		break;
	case ALLOSTERY_LOGIC_COMPLEX:
		p = AdvAllosteryPatchyShapeParticle(_N_patches, type, N_vertexes);
		break;
	}
	if (N_vertexes  != 0)
	{
		p._set_vertexes();
	}
	int position = 0;
	for(typename std::vector<AllostericPatch >::iterator i = all_patches.begin(); i != all_patches.end(); ++i)
	{
		p.add_patch(*i,position);
		position++;
	}
	p.init_allostery();

//	ParticleStateChange test_change(default_state, 2, 0);

//	std::vector<int> t = (*p.allostery_map)[test_change];
	return p;
}


void
//loads patch data from patches.txt and particles.txt
AllostericPatchyInteraction::_load_patchy_particle_files(std::string& patchy_file, std::string& particle_file)
{
	//first process patches
	FILE *fpatch = fopen(patchy_file.c_str(),"r");
	if(!fpatch) // if file not found
	{
		throw oxDNAException("Could not open file %s ",patchy_file.c_str());
	}
	input_file obs_input;
	loadInput(&obs_input,fpatch); // open patch file input

	int no = 0;
	char patch_no[1024]; //instantiate C string buffer - 1024 chars should be enough
	snprintf(patch_no,1020,"patch_%d",no); // write patch index
	std::string patch_string;

	while(  getInputString(&obs_input,patch_no,patch_string,0) == KEY_FOUND )
	{
		AllostericPatch patch = _process_patch_type(patch_string);
		if(no >= _N_patch_types)
		{
			throw oxDNAException("Number of patch types is larger than N_patch_types = %d ",_N_patch_types);
		}
		_patch_types[no] = patch;
		no++;
		snprintf(patch_no,1020,"patch_%d",no);
	}

	fclose(fpatch);

	//now process particles
	FILE *fparticle = fopen(particle_file.c_str(), "r");
	if(!fparticle)
	{
		throw oxDNAException("Could not open file %s ",particle_file.c_str());
	}
	input_file p_input;
	loadInput(&p_input, fparticle);

	int p_no = 0;
	char particle_no[1024];
	snprintf(particle_no, 1020, "particle_%d", p_no);
	std::string particle_string;

	while (getInputString(&p_input, particle_no, particle_string, 0) == KEY_FOUND) {
		switch (this->_allostery_logic_type){
		case ALLOSTERY_LOGIC_NONE:
			break;
		case ALLOSTERY_LOGIC_SIMPLE:
			break;
		case ALLOSTERY_LOGIC_COMPLEX:
            break;
		}
        AllostericPatchyParticle particle = _process_particle_type(particle_string);
		if(p_no >= _N_particle_types)
			throw oxDNAException ("More particle types in particle config file than specified in the input file. Aborting");

		_particle_types[p_no].copy_from(particle);
		p_no++;
		snprintf(particle_no, 1020, "particle_%d", p_no);
	}

	fclose(fparticle);

	OX_LOG(Logger::LOG_INFO, "Loaded %d patch types and %d particle types", no,p_no);
	if(p_no != _N_particle_types || no != _N_patch_types)
		throw oxDNAException ("More (or less) particle or patches types in particle/patchy config file than specified in the input file. Aborting");
}


void AllostericPatchyInteraction::_load_interaction_tensor(std::string &tensor_file)
{
	int p1, p2, patch1, patch2;
	std::ifstream inf(tensor_file.c_str());
	if (!inf.good())
		throw oxDNAException("Cannot open %s",tensor_file.c_str());
	char line[1024];
	inf.getline(line,1023);
	int c = 0;
	while(inf.good()) {

		std::istringstream  s(line);
		if(s.str().size() >= 4)
		{
			s >> p1 >> p2 >> patch1 >> patch2;
			if(p1 >= this->_N_particle_types || p2 >=  this->_N_particle_types || patch1 >= this->_particle_types[p1].N_patches || patch2 >= this->_particle_types[p2].N_patches  )
				throw oxDNAException("Invalid index found in  %s, line: %s",tensor_file.c_str(),line);

			int i1 = this->_particle_types[p1].patches[patch1].id;
			int i2 = this->_particle_types[p2].patches[patch2].id;
			c++;
			this->_interaction_patch_types[i1*_N_patch_types + i2] = this->_interaction_patch_types[i2*_N_patch_types +  i1] = 1;
			printf("DBG: Loaded %d %d %d %d, which is id: %d %d\n",p1,p2,patch1,patch2,i1,i2);
		}

		if(s.fail())
		{
			throw  oxDNAException ("Malformed line in interaction tensor file: %s ",line);
		}
		inf.getline(line,1023);
	}

	OX_LOG(Logger::LOG_INFO, "Loaded %d patch interaction tensor entries from %s", c,tensor_file.c_str());
}



// loads settings from the input file passed as an arguement
void AllostericPatchyInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);

	//getInputInt(&inp, "PATCHY_N", &_N_patches, 1);

	int narrow_type = 0;
	getInputInt(&inp,"narrow_type",&narrow_type,0);
	this->_narrow_type = narrow_type;

	//DETECT SHAPE!
	this->_shape = 0;
	std::string shapestring; //this file contains information about types of patches
	getInputString(&inp, "shape", shapestring, 1);
	if ( shapestring == string("sphere") )
	{
		this->_shape = SPHERE_SHAPE;
	}
	else if ( shapestring == string("icosahedron") )
	{
		this->_shape = ICOSAHEDRON_SHAPE;
	}
	else
	{
		throw oxDNAException("Unsupported shape: %s",shapestring.c_str());
	}

	OX_LOG(Logger::LOG_INFO, "Using excluded volume shape %d (%s)",this->_shape,shapestring.c_str());

	int no_multi = 0;
	if( getInputBoolAsInt(&inp,"no_multipatch",&no_multi,0) == KEY_FOUND)
	{
		this->_no_multipatch = (bool)no_multi;
	}
	else
	{
		this->_no_multipatch = true;
	}
	OX_LOG(Logger::LOG_INFO, "Using no_multipatch option: %d; only makes sense if used with MC2 MCMovePatchyShape!",this->_no_multipatch);

	/*
	int use_torsion = 1;
	if( getInputBoolAsInt(&inp,"use_torsion",&use_torsion,0) == KEY_FOUND)
	{
		this->_use_torsion = (bool)use_torsion;
	}
	else
	{
	   this->_use_torsion = true;
	}
	if(this->_use_torsion)
	{
	    this->_int_map[PATCHY] = &AllostericPatchyInteraction::_patchy_1PONLY_LJ4896_noEXC_interaction;
	    printf("Torsional constraints are on\n");
	}
	else
	{
		this->_int_map[PATCHY]  = &AllostericPatchyInteraction::_patchy_1PONLY_LJ4896_noEXCnoTorsion_interaction;
		printf("Torsional constraints are off\n");
	}
	 */

	//NOT IMPLEMENTED
	int same_type_bonding = 1;
	if( getInputBoolAsInt(&inp,"same_type_bonding",&same_type_bonding,0) == KEY_FOUND)
	{
		this->_same_type_bonding = (bool)same_type_bonding;
	}
	else
	{
		this->_same_type_bonding = true;
	}

	if(this->_same_type_bonding)
	{
		OX_LOG(Logger::LOG_INFO, "Particles of the same type can bond");
	}
	else
	{
		OX_LOG(Logger::LOG_INFO, "Particles of the same type cannot bond");
	}

	//load allosteric logic type
	std::string allostery_logic_type;
	if (getInputString(&inp, "allosteric_logic_type", allostery_logic_type, 0) == KEY_FOUND){
		if (allostery_logic_type == "simple"){
			this->_allostery_logic_type = ALLOSTERY_LOGIC_SIMPLE;
		} else if (allostery_logic_type == "complex") {
			this->_allostery_logic_type = ALLOSTERY_LOGIC_COMPLEX;
		}
		else {
			this->_allostery_logic_type = ALLOSTERY_LOGIC_NONE;
		}
	}
	else {
		this->_allostery_logic_type = ALLOSTERY_LOGIC_NONE;
	}

	getInputInt(&inp,"patch_types_N",&_N_patch_types,1);
	getInputInt(&inp,"particle_types_N",&_N_particle_types,1);

	_patch_types = new AllostericPatch [_N_patch_types];
	_particle_types = new AllostericPatchyParticle [_N_particle_types];

	std::string patchy_file; //this file contains information about types of patches
	getInputString(&inp, "patchy_file", patchy_file, 1);
	std::string particle_file; //this file contains information about types of particles
	getInputString(&inp, "particle_file", particle_file, 1);

	_load_patchy_particle_files(patchy_file,particle_file);

	int interaction_tensor = 0;
	this->_interaction_patch_types = new int [_N_patch_types * _N_patch_types]();

	if( getInputBoolAsInt(&inp,"interaction_tensor",&interaction_tensor,0) == KEY_FOUND)
	{
		this->_interaction_tensor = (bool)interaction_tensor;
	}
	if(this->_interaction_tensor)  //possible interactions are specified by an external file
	{
		std::string interaction_tensor_file; //this file contains information about types of patches
		getInputString(&inp, "interaction_tensor_file", interaction_tensor_file, 1);

		this->_load_interaction_tensor(interaction_tensor_file);
	}
	else  //no external file provided; we assign allowed interactions based on colors of patches!
	{
		for(int i = 0; i < _N_particle_types; i++)
		{
            AllostericPatchyParticle  *p = &_particle_types[i];

			for(int j = i; j < _N_particle_types; j++)
			{
                AllostericPatchyParticle  *q = &_particle_types[j];
				for(int pi = 0; pi < p->N_patches; pi++)
				{
					for(int qj = 0; qj < q->N_patches; qj++)
					{
						int pid = p->patches[pi].id;
						int qid = q->patches[qj].id;
						if(p->patches[pi].color == q->patches[qj].color)
						{
							this->_interaction_patch_types[pid * _N_patch_types + qid] = this->_interaction_patch_types[qid * _N_patch_types + pid] = 1;
						}

						if(i == j && !this->_same_type_bonding)
						{
							this->_interaction_patch_types[pid * _N_patch_types + qid] = this->_interaction_patch_types[qid * _N_patch_types + pid] = 0;
						}
					}
				}

			}
		}
	}

	_interaction_table_types =  new int [_N_particle_types * _N_particle_types]();
	for(int i = 0; i < _N_particle_types; i++)
	{
		for(int j = i; j < _N_particle_types; j++)
		{
			bool possible_bond = false;
            AllostericPatchyParticle  *p = &_particle_types[i];
            AllostericPatchyParticle  *q = &_particle_types[j];
			for(int pi = 0; pi < p->N_patches; pi++)
			{
				for(int qi = 0; qi < q->N_patches; qi++)
				{
					if(_bonding_allowed(p,q,pi,qi))
					{
						possible_bond = true;
						break;
					}

				}
				if(possible_bond)
					break;
			}
			if(possible_bond)
			{
				_interaction_table_types[i*_N_particle_types + j]  = _interaction_table_types[j*_N_particle_types + i] = 1;
				// printf("Particle types %d and %d can interact!\n",i,j);
			}
		}
	}

	//now actually check it for the particles themselves:

	/*
	printf("CRITICAL: Loaded %d particle types:\n",_N_particle_types);
	for(int i = 0; i < _N_particle_types; i++)
	{
	    	printf("CRITICAL: Particle %d has %d patches, which have the following a1,a2: ",i,_particle_types[i].N_int_centers);
	    	for(int c = 0; c < _particle_types[i].N_int_centers; c++)
	    	{
	    		//printf("%d ",dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].color);
	    		printf("%f %f %f xx %f %f %f ",(_particle_types + i)->patches[c].a1_x,(_particle_types + i)->patches[c].a1_y,(_particle_types + i)->patches[c].a1_z,
	    				(_particle_types + i)->patches[c].a2_x,(_particle_types + i)->patches[c].a2_y,(_particle_types + i)->patches[c].a2_z
	    		);
	    	}
	    	printf("\n");
	 }
	 */
	//if(getInputInt(&inp, "PATCHY_N_B", &_N_patches_B, 0) == KEY_FOUND) _is_binary = true;

	float tmp = 1.2;
	getInputFloat(&inp, "PATCHY_rcut", &tmp, 0);
	this->_rcut = (number) tmp;

	tmp = 0.12;
	getInputFloat(&inp, "PATCHY_alpha", &tmp, 0);
	_patch_alpha = (number) tmp;

	tmp = 0.5;
	getInputFloat(&inp, "PATCHY_radius", &tmp, 0);
	_sphere_radius = (number) tmp;

	tmp = -0.1;
	getInputFloat(&inp, "PATCHY_multi_cutoff", &tmp, 0);
	_lock_cutoff = (number) tmp;

	OX_LOG(Logger::LOG_INFO, "(AllostericPatchyInteraction) using radius=%g, alpha=%g, cutoff=%g, multipatch=%d", _sphere_radius, _patch_alpha, _lock_cutoff, _no_multipatch);

}


void AllostericPatchyInteraction::init() {

	_E_cut = powf((number) this->_rcut, -PATCHY_POWER);
	this->_sqr_rcut = SQR(this->_rcut);

	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);

	number r8b10 = powf(PATCHY_CUTOFF, (number) 8.f) / _patch_pow_alpha;
	_patch_E_cut = -1.001f * expf(-(number)0.5f * r8b10 * SQR(PATCHY_CUTOFF));
	OX_LOG(Logger::LOG_INFO, "INFO: setting _patch_E_cut to %f which is %f, with alpha=%f", _patch_E_cut,-1.001f * expf(-(number)0.5f * r8b10 * SQR(PATCHY_CUTOFF)),_patch_alpha);

	// choose narrow type
	// narrow type determines the distance required for binding between two patches(?)
	if(_narrow_type == 0)
	{
		//default option, very wide!
		PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.;
		PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.7;
		PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 3.10559;
		PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 0.46;
		PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 0.133855;
	}
	else if(_narrow_type == 1)
	{
		//narrower!
		PLPATCHY_THETA_T0[0] =  PLPATCHY_THETA_T0[1] = 0.;
		PLPATCHY_THETA_TS[0] =  PLPATCHY_THETA_TS[1] = 0.2555;
		PLPATCHY_THETA_TC[0] =  PLPATCHY_THETA_TC[1] = 1.304631441617743;
		PLPATCHY_THETA_A[0]  =  PLPATCHY_THETA_A[1] = 3.;
		PLPATCHY_THETA_B[0]  =  PLPATCHY_THETA_B[1] = 0.7306043547966398;
	}
	else if(_narrow_type == 2)
	{
		//narrower!
		PLPATCHY_THETA_T0[0] =  PLPATCHY_THETA_T0[1] = 0.;
		PLPATCHY_THETA_TS[0] =  PLPATCHY_THETA_TS[1] = 0.2555;
		PLPATCHY_THETA_TC[0] =  PLPATCHY_THETA_TC[1] = 0.782779;
		PLPATCHY_THETA_A[0]  =  PLPATCHY_THETA_A[1] = 5.;
		PLPATCHY_THETA_B[0]  =  PLPATCHY_THETA_B[1] = 2.42282;
	}
	else if (_narrow_type == 3)
	{
		//narrower:
		PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.;
		PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.17555;
		PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 0.4381832920710734;
		PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 13.;
		PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 8.68949241736805;
	}
	else if(_narrow_type == 4)
	{
		//narrowest:
		PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.f;
		PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.17555;
		PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 0.322741;
		PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 17.65;
		PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 21.0506;
	}
	else{
		throw oxDNAException ("Invalid narrow_type option, has to be between 0 and 4");
	}

	if (_shape == ICOSAHEDRON_SHAPE) {
		_init_icosahedron();
	}
}

// allocates memory for particles, which must be then copies from AllostericPatchyInteraction::_particle_types

void AllostericPatchyInteraction::allocate_particles(BaseParticle **particles, int N) {
	for(int i = 0; i < N; i++) {
		//int i_patches = (i < _N_A) ? _N_patches : _N_patches_B;
		particles[i] = new AllostericPatchyParticle(1);
	}
}



number AllostericPatchyInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}



number AllostericPatchyInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return (number) 0.f;
}



number AllostericPatchyInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	LR_vector computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	number energy;

	energy = this->_patchy_interaction(p, q, r, update_forces);
	energy += this->_exc_vol_interaction(p,q,r,update_forces);

	return energy;
}



/*

void AllostericPatchyInteraction::generate_random_configuration(BaseParticle **particles, int N, number box_side) {
	number old_rcut = this->_rcut;
	this->_rcut = 1;

	this->_create_cells(particles, N, box_side, true);

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		bool inserted = false;
		int cell_index;
		do {
			p->pos = LR_vector(drand48()*box_side, drand48()*box_side, drand48()*box_side);
			cell_index = (int) ((p->pos.x / box_side - floor(p->pos.x / box_side)) * (1.f - FLT_EPSILON) * this->_cells_N_side);
			cell_index += this->_cells_N_side * ((int) ((p->pos.y / box_side - floor(p->pos.y / box_side)) * (1.f - FLT_EPSILON) * this->_cells_N_side));
			cell_index += this->_cells_N_side * this->_cells_N_side * ((int) ((p->pos.z / box_side - floor(p->pos.z / box_side)) * (1.f - FLT_EPSILON) * this->_cells_N_side));

			inserted = true;
			for(int c = 0; c < 27; c ++) {
				int j = this->_cells_head[this->_cells_neigh[cell_index][c]];
				while (j != P_INVALID) {
					BaseParticle *q = particles[j];
					if(p->pos.minimum_image(q->pos, box_side).norm() < SQR(this->_rcut)) inserted = false;
					j = this->_cells_next[q->index];
				}
			}
		} while(!inserted);

		int old_head = this->_cells_head[cell_index];
		this->_cells_head[cell_index] = i;
		this->_cells_index[i] = cell_index;
		this->_cells_next[i] = old_head;

		p->orientation.v1 = Utils::get_random_vector();
		p->orientation.v2 = Utils::get_random_vector();
		p->orientation.v3 = Utils::get_random_vector();
		Utils::orthonormalize_matrix(p->orientation);
	}

	this->_rcut = old_rcut;
	this->_delete_cell_neighs();

}
 */

// loads topology data from a topology file

void AllostericPatchyInteraction::read_topology(int N, int *N_strands, BaseParticle **particles) {
	*N_strands = N;
	int N_types;
	std::ifstream topology(this->_topology_filename, ios::in); // open a file stream to topology file
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[4096]; // bad code, has caused segfault issue in the past. TODO: rewrite
	topology.getline(line, 512);
	sscanf(line, "%*d %d\n", &N_types);
	allocate_particles(particles, N);
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
		if(total_count >= N || type >= _N_particle_types || type < 0)
			throw oxDNAException("The sum of number of species is larger than number of particles, or unknown type encountered. Aborting while processing file %s (%d %d %d)", this->_topology_filename, total_count, _N_particle_types, type);
		int i = total_count;
		particles[i]->copy_from(this->_particle_types[type]);

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
	OX_LOG(Logger::LOG_INFO, "There were %d particles, %d types, and finished allocation, and line was %s, and N_particle types was %d",N,N_types,line,_N_particle_types);
	int patch_index = 0;
	for(int i = 0; i < N; i++)
	{
		//printf("Particle %d has %d patches, which have the following colors: ",i,dynamic_cast<PatchyShapeParticle *>(particles[i])->N_patches);
		particles[i]->set_positions();
		for(int c = 0; c < dynamic_cast<AllostericPatchyParticle *>( particles[i])->N_patches; c++)
		{
			//now assign each patch its unique id
			//printf("%d ",dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].color);
			dynamic_cast<AllostericPatchyParticle *>(particles[i])->patches[c].index = patch_index;
			patch_index++;
			//printf("%d (%f %f %f) ",dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].color, dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].a1.x,dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].a1.y,dynamic_cast<PatchyShapeParticle *>(particles[i])->patches[c].a1.z);
		}
		//printf("\n");
	}


	this->N_patches = patch_index;
}




void AllostericPatchyInteraction::check_input_sanity(BaseParticle **particles, int N) {

}




number AllostericPatchyInteraction::just_two_patch_interaction(AllostericPatchyParticle *p, AllostericPatchyParticle *q, int pi,int  qi,LR_vector *r)
{
	number rnorm = r->norm();
	if(rnorm > this->_sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

    AllostericPatchyParticle *pp = static_cast<AllostericPatchyParticle *>(p);
    AllostericPatchyParticle *qq = static_cast<AllostericPatchyParticle *>(q);

	LR_vector ppatch = p->int_centers[pi];
	if(this->_bonding_allowed(pp,qq,pi,qi) )
	{
		number K = pp->patches[pi].strength;
		LR_vector qpatch = q->int_centers[qi];

		LR_vector patch_dist = *r + qpatch - ppatch;
		number dist = patch_dist.norm();
		LR_vector patch_dist_dir = patch_dist / sqrt(dist);
		number rdist = sqrt(rnorm);
		LR_vector r_dist_dir = *r / rdist;

		if(dist < SQR(PATCHY_CUTOFF)) {
			//distance part of attractive interaction
			number r8b10 = dist*dist*dist*dist / this->_patch_pow_alpha;
			number exp_part = -1.001f * exp(-(number)0.5f * r8b10 * dist);

			//angular part of interaction

			number cosa1 = pp->patches[pi].a1 * r_dist_dir;
			number cosb1 = -qq->patches[qi].a1 * r_dist_dir;
			number cosa2b2 = pp->patches[pi].a2 * qq->patches[qi].a2;

			number  ta1 = LRACOS(cosa1);
			number  tb1 = LRACOS(cosb1);
			number  ta2b2 = LRACOS(cosa2b2);

			number  fa1 =  this->_V_mod(PLPATCH_VM1,ta1);
			number  fb1 =  this->_V_mod(PLPATCH_VM1,tb1) ;
			number  fa2b2 =   this->_V_mod(PLPATCH_VM3,ta2b2);

			number f1 =  K * (exp_part - this->_patch_E_cut);

			return f1;
			/*
				    number angular_part;

				    if(this->_use_torsion)
				    	angular_part =  fa1 * fb1* fa2b2;
				    else
				    	angular_part = fa1 * fb1;

				    energy = f1 * angular_part;
			 */
			//energy += energy_ij;
		}
	}

	return energy;
}





void AllostericPatchyInteraction::_init_icosahedron(void)
{
	_close_vertexes = new int[12 * 5]; // 5 close vertexes for each vertex
	_tworinscribed = (number) (2.*sqrt((1./12.) + (1./(6.*sqrt(5.))))); // twice the radius of inscribed sphere
	// compute the close vertexes thanks to a temporary hard icosahedron
    AllostericPatchyParticle tmp (0,0,12);
	tmp._set_icosahedron_vertexes();

	/*
	for(int i = 0; i < 12; i++)
	{
		printf("Icosahedron vertex %d = %f %f %f \n",i,tmp._vertexes[i][0],tmp._vertexes[i][1],tmp._vertexes[i][2]);
	}
	 */

	tmp.pos = LR_vector ((number)0.f, (number)0.f, (number)0.f);
	tmp.orientation = LR_matrix ((number)1.f, (number)0.f, (number)0.f,
			(number)0.f, (number)1.f, (number)0.f,
			(number)0.f, (number)0.f, (number)1.f);
	tmp.orientationT = tmp.orientation.get_transpose();
	tmp.set_positions();

	for (int i = 0; i < 12; i ++) {
		int myn = 0;
		for (int j = 0; j < 12; j ++) {
			if (i != j and ((tmp.int_centers[i])*(tmp.int_centers[j]) > 0.)) {
				_close_vertexes[5*i + myn] = j;
				myn ++;
			}
		}
		if (myn != 5) throw oxDNAException("something wrong while initializing hardico interaction...");
	}

	// now we need to order the arrays so that we move counterclockwise
	for (int i = 0; i < 12; i ++) { // for each vertex
		bool good = false;
		while (good == false) {
			good = true;
			for (int j = 0; j < 4; j ++) {
				if (tmp.int_centers[_close_vertexes[5*i+j]]*tmp.int_centers[_close_vertexes[5*i + (j + 1)]] < 0.) {
					int save = _close_vertexes[5*i + j + 1];
					for (int k=1; k < 5-j-1; k ++) {
						_close_vertexes[5*i + j + k] = _close_vertexes[5*i + j + k + 1];
					}
					_close_vertexes[5*i+4] = save;
					good = false;
					break;
				}
			}
		}
	}
}



void AllostericPatchyInteraction::_init_patchy_locks(shared_ptr<ConfigInfo> Info)
{
	if (Info == NULL)
	{
		Info = ConfigInfo::instance();
	}

	//printf("!!INITPATCHY: Starting locking, we have %d particles\n",*Info->N);
	Info->lists->global_update();
	for(int pid = 0; pid < Info->N(); pid++)
	{
        AllostericPatchyParticle *p = static_cast< AllostericPatchyParticle *>(Info->particles()[pid]);
		//printf("XXXXXXXXXXXXXXXXXXxParticle has %d partches\n",p->N_patches);
		//fflush(stdout);
		//p->_set_vertexes();
		p->unlock_patches();
	}

	for(int pid = 0; pid < Info->N(); pid++)
	{

        AllostericPatchyParticle *p = static_cast< AllostericPatchyParticle *>(Info->particles()[pid]);
		//printf("Pidf is %d, p->index is %d\n",pid,p->index);

		//std::vector<BaseParticle *> neighs = Info->lists->get_neigh_list(p);;
		for(int n = pid+1; n < Info->N(); n++) {
            AllostericPatchyParticle *qq =   static_cast< AllostericPatchyParticle *>(Info->particles()[n]);
			//printf("qid is %d, qq->index is %d\n",n,qq->index);

			LR_vector r = Info->box->min_image(p,qq);
			for(int ppatch = 0; ppatch < p->N_patches; ppatch++)
			{
				for(int qqpatch = 0; qqpatch < qq->N_patches; qqpatch++)
				{

					number new_ene = this->just_two_patch_interaction(p,qq,ppatch,qqpatch,&r);

					if(new_ene < this->get_patch_cutoff_energy())
					{
						//throw oxDNAException("Locking ");
						p->set_lock(ppatch, qq->index,qqpatch);
						qq->set_lock(qqpatch, p->index,ppatch);
						//printf("!!INITPATCHY: Locking %d (%d) to %d (%d) \n",p->index,ppatch,qq->index,qqpatch);
					}
				}
			}
		}
	}

	/*
	for(int pid = 0; pid < *Info->N; pid++)
	{
		PatchyShapeParticle *p = dynamic_cast< PatchyShapeParticle *>(Info->particles[pid]);
		//p->unlock_patches();

		printf("Afterwards Particle %d: , patches: ",pid);
		for(int kk = 0; kk < p->N_patches; kk++)
		{
			printf(" %d-(%d %d) ",kk,p->patches[kk].locked_to_particle,p->patches[kk].locked_to_patch);
		}
		printf("\n");

	}
	 */
	//throw oxDNAException("Finished init");
}



void AllostericPatchyInteraction::check_patchy_locks(shared_ptr<ConfigInfo> Info)
{
	if (Info == NULL) {
		Info = ConfigInfo::instance();
	}

	//Info->lists->global_update();
	for(int pid = 0; pid < Info->N(); pid++)
	{
        AllostericPatchyParticle *p = static_cast< AllostericPatchyParticle *>(Info->particles()[pid]);
		//std::vector<BaseParticle *> neighs = Info->lists->get_neigh_list(p);
		for(int qid = 0; qid < pid; qid++) {
            AllostericPatchyParticle *qq =   static_cast< AllostericPatchyParticle *>(Info->particles()[qid]);

			LR_vector r = Info->box->min_image(p,qq);
			for(int ppatch = 0; ppatch < p->N_patches; ppatch++)
			{
				for(int qqpatch = 0; qqpatch < qq->N_patches; qqpatch++)
				{

					this->_no_multipatch = false;
					number new_ene = this->just_two_patch_interaction(p,qq,ppatch,qqpatch,&r);
					this->_no_multipatch = true;

					if(new_ene < this->get_patch_cutoff_energy())
					{
						// either the particles are locked to each other
						bool mutual_lock = p->patches[ppatch].locked_to(qid,qqpatch) && qq->patches[qqpatch].locked_to(pid,ppatch);
						// or either of them is locked to somone else
						bool external_lock = p->patches[ppatch].is_locked() || qq->patches[qqpatch].is_locked(); 

						//This can be either already locked:
						if (mutual_lock == false and external_lock == false) {
							printf("particles %d (patch %d) and %d (patch %d) have energy %g\n", pid, ppatch, qid, qqpatch, new_ene);
							int tmpi, tmpj;
							p->patches[ppatch].get_lock(tmpi, tmpj);
							printf("%d(%d) locked to %d(%d)\n", pid, ppatch, tmpi, tmpj);
							qq->patches[qqpatch].get_lock(tmpi, tmpj);
							printf("%d(%d) locked to %d(%d)\n", qid, qqpatch, tmpi, tmpj);
							//throw oxDNAException("Found a case where lock is missing: %d (%d) - %d (%d), %f ",pid,ppatch,qid,qqpatch, new_ene);
							printf("Found a case where lock is missing: %d (%d) - %d (%d), %f ",pid,ppatch,qid,qqpatch, new_ene);
						}
					}
					else //they should not be locked to each other!
					{
						if(p->patches[ppatch].locked_to(qid,qqpatch) || qq->patches[qqpatch].locked_to(pid,ppatch))
						{
							throw oxDNAException("Found a wrong lock, they should be not locked: %d (%d) - %d (%d), %f",pid,ppatch,qid,qqpatch,new_ene);
						}

					}
				}
			}
		}
	}
}

/**
 * Imported from a previous version of BaseInteraction because I don't understand what it does and I'm
 * afraid to touch it
 * @param that
 * @param name
 * @param p
 * @param q
 * @param r
 * @param update_forces
 * @return
 */
number AllostericPatchyInteraction::_pair_interaction_term_wrapper(BaseInteraction *that, int name, BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
    LR_vector computed_r;
    if(r == NULL) {
        computed_r = this->_box->min_image(p->pos, q->pos);
        r = &computed_r;
    }

    // _int_map.find(name) returns an iterator which points to the
    // requested function pointer. We dereference it, apply it to the
    // that, which is a pointer to an instance of a class which inherits
    // from BaseInteraction. We then pass in the required parameters and
    // return the interaction energy.
    typename interaction_map::iterator interaction = _int_map.find(name);
    if(interaction == _int_map.end()) {
        throw oxDNAException("%s, line %d: Interaction term '%d' not found", __FILE__, __LINE__, name);
    }

    return (interaction->second)(p, q, r, update_forces);
}