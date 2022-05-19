/*
 * PDClusterTopology.h
 *
 *  Created on: May 5, 2022
 *      Author: josh evans
 *  Initially a copy of PLClusterTopology.h, with modifications to output more complex cluster data
 *  PD stands for "Patchy Detailed"
 */

#ifndef PDCLUSTERTOPOLOGY_H_
#define PDCLUSTERTOPOLOGY_H_

#define MCMOVE_CUSTOM
#include "../../../../src/Observables/BaseObservable.h"

/**
 * @brief prints out each cluster of bonded (locked) particles
 *
 * @verbatim
particle1_id = <int> (particle 1 id)
particle2_id = <int> (particle 2 id)
@endverbatim
 */
template<typename number>
class PDClusterTopology: public BaseObservable<number> {
protected:

    bool _show_types ;
public:
    PDClusterTopology();
	virtual ~PDClusterTopology();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable<float> * make_PDClusterTopology_float() { return new PDClusterTopology<float>(); }
extern "C" BaseObservable<double> * make_PDClusteropology_double() { return new PDClusterTopology<double>(); }


#endif /* PDCLUSTERTOPOLOGY_H_ */
