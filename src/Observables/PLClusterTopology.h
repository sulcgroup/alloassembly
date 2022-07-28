/*
 * HBEnergy.h
 *
 *  Created on: Feb 14, 2013
 *      Author: petr
 */

//Minimally modified by Josh to remove template args, etc.

#ifndef PLCLUSTERTOP_H_
#define PLCLUSTERTOP_H_

#define MCMOVE_CUSTOM
#include "BaseObservable.h"

/**
 * @brief Prints out all the interactions between all pairs of nucleotides with non-zero interaction energies (default) or between a pair of nucleotides (if specified)
 *
 * @verbatim
particle1_id = <int> (particle 1 id)
particle2_id = <int> (particle 2 id)
@endverbatim
 */
class PLClusterTopology: public BaseObservable {
protected:

    bool _show_types ;
public:
    PLClusterTopology();
    virtual ~PLClusterTopology();

    std::string get_output_string(llint curr_step);

    virtual void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* PAIRENERGY_H_ */
