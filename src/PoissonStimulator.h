/* 
* Copyright 2014-2015 Friedemann Zenke
*
* This file is part of Auryn, a simulation package for plastic
* spiking neural networks.
* 
* Auryn is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* Auryn is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with Auryn.  If not, see <http://www.gnu.org/licenses/>.
*
* If you are using Auryn or parts of it for your work please cite:
* Zenke, F. and Gerstner, W., 2014. Limits to high-speed simulations 
* of spiking neural networks using general-purpose computers. 
* Front Neuroinform 8, 76. doi: 10.3389/fninf.2014.00076
*/

#ifndef POISSONSTIMULATOR_H_
#define POISSONSTIMULATOR_H_

#include "auryn_definitions.h"
#include "System.h"
#include "Logger.h"
#include "Monitor.h"
#include "NeuronGroup.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/poisson_distribution.hpp>

using namespace std;

/*! \brief Stimulator class to inject timeseries of currents to patterns (subpopulations) of neurons 
 * 
 * Instances of this class inject currents that vary over time to subpopulations of the NeuronGroup assigned.
 */

class PoissonStimulator : protected Monitor
{
private:

	static boost::mt19937 gen; 
	boost::poisson_distribution<int> * dist;
	boost::variate_generator<boost::mt19937&, boost::poisson_distribution<int> > * die;


	/*! Vector storing all the current values */
	AurynState * currents;

	/*! Vector storing all new current values */
	AurynState * newcurrents;

	/*! Target membrane */
	auryn_vector_float * target_vector;

	/*! Scale stimulus rate */
	AurynFloat poisson_rate;

	/*! Scale stimulus size */
	AurynFloat poisson_weight;

	/*! Default init method */
	void init(NeuronGroup * target, AurynFloat rate, AurynWeight w );

	void free();

	/*! Returns the lambda parameter of the pmf for Poisson. */
	AurynFloat get_lambda();

protected:

	/*! The target NeuronGroup */
	NeuronGroup * dst;

	
public:
	/*! Default Constructor 
	 * @param[target] The target spiking group. 
	 * @param[rate]   The firing rate of each the Poisson process.
	 * @param[weight] The weight or unit of amount of change on the state variable
	 */
	PoissonStimulator(NeuronGroup * target, AurynFloat rate=100.0, AurynWeight w = 0.1 );

	/*! Default Destructor */
	virtual ~PoissonStimulator();


	/*! Sets the event rate of the underlying Poisson generator */
	void set_rate(AurynFloat rate);

	/*! Returns the  event rate of the underlying Poisson generator */
	AurynFloat get_rate();

	/*! Seeds the random number generator of all PoissonStimulator objects */
	void seed(int s);


	/*! Sets the state that is stimulated with Poisson input.
	 * This must be a valid state vector name (default = mem) */
	void set_target_state( string state_name = "mem" );

	/*! Implementation of necessary propagate() function. */
	void propagate();

};

#endif /*POISSONSTIMULATOR_H_*/
