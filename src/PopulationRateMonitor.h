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

#ifndef POPULATIONRATEMONITOR_H_
#define POPULATIONRATEMONITOR_H_

#include "auryn_definitions.h"
#include "Monitor.h"
#include "System.h"
#include "SpikingGroup.h"
#include <fstream>
#include <iomanip>

using namespace std;

/*! \brief Monitor class to record population firing rates
 * 
 * Instances of this class record the population firing rate of the src SpikingGroup assigned.
 * Binning is done discretely in bins of size bsize that is directly transformed in discrete 
 * AurynTime steps. The default 
 */

class PopulationRateMonitor : protected Monitor
{
private:
	/*! Varible used to count the spike events of the src SpikingGroup */
	NeuronID counter;
	/*! Stepsize = binsize in units of AurynTime (dt) */
	AurynTime ssize;
	/*! Binsize used in seconds */
	AurynDouble invbsize;

protected:
	/*! The source SpikingGroup */
	SpikingGroup * src;
	/*! Default init method */
	void init(SpikingGroup * source, string filename, AurynDouble binsize);
	
public:
	/*! Default Constructor 
	 @param[source] The source spiking group.
	 @param[filename] The filename to write to (should be different for each rank.)
	 @param[binsize] The binsize used for counting in seconds.*/
	PopulationRateMonitor(SpikingGroup * source, string filename, AurynDouble binsize=1e-3);
	/*! Default Destructor */
	virtual ~PopulationRateMonitor();
	/*! Implementation of necessary propagate() function. */
	void propagate();
};

#endif /*POPULATIONRATEMONITOR_H_*/
