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

#include "FileInputGroup.h"

void FileInputGroup::init(string filename)
{
	sys->register_spiking_group(this);

	ftime = 0;
	lastspike = 0;
	therewasalastspike = false;

	active = true;

	if ( evolve_locally() ) {
		spkfile.open(filename.c_str(),ifstream::in);
		if (!spkfile) {
		  cerr << "Can't open input file " << filename << endl;
		  exit(1);
		}
	}
}

FileInputGroup::FileInputGroup(NeuronID n, string filename) : SpikingGroup(n, 0.0 ) // last 0 enforces RankLock
{
	playinloop = false;
	dly = 0;
	off = 0;
	init(filename);
}

FileInputGroup::FileInputGroup(NeuronID n, string filename, 
		bool loop, AurynFloat delay) 
: SpikingGroup( n , 0.0 )
{
	playinloop = loop;
	dly = (AurynTime) (delay/dt);
	off = 0;
	init(filename);
}

FileInputGroup::~FileInputGroup()
{
	spkfile.close();
}


void FileInputGroup::evolve()
{

	if (active) {
		NeuronID i;
		AurynFloat t;

		if (ftime == sys->get_clock() && therewasalastspike) {
			if (localrank(lastspike))
				spikes->push_back(lastspike);
			therewasalastspike = false;
		}

		while (ftime <= sys->get_clock() && spkfile.getline(buffer, 256) ) {
			istringstream line ( buffer ) ;
			line >> t;
			ftime = t/dt+off;
			line >> i;
			if (ftime == sys->get_clock()) {
				if (localrank(lastspike)) 
					spikes->push_back(i);
			} else {
				lastspike = i;
				therewasalastspike = true;
			}
		}

		if ( playinloop && spkfile.eof() ) {
			off = ftime+dly;
			spkfile.clear();
			spkfile.seekg(0,ios::beg);
		}
	}
	else { // keep track of time
		off = sys->get_clock();
		ftime = off;
	}
}
