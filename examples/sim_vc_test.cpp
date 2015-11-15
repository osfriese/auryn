/*
 * Copyright (C) 2015 Tobias Neumann <tobias.neumann@scible.de>
 *
 * Based on example: sim_isp_orig.cpp
 *
 * this is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * this file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "auryn.h"


#define NE 4096 //64     // number of neurons
#define NI 1024 //16     // 4:1

#define FW 64      // fieldWidth
#define PB 32      // pinwheel box width

#define NP 100     // number of spiking inputneurons
#define NSTIM 20

using namespace std;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

int main(int ac, char* av[])
{
    // params for Connection
    // These values are for testing only and not functional with normal 
    // spiking groups. You will have to adjust them to get a good working network
    double wEE = 0.2;   // Weight for EE connection
    double wEI = 0.01;   // Weight for EI connection
    double wIE = 0.01;   // Weight for IE connection
    double wII = 0.01;   // Weight for II connection
    double wLR = 0.02;   // Weight for long-range connection
    double sigma = 32;   // p = exp(-disance/sigma)

    // params for SpikingGroup
    double poisson_rate = 10.;     // mean spinking rate /s???
    double sparseness_afferents = 0.05;

    bool quiet = false;
    double simtime = 1000. ;


    // handle command line options

    string infilename = "";
    string outputfile = "out";
    string netstatfile = "";
    string stimfile = "";
    string strbuf ;

    int errcode = 0;

    // Options

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
//            ("quiet", "quiet mode")
//            ("load", po::value<string>(), "input weight matrix")
            ("out", po::value<string>(), "output filename")
//            ("stimfile", po::value<string>(), "stimulus file")
//            ("eta", po::value<double>(), "learning rate")
//            ("kappa", po::value<double>(), "target rate")
            ("simtime", po::value<double>(), "simulation time")
//            ("active", po::value<bool>(), "toggle learning")
//            ("poisson", po::value<bool>(), "toggle poisson stimulus")
//            ("winh", po::value<double>(), "inhibitory weight multiplier")
//            ("wei", po::value<double>(), "ei weight multiplier")
//            ("chi", po::value<double>(), "chi current multiplier")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }
/*
        if (vm.count("quiet")) {
            quiet = true;
        }

        if (vm.count("load")) {
            cout << "input weight matrix "
                 << vm["load"].as<string>() << ".\n";
            infilename = vm["load"].as<string>();
        }
*/
        if (vm.count("out")) {
            cout << "output filename "
                 << vm["out"].as<string>() << ".\n";
            outputfile = vm["out"].as<string>();
        }
/*
        if (vm.count("stimfile")) {
            cout << "stimfile filename "
                 << vm["stimfile"].as<string>() << ".\n";
            stimfile = vm["stimfile"].as<string>();
        }

        if (vm.count("eta")) {
            cout << "eta set to "
                 << vm["eta"].as<double>() << ".\n";
            eta = vm["eta"].as<double>();
        }

        if (vm.count("kappa")) {
            cout << "kappa set to "
                 << vm["kappa"].as<double>() << ".\n";
            kappa = vm["kappa"].as<double>();
        }
*/
        if (vm.count("simtime")) {
            cout << "simtime set to "
                 << vm["simtime"].as<double>() << ".\n";
            simtime = vm["simtime"].as<double>();
        }
/*
        if (vm.count("active")) {
            cout << "stdp active : "
                 << vm["active"].as<bool>() << ".\n";
            stdp_active = vm["active"].as<bool>();
        }

        if (vm.count("poisson")) {
            cout << "poisson active : "
                 << vm["poisson"].as<bool>() << ".\n";
            poisson_stim = vm["poisson"].as<bool>();
        }


        if (vm.count("winh")) {
            cout << "inhib weight multiplier : "
                 << vm["winh"].as<double>() << ".\n";
            winh = vm["winh"].as<double>();
        }

        if (vm.count("wei")) {
            cout << "ei weight multiplier : "
                 << vm["wei"].as<double>() << ".\n";
            wei = vm["wei"].as<double>();
        }

        if (vm.count("chi")) {
            cout << "chi multiplier : "
                 << vm["chi"].as<double>() << ".\n";
            chi = vm["chi"].as<double>();
        }
*/
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }



    // BEGIN Global definitions
    mpi::environment env(ac, av);
    mpi::communicator world;
    communicator = &world;

    netstatfile = outputfile;
    stringstream oss;
    oss << outputfile << "." << world.rank();
    string basename = oss.str();
    oss << ".log";
    string logfile = oss.str();
    logger = new Logger(logfile,world.rank(),DEBUG);

    sys = new System(&world);
    // END Global definitions



    logger->msg("Setting up neuron groups ...",PROGRESS,true);

    // Create NeuronGroup
    AdExGroup * neurons_e = new AdExGroup(NE);
    AdExGroup * neurons_i = new AdExGroup(NI);

    // initial membrane potentials with a Gaussian: random_mem(mean,sigma)
    neurons_e->random_mem(-60e-3,5e-3);
    neurons_i->random_mem(-60e-3,5e-3);

    // Create SpikingGroup (Test Noise)
    PoissonGroup * poisson = new PoissonGroup(NP,poisson_rate);

    logger->msg("Setting up connections ...",PROGRESS,true);

    // Connect NeuronGroup intern: GeoConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight conWeight, AurynWeight lrcWeight,
    // bool isSourceInhibitory, bool isDestInhibitory, unsigned int fieldWidth, double sigma,
    // TransmitterType transmitter, string name)
    GeoConnection * con_ee = new GeoConnection(neurons_e, neurons_e, wEE, wLR, false,false,FW,sigma,GLUT, "EEConnection");
    GeoConnection * con_ei = new GeoConnection(neurons_e, neurons_i, wEI, wLR, false,true,FW,sigma,GLUT, "EIConnection");
    GeoConnection * con_ii = new GeoConnection(neurons_i, neurons_i, wII, wLR, true,true,FW,sigma,GABA, "IIConnection");
    GeoConnection * con_ie = new GeoConnection(neurons_i, neurons_e, wIE, wLR, true,false,FW,sigma,GABA, "IEConnection");

    // save connections for testing with matlab
    con_ee->write_to_file("log/con_ee.pat");
    con_ei->write_to_file("log/con_ei.pat");
    con_ii->write_to_file("log/con_ii.pat");
    con_ie->write_to_file("log/con_ie.pat");


    // Connect SpikingGroup to NeuronGroup
    SparseConnection * con_exte = new SparseConnection(poisson,neurons_e,0.3,sparseness_afferents,GLUT);
    SparseConnection * con_exti = new SparseConnection(poisson,neurons_i,0.3,sparseness_afferents,GLUT);

    // Create Monitors
    logger->msg("Setting up monitors ...",PROGRESS,true);

    // Spikes neurons_e -> %.e.ras
    strbuf = basename;
    strbuf += ".e.ras";
    SpikeMonitor * smon_e = new SpikeMonitor( neurons_e , strbuf.c_str() );

    // Spikes neurons_i -> %.i.ras
    strbuf = basename;
    strbuf += ".i.ras";
    SpikeMonitor * smon_i = new SpikeMonitor( neurons_i, strbuf.c_str() );
    smon_i->set_offset(NE);

    // membran voltage of neurons_e -> %.volt
//    strbuf = basename;
//    strbuf += ".volt";
//    StateMonitor * vmon = new StateMonitor( neurons_e, record_neuron, "mem", strbuf.c_str() );

    // AMPA of neurons_e -> %.ampa
//    strbuf = basename;
//    strbuf += ".ampa";
//    StateMonitor * amon = new StateMonitor( neurons_e, record_neuron, "g_ampa", strbuf.c_str() );

    // GABA of neurons_e -> %.gaba
//    strbuf = basename;
//    strbuf += ".gaba";
//    StateMonitor * gmon = new StateMonitor( neurons_e, record_neuron, "g_gaba", strbuf.c_str() );

    // tracks population firing rate and breaks a run if it goes out of bound
    // RateChecker(source, min, max, tau (time to compute average))
    RateChecker * chk = new RateChecker( neurons_e , 0.1 , 10000. , 100e-3);




    /*
    // TODO
    for ( int j = 0; j < NE ; j++ ) {
      neurons_e->set_bg_current(j,bg_current);
    }

    for ( int j = 0; j < NI ; j++ ) {
      neurons_i->set_bg_current(j,bg_current);
    }


    // extra stimulus
    if (!stimfile.empty()) {
        char ch;
        NeuronID counter = 0;
        ifstream fin(stimfile.c_str());
        while (!fin.eof() && counter<NE) {
            ch = fin.get();
            if (ch == '1') {
                if (poisson_stim==true) {
                    for (int i = 0 ; i < NP ; ++i)
                        con_exte->set(i,counter,w_ext);     // Connection SpikingGroup -> neurons_e
                } else {
                    neurons_e->set_bg_current(counter,chi*bg_current);
                }
            }
            counter++;
        }
        fin.close();
    }*/

    if (!infilename.empty()) {
        sys->load_network_state(infilename);
    }


    /*
     *      Start Simulation
     */

    logger->msg("Simulating ...",PROGRESS,true);
    //con_ie->stdp_active = stdp_active;

    sys->run(simtime);
    logger->msg("Saving network state ...",PROGRESS,true);
    sys->save_network_state(netstatfile);

    logger->msg("Freeing ...",PROGRESS,true);
    delete sys;

    logger->msg("Exiting ...",PROGRESS,true);
    return errcode;
}
