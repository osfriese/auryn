#ifndef GEOCONNECTION_H
#define GEOCONNECTION_H

#define PI 3.14159265


#include "SparseConnection.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <math.h>
using namespace std;

/*! \brief This class implements geometric short and long-range connections
 *
 * x/y  1 2 3 4 5 6 7 8 9 ...
 *    1 x   x   x   x   x   x
 *    2   o   o   o   o   o
 *    3 x   x   x   x   x   x
 *    4   o   o   o   o   o
 *    5 x   x   x   x   x   x
 *
 * x: exibitory
 * o: inhibitory
 *
 * p = e^(-(x1-x2)^2 + (y1-y2)^2) / sigma )
 *
 *
 *
 */

class GeoConnection : public SparseConnection
{
private:
//    bool sourceIn;  // Are source and/or dest inhibitory Cells?
//    bool destIn;

    // params for grid
    double sourceGap;
    unsigned int sourceWidth;

    double destGap;
    unsigned int destWidth;

    // Size of Pinwheel boxes (width)
    unsigned int oBoxSize;
    unsigned short orientationGap;

    // Size of Cortex Area: fieldSize x fieldSize
//    unsigned int fieldSize;

    // Weights
    AurynWeight cWeight; // normal connection weight
    AurynWeight oWeight; // Weight for long range connection

    // Probability for long-range connections and
    double oProb;
    double sigma;   // for exp function

    void virtual_serialize(boost::archive::binary_oarchive & ar, const unsigned int version );
    void virtual_serialize(boost::archive::binary_iarchive & ar, const unsigned int version );

public:

    GeoConnection(SpikingGroup * source, NeuronGroup * destination, TransmitterType transmitter=GLUT, string name="GeoConnection");

    GeoConnection(SpikingGroup * source, NeuronGroup * destination,
                  AurynWeight conWeight, AurynWeight lrcWeight,
                  bool isSourceInhibitory, bool isDestInhibitory,
                  unsigned int fieldWidth, double sigma,
                  TransmitterType transmitter=GLUT, string name="GeoConnection");

    void set_weights(AurynWeight normalConnectionWeight, AurynWeight longrangeConnectionWeight);

    void set_params(unsigned int fieldWidth, unsigned int pinweheelSize,
                    unsigned short orientationGap,
                    bool isSourceInhibitory, bool isDestInhibitory);

    void set_params(unsigned int sourceWidth, unsigned int destWidth, unsigned int pinweheelSize,
                    unsigned short orientationGap,
                    double sourceGap, double destGap);

    void set_probability(double sigma, double longrangeProbability);

protected:

    /*! Calculates probability for connection based on geometry */
    double getProbability(NeuronID i, NeuronID j, double sigma);

    /*! Calculates x and y position of neurons */
    void get_xy(NeuronID id, double & x, double &y, bool source);

    /*! returns orientation selectivity in degree for any x,y coordinate */
    unsigned short get_orientation(double x, double y);

    /*! true whenever two Neurons have same orientation selectivity */
    bool same_orientation(double xi, double yi, double xj, double yj);

    /*! true whenever two Neurons have same orientation selectivity */
    bool same_orientation(NeuronID i, NeuronID j);

    int connectLikeInVC_Block(NeuronID i, NeuronID j);

    /*! Wrapper for connectLikeInVC_Block for all Neurons */
    void connectLikeInVC(bool skip_diag=false);

    /*! Method wich generates connections
     */
    void connectLikeInVC_Block(NeuronID lo_row,
                            NeuronID hi_row,
                            NeuronID lo_col,
                            NeuronID hi_col,
                            bool skip_diag=false );

    void allocateMEM(SpikingGroup * source, NeuronGroup * destination);
};

#endif // GEOCONNECTION_H

