#include "GeoConnection.h"


void GeoConnection::virtual_serialize(boost::archive::binary_oarchive &ar, const unsigned int version)
{
    Connection::virtual_serialize(ar,version);
    ar & *w;
}

void GeoConnection::virtual_serialize(boost::archive::binary_iarchive &ar, const unsigned int version)
{
    Connection::virtual_serialize(ar,version);
    ar & *w;
}

GeoConnection::GeoConnection(SpikingGroup *source, NeuronGroup *destination, TransmitterType transmitter, string name)
    :SparseConnection(source,destination,transmitter)
{
    allocateMEM(source,destination);
}

GeoConnection::GeoConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight conWeight, AurynWeight lrcWeight,
                             bool isSourceInhibitory, bool isDestInhibitory, unsigned int fieldWidth, double sigma,
                             TransmitterType transmitter, string name)
    :SparseConnection(source,destination,transmitter)
{
    set_weights(conWeight,lrcWeight);
    set_params(fieldWidth,fieldWidth/2,5,isSourceInhibitory,isDestInhibitory);
    set_probability(sigma,0.4);
    allocateMEM(source,destination);
    connectLikeInVC((source==destination));
}



void GeoConnection::set_weights(AurynWeight normalConnectionWeight, AurynWeight longrangeConnectionWeight)
{
    this->cWeight = normalConnectionWeight;
    this->oWeight = longrangeConnectionWeight;
}

void GeoConnection::set_params(unsigned int fieldWidth, unsigned int pinweheelSize, unsigned short orientationGap, bool isSourceInhibitory, bool isDestInhibitory)
{
    this->destWidth = fieldWidth;
    this->sourceWidth = fieldWidth;
    this->orientationGap = orientationGap;
    this->sourceGap = (isSourceInhibitory)? 1.5: 1;
    this->destGap = (isDestInhibitory)? 1.5: 1;
    this->oBoxSize = pinweheelSize;
}

void GeoConnection::set_params(unsigned int sourceWidth, unsigned int destWidth, unsigned int pinweheelSize, unsigned short orientationGap, double sourceGap, double destGap)
{
    this->destWidth = destWidth;
    this->sourceWidth = sourceWidth;
    this->orientationGap = orientationGap;
    this->sourceGap = sourceGap;
    this->destGap = destGap;
    this->oBoxSize = pinweheelSize;
}

void GeoConnection::set_probability(double sigma, double longrangeProbability)
{
    this->sigma = sigma;
    this->oProb = longrangeProbability;
}

double GeoConnection::getProbability(NeuronID i, NeuronID j, double sigma)
{
    double xi,yi,xj,yj;
    get_xy(i,xi,yi,true);
    get_xy(j,xj,yj,false);

    if (same_orientation(xi,yi,xj,yj)){
        return this->oProb;
    }

    double x = pow(xj-xi,2);
    double y = pow(yj-yi,2);

    return exp(-1*(x+y) / sigma);
}

void GeoConnection::get_xy(NeuronID id, double &x, double &y, bool source)
{
    id += 1;    // NeuronID starts with 0
    if (source) {
        unsigned int r = id % sourceWidth;
        x = (r==0)? sourceWidth + sourceGap : r + sourceGap;
        y = floor((id-x)/sourceWidth)+1+sourceGap;
    }else {
        unsigned int r = id % destWidth;
        x = (r==0)? destWidth + destGap : r + destGap;
        y = floor((id-x)/destWidth)+1+destGap;
    }
}

unsigned short GeoConnection::get_orientation(double x, double y)
{
    double b = (double) oBoxSize;
    double m = (double) round(b/2);

    // calculate quadrant
    unsigned int fx = floor((x-1)/b);
    unsigned int fy = floor((y-1)/b);

    // Center Point of quadrant
    double mx = m + fx * b;
    double my = m + fy * b;

    // center has 0 degree
    if(mx==x && my == y) return 0;

    // determine rotation direction
    short clockwise = (fx^fy)? -1 : 1;

    // calculate angle
    double a = (atan2((my-y),(mx-x))*180) / PI -90;
    int o = clockwise * a + 0.5;

    // return value between 0-360 degree
    return (o < 0)? o + 360: o;
}

bool GeoConnection::same_orientation(double xi, double yi, double xj, double yj)
{
    unsigned short oi = get_orientation(xi,yi);
    unsigned short oj = get_orientation(xj,yj);

    if ((oj > oi-orientationGap) && (oj < oi + orientationGap)) {
        return true;
    }
    return false;
}

bool GeoConnection::same_orientation(NeuronID i, NeuronID j)
{
    double xi,yi,xj,yj;
    get_xy(i,xi,yi,true);
    get_xy(j,xj,yj,false);
    return same_orientation(xi,yi,xj,yj);
}

int GeoConnection::connectLikeInVC_Block(NeuronID i, NeuronID j)
{
    unsigned int r = 0;
    unsigned int s = 1;


}

void GeoConnection::connectLikeInVC(bool skip_diag)
{
    if ( dst->evolve_locally() ) { // if there are no local units there is no need for synapses

        stringstream oss;
        oss << "GeoConnection: ("<< get_name() <<"): Randomfill with weight "<< cWeight <<  " and long range weight " << oWeight
            << ". FieldSize is " << sourceWidth << "|" << destWidth << " and PinwheelSize is " << oBoxSize << ".\n";
        logger->msg(oss.str(),DEBUG);

        w->clear();
        connectLikeInVC_Block(0,get_m_rows(),0,get_n_cols(),skip_diag);
    }
    finalize();
}

void GeoConnection::connectLikeInVC_Block(NeuronID lo_row, NeuronID hi_row, NeuronID lo_col, NeuronID hi_col, bool skip_diag)
{
    int r = 0; // these variables are used to speed up building the matrix if the destination is distributed
    int s = 1;

    r = communicator->rank() - dst->get_locked_rank();
    s = dst->get_locked_range();

    // Create uniform int generator
    boost::random::uniform_int_distribution<> dist(1,100);
    boost::variate_generator<boost::mt19937&, boost::random::uniform_int_distribution<> > die(GeoConnection::sparse_connection_gen, dist);

    if (!has_been_allocated)
        throw AurynConnectionAllocationException();

    AurynLong idim = (hi_row - lo_row);
    AurynLong jdim = (hi_col - lo_col) / s;

    if ( (hi_col - lo_col) % s > r ) { // some ranks have one more "carry" neuron
        jdim += 1;
    }

    AurynLong x = 0;
    AurynLong stop = idim*jdim;
    AurynLong count = 0;
    NeuronID i = 0;
    NeuronID j = 0;
    AurynWeight temp_weight = 0;
    while ( x < stop ) {

        // calculate i and j position
        i = lo_row + x / jdim;
        j = lo_col + s*(x % jdim) + r;

        double p = getProbability(i,j,sigma)*100;
        double ps = die();

        if ( (j >= lo_col) && (!skip_diag || i!=j) && p > ps) {

            try {
                temp_weight = (same_orientation(i,j))? oWeight: cWeight;
                if ( push_back(i,j,temp_weight) )
                    count++;
            }
            catch ( AurynMatrixDimensionalityException )
            {
                stringstream oss;
                oss << "GeoConnection: ("
                    << get_name()
                    <<"): Trying to add elements outside of matrix (i="
                    << i
                    << "j="
                    << j
                    << ", "
                    << count
                    << "th element) ";
                logger->msg(oss.str(),ERROR);
                return;
            }
            catch ( AurynMatrixPushBackException )
            {
                stringstream oss;
                oss << "GeoConnection: ("<< get_name()
                    << "): Failed pushing back element. Maybe due to out of order pushing? "
                    << " (" << i << "," << j << ") "
                    << " with count=" << count
                    << " in connect_block_random ( fill_level= " << w->get_fill_level() << " )";
                logger->msg(oss.str(),ERROR);
                return;
            }
            catch ( AurynMatrixBufferException )
            {
                stringstream oss;
                oss << "GeoConnection: ("
                    << get_name()
                    <<"): Buffer full after pushing "
                    << count
                    << " elements."
                    << " There are pruned connections!";
                logger->msg(oss.str(),ERROR);
                return;
            }
        }

        x++;

    }

    stringstream oss;
    oss << "GeoConnection: ("<< get_name() <<"): Finished connectLikeInVC ["
        << lo_row << ", " << hi_row << ", " << lo_col << ", " << hi_col <<  "] " << " (stop count "
        << std::scientific << setprecision(4) << (double) stop
        << ") and successfully pushed " << (double) count <<  " entries. "
        << "Resulting overall sparseness " << 1.*get_nonzero()/src->get_pre_size()/dst->get_post_size()
        << " Size of w: " << get_nonzero();

    logger->msg(oss.str(),DEBUG);
}

void GeoConnection::allocateMEM(SpikingGroup *source, NeuronGroup *destination)
{
    AurynLong size = source->get_size() * destination->get_size();
    allocate_manually(size);
}



