module Astro
{

sequence<double> agileEvt;
const int agileEvtSize = 6;

sequence<double> agileLog;
const int agileLogSize = 13;

struct AgileEvtKey
{
double time;
double ra;
double dec;
double energy;
double theta;
double evstatus;
};

struct AgileLogKey
{
double time;
double livetime;
double logStatus;
double mode;
double phase;
};

sequence<int> SimpleSeq;
sequence<SimpleSeq> Matrix;

sequence<double> Ra;
sequence<double> Dec;

sequence<AgileEvtKey> SeqEvtKey;

interface AstroInterface
{

Matrix calculateMapKey(SeqEvtKey keys);
Matrix calculateMapVector(Ra raVector, Dec decVector);
void shutdown();

};

};