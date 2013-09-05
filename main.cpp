/**
 * The EQeq Method
 * Authors: Christopher E. Wilmer
 *          Randall Q. Snurr (advisor)
 *          Hansung Kim (car output)
 *          Patrick Fuller (streaming functionality)
 *          Louis Knapp (json output)
 */

#include <iostream>        // To read files
#include <fstream>         // To output files
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <map>             // For string enumeration (C++ specific)
#include <cmath>           // For basic math functions
#include <cstdlib>
using namespace std;

#define TABLE_OF_ELEMENTS_SIZE 84
#define PI 3.1415926535897932384626433832795    // 32 digits of PI

// This is a clumsy way to enable switch statements with atom labels in C++
enum StringAtomLabels {
    ev_H, ev_He, ev_Li, ev_Be, ev_B, ev_C, ev_N, ev_O, ev_F, ev_Ne, ev_Na,
    ev_Mg, ev_Al, ev_Si, ev_P, ev_S, ev_Cl, ev_Ar, ev_K, ev_Ca, ev_Sc, ev_Ti,
    ev_V , ev_Cr, ev_Mn, ev_Fe, ev_Co, ev_Ni, ev_Cu, ev_Zn, ev_Ga, ev_Ge,
    ev_As, ev_Se, ev_Br, ev_Kr, ev_Rb, ev_Sr, ev_Y , ev_Zr, ev_Nb, ev_Mo,
    ev_Tc, ev_Ru, ev_Rh, ev_Pd, ev_Ag, ev_Cd, ev_In, ev_Sn, ev_Sb, ev_Te,
    ev_I, ev_Xe, ev_Cs, ev_Ba, ev_La, ev_Ce, ev_Pr, ev_Nd, ev_Pm, ev_Sm,
    ev_Eu, ev_Gd, ev_Tb, ev_Dy, ev_Ho, ev_Er, ev_Tm, ev_Yb, ev_Lu, ev_Hf,
    ev_Ta, ev_W , ev_Re, ev_Os, ev_Ir, ev_Pt, ev_Au, ev_Hg, ev_Tl, ev_Pb,
    ev_Bi, ev_Po
};

// Map to associate the strings with the enum values
std::map<std::string, StringAtomLabels> s_mapStringAtomLabels;

class Coordinates {
    public:
        // Constructor
        Coordinates();

        double x;
        double y;
        double z;
};

class IonizationDatum {
    public:
        IonizationDatum();

        // TODO: Mass, radii and other properties can be added here if that would help for some reason
        string Label;
        vector<bool> isDataAvailable; // True or false
        vector<double> ionizationPotential; // The first 8 ionization potentials and the electron affinity
        int chargeCenter;
};

// EQeq function headers (alphabetical order)
void DetermineReciprocalLatticeVectors();
double GetJ(int i, int j);
void InitializeStringAtomLabelsEnumeration();
void LoadIonizationData(const char *filename);
void LoadChargeCenters(const char *filename);
void LoadCIFData(string data);
void LoadCIFFile(string filename); // Reads in CIF files, periodicity can be switched off
string OutputCIFData();
void OutputCIFFormatFile(string filename);
void OutputFile(string filename, string data);
string OutputPDBData();
void OutputPDBFormatFile(string filename);
string OutputMOLData();
void OutputMOLFormatFile(string filename); // Outputs 'RASPA' MOL file
string OutputCARData();
void OutputJSONFormatFile(string filename); // Outputs list of partial charges
string OutputJSONData();
void OutputCARFormatFile(string filename);
void Qeq();
void RoundCharges(int digits); // Make *slight* adjustments to the charges for nice round numbers

// Algebra helper functions (alphaAnglebetical order)
vector<double> Cross(vector<double> a, vector<double> b);
double Dot(vector<double> a, vector<double> b);
double Mag(vector<double> a);
double Round(double num);
vector<double> Scalar(double a, vector<double> b);
vector<double> SolveMatrix(vector<vector<double> > A, vector<double> b);

// Global variables
bool isPeriodic = true;
bool useEwardSums = true; // will use direct sums if false
double aLength; double bLength; double cLength;
double alphaAngle; double betaAngle; double gammaAngle;
double unitCellVolume;
vector<double> aV(3); vector<double> bV(3); vector<double> cV(3); // Real-space vectors
vector<double> hV(3); vector<double> jV(3); vector<double> kV(3); // Reciprocal-lattice vectors
int numAtoms; // To be read from input file
double Qtot; // To be read in from file
vector<Coordinates> Pos; // Array of atom positions
vector<double> J; // Atom "hardness"
vector<double> X; // Atom electronegativity
vector<double> Q; // Partial atomic charge
vector<string> Label; // Atom labels (e.g., "C1" "C2" "ZnCation" "dummyAtom")
vector<string> Symbol; // Atom symbols (e.g., "C" "O" "Zn")
vector<IonizationDatum> IonizationData(TABLE_OF_ELEMENTS_SIZE);

// Parameters and constants
double k = 14.4; // Physical constant: the vacuum permittivity 1/(4pi*epsi) [units of Angstroms * electron volts]
double eta = 50; // Ewald splitting parameter
double lambda = 1.2; // Coulomb scaling parameter
float hI0 = -2.0; // Default value used in paper
float hI1 = 13.598; // This is the empirically mesaured 1st ionization energy of hydrogen
int chargePrecision = 3; // Number of digits to use for point charges
int mR = 2;  int mK = 2;
int aVnum = mR; int bVnum = mR; int cVnum = mR; // Number of unit cells to consider in per. calc. ("real space")
int hVnum = mK; int jVnum = mK; int kVnum = mK; // Number of unit cells to consider in per. calc. ("frequency space")
string method = "ewald";

// Function to compile down to C. Used by the Python wrapper.
extern "C" {
char *run(const char *data, const char *outputType, double _lambda, float _hI0,
          int _chargePrecision, const char *_method, int _mR, int _mK,
          double _eta, const char *ionizationDataFilename,
          const char *chargeCentersFilename);
}

/*****************************************************************************/
/*****************************************************************************/

int main (int argc, char *argv[]) {
    if (argc <= 1) { cerr << "Error, invalid input!" << endl; exit(1); }
    if (argc > 2) lambda = atof(argv[2]); // The dielectric screening parameter (optional, default value above)
    if (argc > 3) hI0 = atof(argv[3]); // The electron affinity of hydrogen (optional, default value above)
    if (argc > 4) chargePrecision = atoi(argv[4]); // Num of digits to use for charges (optional, default value above)
    if (argc > 5) method = argv[5];
    if (argc > 6) mR = atoi(argv[6]);
    if (argc > 7) mK = atoi(argv[7]);
    if (argc > 8) eta = atof(argv[8]);

    // The only mandatory parameter is the input file/stream parameter
    run(argv[1], (const char *)"files", lambda, hI0, chargePrecision, method.c_str(),
        mR, mK, eta, (const char *)"ionizationdata.dat", (const char *)"chargecenters.dat");
    return 0;
}
/*****************************************************************************/
char *run(const char *data, const char *outputType, double _lambda, float _hI0,
                int _chargePrecision, const char *_method, int _mR, int _mK,
                double _eta, const char *ionizationDataFilename,
                const char *chargeCentersFilename) {
    char *output;
    string input, type, inputFilename, outString;
    input.assign(data);
    type.assign(outputType);
    method.assign(_method);

    // Converts to lowercase
    transform(type.begin(), type.end(), type.begin(), ::tolower);
    transform(method.begin(), method.end(), method.begin(), ::tolower);

    // Set global variables (underscore name mangling is reverse of convention :-( )
    lambda = _lambda;
    hI0 = _hI0;
    chargePrecision = _chargePrecision;
    mR = _mR;
    mK = _mK;
    eta = _eta;
    if (method.compare("nonperiodic") == 0) {
        isPeriodic = false;
    } else {
        useEwardSums = (method.compare("ewald") == 0);
    }

    InitializeStringAtomLabelsEnumeration();  // Part of the clumsy way to enable string-based switch statements
    LoadIonizationData(ionizationDataFilename);
    LoadChargeCenters(chargeCentersFilename);

    // EQeq uses globals. The one below was causing errors when EQeq was run
    // more than once. This is a quick fix. Longer term, remove globals.
    Pos.clear();
    J.clear();
    X.clear();
    Label.clear();
    Symbol.clear();

    // Quick hack. If the string ends in ".cif", it's a file. Else, it's data.
    if (input.substr(input.length() - 4) == ".cif") {
        LoadCIFFile(input);
        inputFilename = input;
    } else {
        LoadCIFData(input);
        inputFilename = "streamed";
    }

    cerr << "==================================================" << endl;
    cerr << "===== Calculating charges... please wait. ========" << endl;
    Qtot = 0; // Can be non-zero for non-periodic structures
    Qeq();
    RoundCharges(chargePrecision);
    cerr << "===== ... done!                           ========" << endl;
    cerr << "==================================================" << endl;

    char buffer[50];
    sprintf(buffer,"_EQeq_%s_%4.2f_%4.2f", method.c_str(), lambda, hI0);

    // This is the standard behavior from the command line
    if (type.compare("files") == 0) {
        OutputCIFFormatFile(inputFilename + buffer + ".cif");
        OutputMOLFormatFile(inputFilename + buffer + ".mol");
        OutputPDBFormatFile(inputFilename + buffer + ".pdb");
        OutputCARFormatFile(inputFilename + buffer + ".car");
        OutputJSONFormatFile(inputFilename + buffer + ".json");
        return 0;
    // These options allow output streaming of string data
    } else if (type.compare("cif") == 0) {
        outString = OutputCIFData();
    } else if (type.compare("pdb") == 0) {
        outString = OutputPDBData();
    } else if (type.compare("mol") == 0) {
        outString = OutputMOLData();
    } else if (type.compare("car") == 0) {
        outString = OutputCARData();
    } else if (type.compare("json") == 0) {
        outString = OutputJSONData();
    } else {
        cerr << "Output type \"" << outputType << "\" not supported!" << endl;
        exit(1);
        return 0;
    }
    output = new char[outString.length() + 1];
    strcpy(output, outString.c_str());
    return output;
}
/*****************************************************************************/
/*****************************************************************************/
Coordinates::Coordinates() {
  x = 0; y = 0; z = 0; // default coordinates
}
/*****************************************************************************/
IonizationDatum::IonizationDatum() {
    isDataAvailable.resize(9,false);
    ionizationPotential.resize(9,0);
    chargeCenter = 0;
}
/*****************************************************************************/
void DetermineReciprocalLatticeVectors() {
    vector<double> crs;
    double pf; // pf => PreFactor

    crs = Cross(bV, cV);
    pf = 2*PI / Dot(aV, crs);
    hV[0] = pf * crs[0];
    hV[1] = pf * crs[1];
    hV[2] = pf * crs[2];

    crs = Cross(cV, aV);
    pf = 2*PI / Dot(bV, crs);
    jV[0] = pf * crs[0];
    jV[1] = pf * crs[1];
    jV[2] = pf * crs[2];

    crs = Cross(aV, bV);
    pf = 2*PI / Dot(cV, crs);
    kV[0] = pf * crs[0];
    kV[1] = pf * crs[1];
    kV[2] = pf * crs[2];
}
/*****************************************************************************/
void InitializeStringAtomLabelsEnumeration() {
    s_mapStringAtomLabels["H "] = ev_H;    // 1
    s_mapStringAtomLabels["He"] = ev_He;// 2
    s_mapStringAtomLabels["Li"] = ev_Li;// 3
    s_mapStringAtomLabels["Be"] = ev_Be;// 4
    s_mapStringAtomLabels["B "] = ev_B;    // 5
    s_mapStringAtomLabels["C "] = ev_C;    // 6
    s_mapStringAtomLabels["N "] = ev_N;    // 7
    s_mapStringAtomLabels["O "] = ev_O;    // 8
    s_mapStringAtomLabels["F "] = ev_F;    // 9
    s_mapStringAtomLabels["Ne"] = ev_Ne;//10
    s_mapStringAtomLabels["Na"] = ev_Na;//11
    s_mapStringAtomLabels["Mg"] = ev_Mg;//12
    s_mapStringAtomLabels["Al"] = ev_Al;//13
    s_mapStringAtomLabels["Si"] = ev_Si;//14
    s_mapStringAtomLabels["P "] = ev_P;    //15
    s_mapStringAtomLabels["S "] = ev_S;    //16
    s_mapStringAtomLabels["Cl"] = ev_Cl;//17
    s_mapStringAtomLabels["Ar"] = ev_Ar;//18
    s_mapStringAtomLabels["K "] = ev_K ;//19
    s_mapStringAtomLabels["Ca"] = ev_Ca;//20
    s_mapStringAtomLabels["Sc"] = ev_Sc;//21
    s_mapStringAtomLabels["Ti"] = ev_Ti;//22
    s_mapStringAtomLabels["V "] = ev_V ;//23
    s_mapStringAtomLabels["Cr"] = ev_Cr;//24
    s_mapStringAtomLabels["Mn"] = ev_Mn;//25
    s_mapStringAtomLabels["Fe"] = ev_Fe;//26
    s_mapStringAtomLabels["Co"] = ev_Co;//27
    s_mapStringAtomLabels["Ni"] = ev_Ni;//28
    s_mapStringAtomLabels["Cu"] = ev_Cu;//29
    s_mapStringAtomLabels["Zn"] = ev_Zn;//30
    s_mapStringAtomLabels["Ga"] = ev_Ga;//31
    s_mapStringAtomLabels["Ge"] = ev_Ge;//32
    s_mapStringAtomLabels["As"] = ev_As;//33
    s_mapStringAtomLabels["Se"] = ev_Se;//34
    s_mapStringAtomLabels["Br"] = ev_Br;//35
    s_mapStringAtomLabels["Kr"] = ev_Kr;//36
    s_mapStringAtomLabels["Rb"] = ev_Rb;//37
    s_mapStringAtomLabels["Sr"] = ev_Sr;//38
    s_mapStringAtomLabels["Y "] = ev_Y ;//39
    s_mapStringAtomLabels["Zr"] = ev_Zr;//40
    s_mapStringAtomLabels["Nb"] = ev_Nb;//41
    s_mapStringAtomLabels["Mo"] = ev_Mo;//42
    s_mapStringAtomLabels["Tc"] = ev_Tc;//43
    s_mapStringAtomLabels["Ru"] = ev_Ru;//44
    s_mapStringAtomLabels["Rh"] = ev_Rh;//45
    s_mapStringAtomLabels["Pd"] = ev_Pd;//46
    s_mapStringAtomLabels["Ag"] = ev_Ag;//47
    s_mapStringAtomLabels["Cd"] = ev_Cd;//48
    s_mapStringAtomLabels["In"] = ev_In;//49
    s_mapStringAtomLabels["Sn"] = ev_Sn;//50
    s_mapStringAtomLabels["Sb"] = ev_Sb;//51
    s_mapStringAtomLabels["Te"] = ev_Te;//52
    s_mapStringAtomLabels["I "] = ev_I ;//53
    s_mapStringAtomLabels["Xe"] = ev_Xe;//54
    s_mapStringAtomLabels["Cs"] = ev_Cs;//55
    s_mapStringAtomLabels["Ba"] = ev_Ba;//56
    s_mapStringAtomLabels["La"] = ev_La;//57
    s_mapStringAtomLabels["Ce"] = ev_Ce;//58
    s_mapStringAtomLabels["Pr"] = ev_Pr;//59
    s_mapStringAtomLabels["Nd"] = ev_Nd;//60
    s_mapStringAtomLabels["Pm"] = ev_Pm;//61
    s_mapStringAtomLabels["Sm"] = ev_Sm;//62
    s_mapStringAtomLabels["Eu"] = ev_Eu;//63
    s_mapStringAtomLabels["Gd"] = ev_Gd;//64
    s_mapStringAtomLabels["Tb"] = ev_Tb;//65
    s_mapStringAtomLabels["Dy"] = ev_Dy;//66
    s_mapStringAtomLabels["Ho"] = ev_Ho;//67
    s_mapStringAtomLabels["Er"] = ev_Er;//68
    s_mapStringAtomLabels["Tm"] = ev_Tm;//69
    s_mapStringAtomLabels["Yb"] = ev_Yb;//70
    s_mapStringAtomLabels["Lu"] = ev_Lu;//71
    s_mapStringAtomLabels["Hf"] = ev_Hf;//72
    s_mapStringAtomLabels["Ta"] = ev_Ta;//73
    s_mapStringAtomLabels["W "] = ev_W ;//74
    s_mapStringAtomLabels["Re"] = ev_Re;//75
    s_mapStringAtomLabels["Os"] = ev_Os;//76
    s_mapStringAtomLabels["Ir"] = ev_Ir;//77
    s_mapStringAtomLabels["Pt"] = ev_Pt;//78
    s_mapStringAtomLabels["Au"] = ev_Au;//79
    s_mapStringAtomLabels["Hg"] = ev_Hg;//80
    s_mapStringAtomLabels["Tl"] = ev_Tl;//81
    s_mapStringAtomLabels["Pb"] = ev_Pb;//82
    s_mapStringAtomLabels["Bi"] = ev_Bi;//83
    s_mapStringAtomLabels["Po"] = ev_Po;//84
}
/*****************************************************************************/
double GetJ(int i, int j) {
    // Note to reader - significant consolidation of code may be possible in this function
    if (isPeriodic == false) {
        //////////////////////////////////////////////////////////////////////
        //  NonPeriodic                                                     //
        //////////////////////////////////////////////////////////////////////
        if (i == j) {
            return J[i]; // Return the hardness/idempotential
        } else {
            double dx = Pos[i].x - Pos[j].x;
            double dy = Pos[i].y - Pos[j].y;
            double dz = Pos[i].z - Pos[j].z;
            double RabSq = dx*dx + dy*dy + dz*dz;
            double Rab = sqrt(RabSq);

            double Jij = sqrt(J[i] * J[j]);
            double a = Jij / k;
            double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab); // Other functional forms are OK too

            double Jab = lambda * (k/2) * ((1/Rab) + orbitalOverlapTerm);

            return Jab;
        }
        //////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////
    } else
    if (isPeriodic == true) {
        aVnum = mR; bVnum = mR; cVnum = mR; // Number of unit cells to consider in per. calc. (in "real space")
        hVnum = mK; jVnum = mK; kVnum = mK; // Number of unit cells to consider in per. calc. (in "frequency space")
        if (useEwardSums == false) {
            //////////////////////////////////////////////////////////////////////
            // Direct sums                                                      //
            //////////////////////////////////////////////////////////////////////
            if (i == j) {

                double sigmaStar = 0;
                for (int u = -aVnum; u <= aVnum; u++) {
                    for (int v = -bVnum; v <= bVnum; v++) {
                        for (int w = -cVnum; w <= cVnum; w++) {
                            if (!((u==0)&&(v==0)&&(w==0))) {
                                double dx = u*aV[0] + v*bV[0] + w*cV[0];
                                double dy = u*aV[1] + v*bV[1] + w*cV[1];
                                double dz = u*aV[2] + v*bV[2] + w*cV[2];
                                double Rab = sqrt(dx*dx + dy*dy + dz*dz);
                                double RabSq = dx*dx + dy*dy + dz*dz;

                                double Jij = sqrt(J[i] * J[j]);
                                double a = Jij / k;
                                double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab);
                                // Other functional forms (for orbital overlap) are OK too

                                sigmaStar += (1/Rab) + orbitalOverlapTerm;
                            }
                        }
                    }
                }
                return J[i] + lambda * (k/2)*sigmaStar;

            } else {

                double sigma = 0;
                for (int u = -aVnum; u <= aVnum; u++) {
                    for (int v = -bVnum; v <= bVnum; v++) {
                        for (int w = -cVnum; w <= cVnum; w++) {
                            double dx = Pos[i].x - Pos[j].x + u*aV[0] + v*bV[0] + w*cV[0];
                            double dy = Pos[i].y - Pos[j].y + u*aV[1] + v*bV[1] + w*cV[1];
                            double dz = Pos[i].z - Pos[j].z + u*aV[2] + v*bV[2] + w*cV[2];
                            double Rab = sqrt(dx*dx + dy*dy + dz*dz);
                            double RabSq = dx*dx + dy*dy + dz*dz;

                            double Jij = sqrt(J[i] * J[j]);
                            double a = Jij / k;
                            double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab);
                            // Other functional forms (for orbital overlap) are OK too

                            sigma += (1/Rab) + orbitalOverlapTerm;
                        }
                    }
                }

                return lambda * (k/2) * sigma;
            }
        } else {
            //////////////////////////////////////////////////////////////////////
            // Ewald sums                                                       //
            //////////////////////////////////////////////////////////////////////
            if (i == j) {
                // Orbital energy term
                double orbital = 0;
                for (int u = -aVnum; u <= aVnum; u++) {
                    for (int v = -bVnum; v <= bVnum; v++) {
                        for (int w = -cVnum; w <= cVnum; w++) {
                            if ((u==0) && (v==0) && (w==0)) {
                                // do nothing
                            } else {
                                double dx = u*aV[0] + v*bV[0] + w*cV[0];
                                double dy = u*aV[1] + v*bV[1] + w*cV[1];
                                double dz = u*aV[2] + v*bV[2] + w*cV[2];
                                double Rab = sqrt(dx*dx + dy*dy + dz*dz);
                                double RabSq = dx*dx + dy*dy + dz*dz;

                                double Jij = sqrt(J[i] * J[j]);
                                double a = Jij / k;
                                double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab);
                                // Other functional forms (for orbital overlap) are OK too

                                orbital += orbitalOverlapTerm;
                            }
                        }
                    }
                }

                // Real-space Coulomb component
                double alphaStar = 0;
                for (int u = -aVnum; u <= aVnum; u++) {
                    for (int v = -bVnum; v <= bVnum; v++) {
                        for (int w = -cVnum; w <= cVnum; w++) {
                            if ((u==0) && (v==0) && (w==0)) {
                                // do nothing
                            } else {
                                double dx = u*aV[0] + v*bV[0] + w*cV[0];
                                double dy = u*aV[1] + v*bV[1] + w*cV[1];
                                double dz = u*aV[2] + v*bV[2] + w*cV[2];
                                double Rab = sqrt(dx*dx + dy*dy + dz*dz);

                                alphaStar += erfc( Rab / eta ) / Rab;
                            }
                        }
                    }
                }

                // K-space component
                double betaStar = 0;
                double h = 0; double b = 0;
                vector<double> RLV(3); // reciprocal lattice vector
                for (int u = -hVnum; u <= hVnum; u++) {
                    for (int v = -jVnum; v <= jVnum; v++) {
                        for (int w = -kVnum; w <= kVnum; w++) {
                            if ((u==0) && (v==0) && (w==0)) {
                                // do nothing
                            } else {
                                RLV[0] = u*hV[0] + v*jV[0] + w*kV[0];
                                RLV[1] = u*hV[1] + v*jV[1] + w*kV[1];
                                RLV[2] = u*hV[2] + v*jV[2] + w*kV[2];

                                h = Mag(RLV);
                                b = 0.5 * h * eta;

                                //beta += cos( RLV[0]*dx + RLV[1]*dy + RLV[2]*dz ) / (h*h) * exp(-b*b);
                                betaStar += 1 / (h*h) * exp(-b*b);
                            }
                        }
                    }
                }
                betaStar *= 4*PI / unitCellVolume;

                return J[i] + lambda * (k/2) * (alphaStar + betaStar + orbital - 2/(eta*sqrt(PI)));

            } else {
                // Orbital energy term
                double orbital = 0;
                for (int u = -aVnum; u <= aVnum; u++) {
                    for (int v = -bVnum; v <= bVnum; v++) {
                        for (int w = -cVnum; w <= cVnum; w++) {
                            double dx = Pos[i].x - Pos[j].x + u*aV[0] + v*bV[0] + w*cV[0];
                            double dy = Pos[i].y - Pos[j].y + u*aV[1] + v*bV[1] + w*cV[1];
                            double dz = Pos[i].z - Pos[j].z + u*aV[2] + v*bV[2] + w*cV[2];
                            double Rab = sqrt(dx*dx + dy*dy + dz*dz);
                            double RabSq = dx*dx + dy*dy + dz*dz;

                            double Jij = sqrt(J[i] * J[j]);
                            double a = Jij / k;
                            double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab);
                            // Other functional forms (for orbital overlap) are OK too

                            orbital += orbitalOverlapTerm;
                        }
                    }
                }

                // Real-space Coulomb component
                double alpha = 0;
                for (int u = -aVnum; u <= aVnum; u++) {
                    for (int v = -bVnum; v <= bVnum; v++) {
                        for (int w = -cVnum; w <= cVnum; w++) {
                            double dx = Pos[i].x - Pos[j].x + u*aV[0] + v*bV[0] + w*cV[0];
                            double dy = Pos[i].y - Pos[j].y + u*aV[1] + v*bV[1] + w*cV[1];
                            double dz = Pos[i].z - Pos[j].z + u*aV[2] + v*bV[2] + w*cV[2];
                            double Rab = sqrt(dx*dx + dy*dy + dz*dz);

                            alpha += erfc( Rab / eta ) / Rab;
                        }
                    }
                }

                // K-space component
                double beta = 0;
                double h = 0; double b = 0;
                vector<double> RLV(3); // reciprocal lattice vector
                for (int u = -hVnum; u <= hVnum; u++) {
                    for (int v = -jVnum; v <= jVnum; v++) {
                        for (int w = -kVnum; w <= kVnum; w++) {
                            if ((u==0) && (v==0) && (w==0)) {
                                // do nothing
                            } else {
                                RLV[0] = u*hV[0] + v*jV[0] + w*kV[0];
                                RLV[1] = u*hV[1] + v*jV[1] + w*kV[1];
                                RLV[2] = u*hV[2] + v*jV[2] + w*kV[2];

                                h = Mag(RLV);
                                b = 0.5 * h * eta;

                                double dx = Pos[i].x - Pos[j].x;
                                double dy = Pos[i].y - Pos[j].y;
                                double dz = Pos[i].z - Pos[j].z;

                                beta += cos( RLV[0]*dx + RLV[1]*dy + RLV[2]*dz ) / (h*h) * exp(-b*b);
                            }
                        }
                    }
                }
                beta *= 4*PI / unitCellVolume;

                return lambda * (k/2) * (alpha + beta + orbital);
            }
        }
    } else {
        cerr << "Serious error specifying periodic boundary conditions. Exiting" << endl;
        exit(1);
    }
}
/*****************************************************************************/
void LoadChargeCenters(const char *filename) {
    // Loads charge centers to be used, atoms are assumed to be

    ifstream fileInput(filename,ios::in);
    string tmp, tStr;
    int sInd, Z;

    if(!fileInput) { // Error checking
        printf("%s is not a valid filename\n\n", filename);
        exit(1);
    }

    while(!fileInput.eof()) {
        getline(fileInput, tmp); // Read line-by-line

        // Read atom symbol
        sInd = tmp.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", 0);

        // Handles blank lines and such (they pop up on UNIX machines)
        if (sInd == string::npos) {
            continue;
        }
        tStr = tmp.substr(sInd, 2);
        if (tStr[1] == '\t') tStr[1] = ' '; // Convert tabs to spaces
        Z = s_mapStringAtomLabels[tStr]; // Get Z number from label

        // Read charge center (must be a positive integer)
        sInd = tmp.find_first_of("0123456789", sInd);
        tStr = tmp.substr(sInd, 1);
        IonizationData[Z].chargeCenter = atoi(tStr.c_str());
    }

}
/*****************************************************************************/
void LoadIonizationData(const char *filename) {
    // Loads ionization data into a global vector called IonizationData

    ifstream fileInput(filename, ios::in);
    string data, tmp, cStr, tStr;
    int sInd = 0;
    int eInd = 0;

    if(!fileInput) { // Error checking
        printf("%s is not a valid filename\n\n", filename);
        exit(1);
    }

    while(!fileInput.eof()) { // Read file into a gigantic string
        getline(fileInput, tmp);
        data += tmp;
        data += "\n";
    }

    for (int i = 0; i < TABLE_OF_ELEMENTS_SIZE; i++) {
        // The Atom Label
        sInd = data.find("\t",sInd) + 1;
        eInd = data.find("\t",sInd);
        cStr = data.substr(sInd, eInd - sInd);
        if (cStr.length() == 1) cStr += " ";
        IonizationData[i].Label = cStr;

        // The Data Status
        sInd = data.find("\t",sInd) + 1;
        eInd = data.find("\t",sInd);
        cStr = data.substr(sInd, eInd - sInd);

        // The Electron Affinity
        sInd = data.find("\t",sInd) + 1;
        eInd = data.find("\t",sInd);
        cStr = data.substr(sInd, eInd - sInd);
        if ((cStr.find("<0.5") >= 0) && (cStr.find("<0.5") < 100)) {
            IonizationData[i].isDataAvailable[0] = true;
            IonizationData[i].ionizationPotential[0] = 0.5;
        } else
        if ((cStr.find("na") >= 0) && (cStr.find("na") < 100)){
            IonizationData[i].isDataAvailable[0] = false;
            IonizationData[i].ionizationPotential[0] = 0.0;
        } else {
            IonizationData[i].isDataAvailable[0] = true;
            IonizationData[i].ionizationPotential[0] = atof( cStr.c_str() );
        }

        int trueCount = 0;
        for (int j  = 1; j < 9; j++) {
            // The J'th Ionization Potential
            sInd = data.find("\t",sInd) + 1;
            eInd = data.find("\t",sInd);
            cStr = data.substr(sInd, eInd - sInd);
            if ((cStr.find("na") >= 0) && (cStr.find("na") < 100)){
                IonizationData[i].isDataAvailable[j] = false;
                IonizationData[i].ionizationPotential[j] = 0.0;
            } else
            if ((cStr.find("np") >= 0) && (cStr.find("np") < 100)){
                IonizationData[i].isDataAvailable[j] = false;
                IonizationData[i].ionizationPotential[j] = 0.0;
            } else {
                IonizationData[i].isDataAvailable[j] = true;
                IonizationData[i].ionizationPotential[j] = atof( cStr.c_str() );
                trueCount++;
            }
        }

        sInd = data.find("\n",sInd) + 1; // Go to next line
    }
}
/*****************************************************************************/
void LoadCIFFile(string filename) {
    // Two string index variables used for generating substrings from larger strings

    ifstream fileInput(filename.c_str(),ios::in);
    string data, tmp;

    if(!fileInput) { // Error checking
        printf("%s is not a valid filename\n\n", filename.c_str());
        exit(1);
    }

    while(!fileInput.eof()) { // Read file into a gigantic string
        getline(fileInput, tmp);
        data += tmp;
        data += "\n";
    }
    LoadCIFData(data);
}
/*****************************************************************************/
void LoadCIFData(string data) {
    string cStr; // current string
    string tStr; // temp string
    int sInd = 0, eInd = 0, iInd = 0;

    // Read in unit cell dimensions
    sInd = data.find("_cell_length_a") + 15;
    eInd = data.find("\n", sInd);
    cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
    aLength = atof( cStr.c_str() );

    sInd = data.find("_cell_length_b") + 15;
    eInd = data.find("\n", sInd);
    cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
    bLength = atof( cStr.c_str() );

    sInd = data.find("_cell_length_c") + 15;
    eInd = data.find("\n", sInd);
    cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
    cLength = atof( cStr.c_str() );

    // Read in unit cell angles
    sInd = data.find("_cell_angle_alpha") + 18;
    eInd = data.find("\n", sInd);
    cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
    alphaAngle = atof( cStr.c_str() );

    sInd = data.find("_cell_angle_beta") + 17;
    eInd = data.find("\n", sInd);
    cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
    betaAngle = atof( cStr.c_str() );

    sInd = data.find("_cell_angle_gamma") + 18;
    eInd = data.find("\n", sInd);
    cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
    gammaAngle = atof( cStr.c_str() );

    // Convert to radians
    alphaAngle *= (PI / 180.0);
    betaAngle *= (PI / 180.0);
    gammaAngle *= (PI / 180.0);

    // Initialize unit cell vectors from |a|,|b|,|c| and alphaAngle, betaAngle, gammaAngle information
    // Here we are applying the A along x-axis, B in xy plane convention
    aV[0] = aLength; aV[1] = 0; aV[2] = 0;
    bV[0] = bLength*cos(gammaAngle); bV[1] = bLength*sin(gammaAngle); bV[2] = 0;
    cV[0] = cLength*cos(betaAngle);
    cV[1] = (cLength*bLength*cos(alphaAngle) - bV[0]*cV[0])/bV[1];
    cV[2] = sqrt(cLength*cLength - cV[0]*cV[0] - cV[1]*cV[1]);

    if (useEwardSums == true) DetermineReciprocalLatticeVectors();

    // Unitcell Volume
    vector<double> crs;
    crs = Cross(bV,cV);
    unitCellVolume = fabs( aV[0]*crs[0] + aV[1]*crs[1] + aV[2]*crs[2] ); // Volume of a parallelipiped

    // Find first line that does not contain underscore
    bool underscoreFound = true;
    int eInd2 = eInd; // we need another index
    while (underscoreFound == true) {
        sInd = eInd2; // End of the previous line
        eInd2 = data.find("\n", eInd2 + 1); // End of the next line
        cStr = data.substr(sInd, eInd2 - sInd); // The line
        if ((cStr.find("_",0) >=0 && cStr.find("_",0) < cStr.size()) || cStr.length() < 20) {
            underscoreFound = true; // Under score found, skip to the next line
        } else {
            // Underscore not found, we are on a legitimate line of data
            underscoreFound = false;
        }
    }

    cerr << "==================================================" << endl;
    cerr << "========= Atom types - X & J values used =========" << endl;
    cerr << "==================================================" << endl;

    Coordinates tempAtom;
    while (underscoreFound == false) {
        if ((cStr.find("_",0) >=0) && (cStr.find("_",0) < cStr.size())) {
            underscoreFound = true; // Under score found, skip to the next line
        } else {
            // Underscore not found, we are on a legitimate line of data
            underscoreFound = false;

            //Read atom label
            iInd = cStr.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789", 0);
            sInd = cStr.find_first_of(" \t", iInd);

            // Handles cases where EOF is hit before underscore found
            if (sInd == string::npos) {
                break;
            }

            tStr = cStr.substr(iInd, sInd-1);
            Label.push_back(tStr);

            // Read atom symbol
            sInd = cStr.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", sInd);
            eInd = sInd + 1;
            tStr = cStr.substr(sInd, eInd - sInd + 1);
            Symbol.push_back(tStr);

            // Find first "x" coordinate
            sInd = cStr.find(".",sInd) - 2;
            eInd = cStr.find_first_of(" \t",sInd + 2);
            tStr = cStr.substr(sInd, eInd - sInd);
            tempAtom.x = atof( tStr.c_str() );    // X Position

            // Find first "y" coordinate
            sInd = cStr.find(".",eInd) - 2;
            eInd = cStr.find_first_of(" \t",sInd + 2);
            tStr = cStr.substr(sInd, eInd - sInd);
            tempAtom.y = atof( tStr.c_str() );    // Y Position

            // Find first "z" coordinate
            sInd = cStr.find(".",eInd) - 2;
            eInd = cStr.find_first_of(" \t\n",sInd + 2);
            tStr = cStr.substr(sInd, eInd - sInd);
            tempAtom.z = atof( tStr.c_str() );    // Z Position

            // Change from fractional to cartesian:
            tempAtom.x = tempAtom.x * aV[0] + tempAtom.y * bV[0] + tempAtom.z * cV[0];
            tempAtom.y = tempAtom.x * aV[1] + tempAtom.y * bV[1] + tempAtom.z * cV[1];
            tempAtom.z = tempAtom.x * aV[2] + tempAtom.y * bV[2] + tempAtom.z * cV[2];

            Pos.push_back(tempAtom);

            int i = Symbol.size() - 1;
            int Z = s_mapStringAtomLabels[Symbol[i]]; // Get Z number from label

            if (Symbol[i] == "H ") {
                X.push_back(0.5*(hI1 + hI0));
                J.push_back(hI1 - hI0);
            } else {
                int cC = IonizationData[Z].chargeCenter;
                X.push_back(0.5*(IonizationData[Z].ionizationPotential[cC+1] +
                    IonizationData[Z].ionizationPotential[cC]));
                J.push_back(IonizationData[Z].ionizationPotential[cC+1] -
                    IonizationData[Z].ionizationPotential[cC]);
                X[i] -= cC*(J[i]);
            }

            bool beenDone = false;
            for (int j = 0; j < i; j++) {
                if (Symbol[i] == Symbol[j]) beenDone = true;
            }
            if (beenDone == false) {
                cerr << Symbol[i] << "\t";
                cerr << "Z: " << Z+1 << "\t";
                cerr << "Ch. Cent: " << IonizationData[Z].chargeCenter << "\t";
                cerr << "X: " << X[i] << "\t";
                cerr << "J: " << J[i] << "\t" << endl;
            }
        }
        sInd = eInd2; // End of the previous line
        eInd2 = data.find("\n", eInd2 + 1); // End of the next line
        if (eInd2 == -1) {
            break;
        }
        cStr = data.substr(sInd, eInd2 - sInd); // The line
    }
    numAtoms = Pos.size();
    Q.resize(numAtoms, 0); // initialize charges to zero
}
/*****************************************************************************/
void OutputFile(string filename, string data) {
    FILE *out;
    out = fopen(filename.c_str(), "wt");
    fprintf(out, "%s", data.c_str());
    fclose(out);
}
/*****************************************************************************/
void OutputCIFFormatFile(string filename) {
    OutputFile(filename, OutputCIFData());
}
/*****************************************************************************/
void OutputPDBFormatFile(string filename) {
    OutputFile(filename, OutputPDBData());
}
/*****************************************************************************/
void OutputMOLFormatFile(string filename) {
    OutputFile(filename, OutputMOLData());
}
/*****************************************************************************/
void OutputCARFormatFile(string filename) {
    OutputFile(filename, OutputCARData());
}
/*****************************************************************************/
void OutputJSONFormatFile(string filename) {
    OutputFile(filename, OutputJSONData());
}
/*****************************************************************************/
string OutputCIFData() {
    ostringstream stringStream;
    stringStream << "data_functionalizedCrystal" << endl;
    stringStream << "_audit_creation_method\t" << "'EQeq! by Chris Wilmer'" << endl;
    stringStream << "_symmetry_space_group_name_H-M\t" << "'P1'" << endl;
    stringStream << "_symmetry_Int_Tables_number\t" << "1" << endl;
    stringStream << "_symmetry_cell_setting\t" << "triclinic" << endl;
    stringStream << "loop_" << endl;
    stringStream << "_symmetry_equiv_pos_as_xyz" << endl;
    stringStream << "  x,y,z" << endl;
    stringStream << "_cell_length_a\t" << aLength << endl;
    stringStream << "_cell_length_b\t" << bLength << endl;
    stringStream << "_cell_length_c\t" << cLength << endl;
    stringStream << "_cell_angle_alpha\t" << alphaAngle * (180 / PI) << endl;
    stringStream << "_cell_angle_beta\t" << betaAngle * (180 / PI) << endl;
    stringStream << "_cell_angle_gamma\t" << gammaAngle * (180 / PI) << endl;
    stringStream << "loop_" << endl;
    stringStream << "_atom_site_label" << endl;
    stringStream << "_atom_site_type_symbol" << endl;
    stringStream << "_atom_site_fract_x" << endl;
    stringStream << "_atom_site_fract_y" << endl;
    stringStream << "_atom_site_fract_z" << endl;
    stringStream << "_atom_site_charge" << endl;

    // For all atoms
    for (int i = 0; i < numAtoms ; i++) {
        // Determine the fractional coordinates
        double dx = Pos[i].x;
        double dy = Pos[i].y;
        double dz = Pos[i].z;

        // Convert to fractional coordinates (below is the "inverse transform matrix")
        double a = (bV[2]*cV[1]*dx - bV[1]*cV[2]*dx - bV[2]*cV[0]*dy + bV[0]*cV[2]*dy + bV[1]*cV[0]*dz - bV[0]*cV[1]*dz)/
                   (aV[2]*bV[1]*cV[0] - aV[1]*bV[2]*cV[0] - aV[2]*bV[0]*cV[1] +
                    aV[0]*bV[2]*cV[1] + aV[1]*bV[0]*cV[2] - aV[0]*bV[1]*cV[2]);
        double b = (aV[2]*cV[1]*dx - aV[1]*cV[2]*dx - aV[2]*cV[0]*dy + aV[0]*cV[2]*dy + aV[1]*cV[0]*dz - aV[0]*cV[1]*dz)/
                   (-(aV[2]*bV[1]*cV[0]) + aV[1]*bV[2]*cV[0] + aV[2]*bV[0]*cV[1] -
                   aV[0]*bV[2]*cV[1] - aV[1]*bV[0]*cV[2] + aV[0]*bV[1]*cV[2]);
        double c = (aV[2]*bV[1]*dx - aV[1]*bV[2]*dx - aV[2]*bV[0]*dy + aV[0]*bV[2]*dy + aV[1]*bV[0]*dz - aV[0]*bV[1]*dz)/
                   (aV[2]*bV[1]*cV[0] - aV[1]*bV[2]*cV[0] - aV[2]*bV[0]*cV[1] +
                   aV[0]*bV[2]*cV[1] + aV[1]*bV[0]*cV[2] - aV[0]*bV[1]*cV[2]);

        stringStream<< "Mof_" << Symbol[i] << "\t" << Symbol[i] << "\t";
        stringStream << a << "\t" << b << "\t" << c << "\t" << Q[i] << endl;
    }

    stringStream << "_end" << endl;
    return stringStream.str();
}
/*****************************************************************************/
string OutputPDBData() {
    ostringstream stringStream;
    char buf[200];

    stringStream << "TITLE       YourMoleculeNameHere            " << endl;
    stringStream << "REMARK   4" << endl;
    stringStream << "REMARK   4      COMPLIES WITH FORMAT V. 2.2, 16-DEC-1996" << endl;
    if (isPeriodic) {
        sprintf(buf, "CRYST1    %5.2f    %5.2f    %5.2f  %3.2f  %3.2f  %3.2f P1\n",
            aLength,bLength,cLength,alphaAngle*180/PI,betaAngle*180/PI,gammaAngle*180/PI);
        stringStream << buf;
    }
    for (int i = 0; i < numAtoms; i++) {
        sprintf(buf, "ATOM    %3d %s   MOL A   0     % 7.3f % 7.3f % 7.3f % 5.2f                %s\n",
            i+1,Symbol[i].c_str(),Pos[i].x,Pos[i].y,Pos[i].z,Q[i],Symbol[i].c_str());
        stringStream << buf;
    }
    return stringStream.str();
}
/*****************************************************************************/
string OutputJSONData() {
    ostringstream stringStream;
    stringStream << "[";
    for (int i = 0; i < numAtoms - 1; i++) {
        stringStream << Q[i] << ",";
    }
    stringStream << Q[numAtoms - 1] << "]";
    return stringStream.str();
}
/*****************************************************************************/
string OutputMOLData() {
    ostringstream stringStream;
    char buf[200];

    stringStream << " Molecule_name: hypotheticalMOF" << endl << endl;
    stringStream << "  Coord_Info: Listed Cartesian None" << endl;
    stringStream << "        " << numAtoms << endl;

    for (int i = 0; i < numAtoms; i++) {
        sprintf(buf, "  %4d  % 8.4f % 8.4f % 8.4f  Mof_%s   % 6.3f  0  0\n",
                i + 1,Pos[i].x, Pos[i].y, Pos[i].z, Symbol[i].c_str(), Q[i]);
        stringStream << buf;
    }

    stringStream << endl << endl << endl;
    stringStream << "  Fundcell_Info: Listed" << endl;
    sprintf(buf, "        %8.4f      %8.4f      %8.4f\n", aLength, bLength, cLength);
    stringStream << buf;
    sprintf(buf, "        %8.4f      %8.4f      %8.4f\n", alphaAngle * 180 / PI,
            betaAngle * 180 / PI, gammaAngle * 180 / PI);
    stringStream << buf;
    stringStream << "        0.00000        0.00000       0.00000" << endl;
    sprintf(buf, "        %8.4f      %8.4f      %8.4f\n", aLength, bLength, cLength);
    stringStream << buf;

    return stringStream.str();
}
/*****************************************************************************/
string OutputCARData() {
    ostringstream stringStream;
    char buf[200];
    
    stringStream << "!BIOSYM archive 3" << endl;
    stringStream << "PBC=ON" << endl;
    stringStream << "Generating .CAR file" << endl;
    stringStream << "!DATE" << endl;
    sprintf(buf, "PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%7s\n", aLength, bLength, cLength,
            alphaAngle * 180 / PI, betaAngle * 180 / PI, gammaAngle * 180 / PI, "(P1)");
    stringStream << buf;

    for (int i = 0; i < numAtoms; i++) {
        sprintf(buf, "%-5s %14.9f %14.9f %14.9f XXXX %-7d %-7s %2s %6.3f\n", Symbol[i].c_str(),
                Pos[i].x, Pos[i].y, Pos[i].z, 1, "xx", Symbol[i].c_str(), Q[i]);
        stringStream << buf;
    }
    stringStream << "end" << endl << "end";

    return stringStream.str();
}
/*****************************************************************************/
void Qeq() {
    int i, j; // generic counter;

    // Formulate problem in the form of A x = b
    vector<double> dummyRow(numAtoms, 0); // Is doing this necessary?
    vector<vector<double> > A(numAtoms, dummyRow);
    vector<double> b(numAtoms,0);

    // First row of A is all ones
    for (int i = 0; i < numAtoms; i++) {
        A[0][i] = 1;
    }

    // First element in b is the total charge
    b[0] = Qtot;

    // Rest of elements in b are the differences in electronegativity
    for (int i = 1; i < numAtoms; i++) {
        b[i] = X[i] - X[i-1];
    }

    // Fill in 2nd to Nth rows of A
    for (int i = 1; i < numAtoms; i++) {
        cerr << "." << flush;
        for (int j = 0; j < numAtoms; j++) {
            A[i][j] = GetJ(i-1, j) - GetJ(i, j);
        }
    }

    Q = SolveMatrix(A,b);
}
/*****************************************************************************/
void RoundCharges(int digits) {

    double qsum = 0;
    double factor = pow((double)10,digits);

    for(int i=0; i < numAtoms; i++) {
        Q[i] = Round(Q[i]*factor)/factor;
        qsum += Q[i];
    }

    if (qsum == 0) { // Great, rounding worked on the first try!
        // do nothing
    } else { // There is a small excess charge from rounding, adjust it
        int numAtomsToAdjust = (int)(abs(qsum * factor) + 0.5); // Weird double-to-int conversion tricks
        cerr << " adjusting the charge of " << numAtomsToAdjust << " atoms!" << endl;

        int sign; if (qsum > 0) sign = -1; else sign = 1;
        for (int i=0; i < numAtomsToAdjust; i++) { // Adjust
            Q[i] += sign*(1/factor);
        }
    }

}
/*****************************************************************************/
vector<double> Cross(vector<double> a, vector<double> b) {

    vector<double> c(3);

    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

    return c;
}
/*****************************************************************************/
double Dot(vector<double> a, vector<double> b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
/*****************************************************************************/
double Mag(vector<double> a) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
/*****************************************************************************/
double Round(double num) {
    return (num > 0.0) ? floor(num + 0.5) : ceil(num - 0.5);
}
/*****************************************************************************/
vector<double> Scalar(double a, vector<double> b) {
    vector<double> c(3);
    c[0] = a*b[0]; c[1] = a*b[1]; c[2] = a*b[2];
    return c;
}
/*****************************************************************************/
vector<double> SolveMatrix(vector<vector<double> > A, vector<double> b) {
    // Assumptions: A x = b, A is a MxN matrix, M = rows, N = cols, x is a vector, b is vector
    // matrix has more rows than columns
    // number of rows of matrix is equal to size of the vector x

    // Initialize x = b
    vector<double> x;
    x = b;

    int i, j, k;
    int N = A.size();
    int M = A[0].size();

    vector<double> d (N);

    /* Perform Householder transformation */
    for (i = 0; i < N; i++) {
        const double aii = A[i][i];
        double alef, f, ak;
        double max_norm = 0.0;
        double r = 0.0;

        for (k = i; k < M; k++) {
          r += A[k][i] * A[k][i];
        }

        if (r == 0) {
          cerr << "Error! Matrix is rank deficient." << endl;
          // return -1;
        }

        if (A[i][i] < 0)
          alef = (-1)*sqrt(r);
        else
          alef = sqrt(r);

        ak = 1.0 / (r + alef * A[i][i]);

        A[i][i] +=  alef;

        d[i] = -alef;

        for (k = i + 1; k < N; k++) {
          double norm = 0.0;
          f = 0.0;

          for (j = i; j < M; j++) {
            norm += A[j][k] * A[j][k];
            f += A[j][k] * A[j][i];
          }

          max_norm = max(max_norm, norm);

          f *= ak;

          for (j = i; j < M; j++) {
            A[j][k] -= f * A[j][i];
          }
        }

        if (fabs(alef) < 0.00001) {
          cerr << "Apparent singularity in matrix." << endl;
          // return -1;
        }

        f = 0.0;

        for (j = i; j < M; j++) {
          f += x[j] * A[j][i];
        }

        f *= ak;

        for (j = i; j < M; j++) {
          x[j] -= f * A[j][i];
        }
    }

    /* Perform back-substitution */

    for (i = N-1; i >= 0; i--) {
      double sum = 0.0;

      for (k = i + 1; k < N; k++) {
        sum += A[i][k] * x[k];
      }

      x[i] = (x[i] - sum) / d[i] ;
    }

    return x;
}
/*****************************************************************************/
