#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H 1

#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>

#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TH1F.h>
#include <TPostScript.h>
#include <TPDF.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TSpline.h>
#include <TPaveStats.h>
#include <TLatex.h>
// // 	
#define DEBUG 1
#undef DEBUG
#define DEBUGMSG 1
#undef DEBUGMSG
#define SLOWFILES 1
// #undef SLOWFILES

#define PeakparamDef(name,id)
using namespace std;

class PEAKPARAM : public TObject {
public:
   int maxtime_pos;
   int stime_pos;
   int ftime_pos;
   int e_peak_end_pos;
   int sig_start_pos;
   int sig_end_pos;
   int tot_sig_end_pos;
   
   double maxtime;
   double ampl;
   double e_peak_end_ampl;
   double sampl;
   double fampl;
   //double t20;
   //double st20;
   double tfit20;
   double tnaive20;
   double te_peak_end;
   
   //double sechargefixed;
   //double secharge;
   double echargefixed;
   //double echarge;
   double echargefit;
   double totchargefixed;
   double risetime;  ///10% - 90%
   double risecharge;
   double width;
   double tot[10]; ///
   double sigmoidR[4];
   double sigmoidF[4];
   double sigmoidtot[6];

   double charge;
   double scharge;
   double t10;
   double tb10;
   double t90;

   double ttrig;
   double bslch;
   double rms;
   double bsl;
 
  PEAKPARAM() {}

//  PeakparamDef(PEAKPARAM,1)
  ClassDef(PEAKPARAM,1)
};
 
typedef struct {
//    float corr;
   int npeaks;
   int maxtime;
   
   int stime;
   int ftime;
   double *ampl;
  
   double *sampl;
   double *fampl;
   double *t20;
   double *st20;
   double *tfit20;

   double *sechargefixed;
   double *secharge;
   double *echargefixed;
   double *echarge;
   double *echargefit;
   double *totchargefixed;
   double *risetime;  ///10% - 90%
   double *risecharge;
   double *width;
   double *tot[10]; ///
   double *sigmoidR[3];
   double *sigmoidF[3];

   double *charge;
   double *scharge;
   double *t10;
   double *tb10;
   double *t90;

   double *ttrig;
   double *bslch;
   } PICPARAM;


typedef struct {
//    float corr;
   double bsl;
   double rms;
   double ftime;
   double stime;
   double maxtime;
   double maxtimeS;
   double time50;
   double ampl;
   double sampl;
   double intg;
   double charge;
   double chargeIon;
   double chi2;
   double t1;
   double t2;
   double t3;
   double a;  /// fit function f = a*x + b
   double b;
   double ttrig;
   } DPARAM;

typedef struct {
//    float corr;
   double bsl;
   double rms;
   int stime;
   int ftime;
   int maxtime;
   double tot;
   double ampl;
   double charge;
   double risecharge;
   double t10;
   double tb10;
   double t90;
   double tmax;
   double ttrig;
   double width;
   double sampl;
   double fampl;
   double bslch;
   } IPARAM;

typedef struct {
//    float corr;
   int npeaks;
   double *tot;
   double *ampl;
   double *charge;
   double *risecharge;
   double *t10;
   double *tb10;
   double *t90;
   double *ttrig;
   double *width;
   double *sampl;
   double *fampl;
   double *bslch;
   } PPARAM;

typedef struct {
//    float corr;
   int detNo;
   int preamNo;
   int runNo;
   double ampl;
   double sampl;
   double charge;
   double scharge;
   double rate;
   double srate;
   double grate;
   double sgrate;
   double tot;
   double stot;
   double risetime;
   double srisetime;
   double width;
   double swidth;
   double chovampl;
   double schovampl;
   } RUNPAR;

typedef struct {
//    float corr;
   int runNo;
   int poolNo;
   int preamNo[4];
   double tsigma[4];
   double trms[4];
   double amplpolya[4][3];
   double echargepolya[4][3];
   double chargepolya[4][3];
   double risetime[4];
   double width[4];
   double chovampl[4];
   } RUNPAR4;

typedef struct {
//    float corr;
   int runNo;
   int poolNo;
   int srsCh;
   float V1[4];
   float V2[4]; 
   float Z[4]; 
   char DetName[4][10];
   char Photocathode[4][12];
   char Amplifier[4][20];
   int AmplifierNo[4];
} OSCSETUP;

typedef struct {
//    float corr;
   int srstriggerctr[20];   /// [i] = trackNumber
   int srstimestamp[20];    /// [i] = trackNumber
   int ntracks;
   float trackchi2[20];    /// [i] = trackNumber
   float hits[20][20][3];   /// [i] = trackNumber  [j] = position  [k] = x,y,z
   float distnextcluster[20][6];  /// [i] = trackNumber [j] = GEM xy
   float totchanextcluster[20][6];  /// [i] = trackNumber [j] = GEM xy
} TRACKDATA;

typedef struct{
    double ampl;
    int pos;
    float x;
} GLOBALMAXIMUM;

typedef struct {
   double ampl; 
   int pos; 
   double x; 
   float dt;   
} STARTPOINT;

typedef struct {
    double ampl;
    int pos;
    double x;
    float dt;
    int arrpoints;
} ENDPOINT;


typedef struct {
    double ampl; 
    int pos; 
    double x;
} SEARCHEPEAKENDPOINT;

const double rootConv = 2871763200.;
const double unixConv = 2082844800.;

const int MAX_N_FILES=11000;
//const char *CODEDIR="/drf/projets/picosecond/Saclay/cpp_analysis/code";
//const char *BASEDIRNAME="/drf/projets/picosecond/Saclay/data/2022_October_h4";
//const char *WORKDIR="/drf/projets/picosecond/Saclay/data/2022_October_h4/wdir";
//const char *PLOTDIR="/drf/projets/picosecond/Saclay/data/2022_October_h4/plots";
//const char *DATADIRNAME="/drf/projets/picosecond/Saclay/data/2022_October_h4/dataTrees";
//const char *TRACKDIRNAME="/drf/projets/picosecond/Saclay/data/2022_October_h4/tracking";
//const char *OUTDIRNAME="/drf/projets/picosecond/Saclay/data/2022_October_h4/processedTrees";
//const char *PARAMDIRNAME="/drf/projets/picosecond/Saclay/data/2022_October_h4/processedTrees/ParameterTrees";
// const char *TRACKDIRNAME="/eos/project-p/picosec/testbeam/2022_October_h4/tracker/reconstructed";

 //const char *CODEDIR="/eos/project-p/picosec/analysis/Saclay/code";
 //const char *BASEDIRNAME="/eos/project-p/picosec/analysis/Saclay/data/2023_April_h4";
 //const char *WORKDIR="/eos/project-p/picosec/analysis/Saclay/data/2023_April_h4/wdir";
 //const char *PLOTDIR="/eos/project-p/picosec/analysis/Saclay/data/2023_April_h4/plots";
 //const char *DATADIRNAME="/eos/project-p/picosec/analysis/Saclay/data/2023_April_h4/dataTrees";
 //const char *TRACKDIRNAME="/eos/project-p/picosec/testbeam/2023_April_h4/tracker/reconstructed";
 //const char *OUTDIRNAME="/eos/project-p/picosec/analysis/Saclay/data/2023_April_h4/processedTrees";
 //const char *PARAMDIRNAME="/eos/project-p/picosec/analysis/Saclay/data/2023_April_h4/processedTrees/ParameterTrees";

const char *CODEDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/cpp_analysis/code";
const char *BASEDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4";
const char *WORKDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/wdir";
const char *PLOTDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots";
const char *DATADIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/dataTrees";
const char *TRACKDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/tracking";
const char *OUTDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees";
const char *PARAMDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees/ParameterTrees";

// const char *CODEDIR="/home/dylan/CLionProjects/picosec_cpp_analysis/code";
// const char *BASEDIRNAME="/home/dylan/Desktop/picosec/data/2022_October_h4";
// const char *WORKDIR="/home/dylan/Desktop/picosec/data/2022_October_h4/wdir";
// const char *PLOTDIR="/home/dylan/Desktop/picosec/data/2022_October_h4/plots";
// const char *DATADIRNAME="/home/dylan/Desktop/picosec/data/2022_October_h4/dataTrees";
// const char *TRACKDIRNAME="/home/dylan/Desktop/picosec/data/2022_October_h4/tracking";
// const char *OUTDIRNAME="/home/dylan/Desktop/picosec/data/2022_October_h4/processedTrees";
// const char *PARAMDIRNAME="/home/dylan/Desktop/picosec/data/2022_October_h4/processedTrees/ParameterTrees";

//const char *CODEDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/cpp_analysis/code";
//const char *BASEDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4";
//const char *WORKDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4/wdir";
//const char *PLOTDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4/plots";
//const char *DATADIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4/dataTrees";
//const char *TRACKDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4/tracking";
//const char *OUTDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4/processedTrees";
//const char *PARAMDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_April_h4/processedTrees/ParameterTrees";

//const char *CODEDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/cpp_analysis/code";
//const char *BASEDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4";
//const char *WORKDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4/wdir";
//const char *PLOTDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4/plots";
//const char *DATADIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4/dataTrees";
//const char *TRACKDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4/tracking";
//const char *OUTDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4/processedTrees";
//const char *PARAMDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2023_August_h4/processedTrees/ParameterTrees";


const char *RTYPE="TESTBEAM";
const int MINRUN = 1;
const int MAXRUN = 1000;
const int TOTRUNS = MAXRUN-MINRUN+1;

const double Threshold = 0.0073;

//const int N_INTEGRATION_POINTS = 10; default
const int N_INTEGRATION_POINTS = 20;

//const double CIVIDEC_PULSE_DURATION = 75; // [ns]
const double CIVIDEC_PULSE_DURATION = 150; // [ns]
const double CIVIDEC_PEAK_DURATION = 6; // [ns]


const double SIGMOID_EXTENTION = 3; //[ns]
/// definitions for the bin2tree.cxx ___________
#ifndef PATH_NAMES_DATA_CODE
#define PATH_NAMES_DATA_CODE 1
//const char *DATA_PATH_NAME="/eos/project-p/picosec/testbeam/2022_October_h4/";
    const char *DATA_PATH_NAME ="/sw/akallits/PicoAnalysis/TestBeams/2022_October_h4";
    //const char *DATA_PATH_NAME ="/sw/akallits/PicoAnalysis/TestBeams/2023_August_h4";

    //const char *DATA_PATH_NAME="/eos/project-p/picosec/testbeam/2023_April_h4/";

    const char *OUT_DIR_NAME=BASEDIRNAME;  /// Must be the same with above
    const char *WORK_DIR_NAME=WORKDIR;   /// Must be the same with above
    const char *RUN_TYPE="TESTBEAMraw";

    const int RUNMIN=MINRUN;
    const int RUNMAX=MAXRUN;

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define KRESET "\x1B[0m"

#endif
///_______________________________________________



const std::string RED("\033[0;31m"); 
const std::string RED_U_GBKG("\033[4;31;42m"); 
const std::string GREEN("\033[0;32m"); 
const std::string YELLOW("\033[0;33m"); 
const std::string BLUE("\033[0;34m"); 
const std::string MAGENTA("\033[0;35m"); 
const std::string CYAN("\033[0;36m"); 
const std::string INVERSE_ON("\033[7m"); 
const std::string INVERSE_OFF("\033[27m"); 
const std::string RESET_COLOR("\033[0m");	
const std::string endlr("\n\033[0m");	

inline int month2int_littleendian(const char *month)
{
  const char hash[] = { 3, 12, 8, 2, 1, 11, 7, 5, 0, 10, 4, 9, 6 };
  return hash[(* (const int32_t *) month & ~0x20202020 + 146732) % 13];
}


inline void replaceEOL(char *fnames)
{
    int nchar=0;
    while (fnames[nchar]!='\n')
    {
      nchar++;
    }
    fnames[nchar]='\0';
}  

inline void replaceEOLT(char *fnames)
{
    int nchar=0;
    while (fnames[nchar]!='\n')
    {
      nchar++;
    }
    fnames[nchar]='\t';
//     fnames[nchar+1]='\0';
}  

inline int FindZeros(int n, double* data)
{
   int n0=0;
   for (int i=1;i<n;i++)
   {
     if (data[i]==0.0)
       if (data[i-1]==0)
	 n0 ++; 
//      cout<<"["<<data[i]<<"/"<<(data[i]==0)<<"] ";
   }
   return n0;
}

inline void WhichAmplifier(int amplifierNo, char* amplifier)
{
   if (amplifierNo==1)
    sprintf(amplifier,"cividec");
  else if (amplifierNo==2)
    sprintf(amplifier,"ortec");
  else if (amplifierNo==3)
    sprintf(amplifier,"ATHL");
  else if (amplifierNo==4)
    sprintf(amplifier,"CERN");
  else 
    sprintf(amplifier,"other");
}

inline int WhichAmplifier(const char* amplifier)
{
  int amplifierNo=-1;
  cout<<MAGENTA<<"Searching for amplifier ==>"<<amplifier<<"<=="<<endlr;
    
  if (strcmp(amplifier,"cividec")==0)
  {
      amplifierNo=1;
//       cout<<RED<<"A cividec is used"<<endl;
  }
  else if (strcmp(amplifier,"ortec")==0)
      amplifierNo=2;
  else if (strcmp(amplifier,"ATHL")==0)
      amplifierNo=3;
  else if (strcmp(amplifier,"CERN")==0)
      amplifierNo=4;
  else 
      amplifierNo=1234;
  return (amplifierNo);
}

int FilterHisto(TH1* , double );

int SubtractBaseline(int , double* , double* , double );

int GetOsciloscopeSetup(int , int , const char* , OSCSETUP*);

double IntegrateA(int , double* , double* ,double );

//void FindGlobalMaximum(double, int, float);
double FindBaselineLevel(double*, int);
double fermi_dirac_general(double*, double* );
double fermi_dirac(double*, double *);
double fermi_dirac_generalsub(double *, double *);
double slope_at_x(double*, int , int , double);


GLOBALMAXIMUM FindGlobalMaximum(double *, int , double *);

STARTPOINT FindStartPoint(double *, double *, double , int , float , STARTPOINT );

ENDPOINT FindEndPoint(double *, double *, double ,  int , float,  int , GLOBALMAXIMUM );

SEARCHEPEAKENDPOINT SearchEpeakEndPoint(double *, double *, GLOBALMAXIMUM, ENDPOINT, int );

#endif
