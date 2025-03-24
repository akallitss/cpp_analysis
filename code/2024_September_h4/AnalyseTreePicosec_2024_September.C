#pragma once

#if defined __GNUC__
#   pragma GCC system_header
#elif defined __SUNPRO_CC
#   pragma disable_warn
#elif defined _MSC_VER
//# pragma warning(push, 0)
#   pragma warning(push)
#   pragma warning(disable : 4800)
//# pragma warning(disable : ...)
#endif
/////////up to here for disabling warnings
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
 #include <TH1F.h>
// #include <TH2F.h>
// #include <TH2D.h>
// #include <TH3F.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TExec.h>
#include <TLine.h>
#include <TObjString.h>
#include<TMultiGraph.h>
#include<TClonesArray.h>
#include "MyFunctions_2024_September.C"
#include "../RMS_Baseline_Calculator/RMSBaselineCalculator_2024_September.cpp"
#include <iomanip>
#include <iostream>



int AnalyseTreePicosec_2024_September(int runNo=15, int poolNo=2, int draw=0, double threshold = 10.0, double peTh = 12.0, bool isCalibration = 0, string filetype = "")
{
  cout<<"RunNo= " <<runNo<< " poolNo = "<<poolNo<<endl;

  gROOT->LoadMacro("MyFunctions_2024_September.C");
  
  bool activeDraw[]={0,0,0,0};
  if (draw)
  {
     int whichchannel=0;
     cout<<RED<<"  which channel is active for drawing? \n 1-4 or 0 for all --------> " ;
     cin>>whichchannel;
     if (whichchannel>=1 && whichchannel<=4)
       activeDraw[whichchannel-1]=kTRUE;
     else
       for (int i=0;i<4;i++)
       activeDraw[i]=kTRUE;
  }
  

  threshold/=1000.;   /// make it in V
  peTh/=1000.;

   const double mV=1000.;   /// set mV=1. to make everything in volts!


  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");
  gStyle->SetLabelFont(132,"X");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetTitleSize(0.045,"X");
  gStyle->SetTitleSize(0.045,"Y");
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetOptStat(1001111);
  gStyle->SetNdivisions(507);
  
  int exttrig=0;
//   double peTh = 0.010;
//   peTh = threshold;
  
  int preamNo = 3;
  
  int vm[4];
  int vd[4];
  int dv[4];

  if(threshold>0)//always negative
    threshold *= -1;
  
  if(peTh>0)
    peTh *= -1;

  peTh = fabs(peTh);

  char dirname[1000];
  char anadirname[1000];
  char command[1500];
  char fname[500];
  char ofname[500];
  char fname2[500];
  char fnames[10000][200];
  char fnames2[10000][200];
  char detector[100];
  
  
  int runb = 0;
  char basedirname[1000];
  char workdirname[1000];

  char treename[200];
  char treetitle[200];

  char afile[500]; 
  char ftypetmp[500];
  char fnametmp[500];
  

  sprintf(basedirname,"%s",OUTDIRNAME);
  sprintf(anadirname,"%s",PARAMDIRNAME);
  sprintf(workdirname,"%s",WORKDIR);
  // // // 

  if (runNo>MINRUN && runNo<=MAXRUN)
  {
      sprintf(afile,"%s/tmpfile.tmp",workdirname);
      FILE *ftmp=fopen(afile,"w");
      if (ftmp == NULL)
      {
        cout<<afile<<" can not be created. Probablly the directory '"<<basedirname<<"' does not exit. Exiting..."<<endl;
        exit (-12);
      }
      fclose(ftmp);
      if (poolNo==0)
         sprintf(command,"cd %s\nls -d Run%03d-GDD*%s*%s*_tree.root > %s 2>/dev/null",basedirname,runNo,filetype.c_str(),RTYPE,afile);
      else
         sprintf(command,"cd %s\nls -d Run%03d-Pool%d*%s*%s*_tree.root > %s 2>/dev/null",basedirname,runNo,poolNo,filetype.c_str(),RTYPE,afile);
      int tst=system(command);
      //    cout<<command<<endl<<"returned: "<<tst<<endl;
      if (tst !=0)
      {
          cout<<command<<endl<<"returned: "<<tst<<endl;
          cout<<"Probably the tree of run No "<<runNo<<" , pool No "<<poolNo<<" was not found in directory "<<basedirname<<endl<<"Exiting..."<<endl;
          return tst;
      }
      ftmp=fopen(afile,"r");
      if (fgets(fnametmp,200,ftmp) == NULL)
	  {
	     cout<<"Failed to read the filename for run "<<runNo<<" , pool No "<<poolNo<<" at "<<basedirname<<endl<<"Exiting..."<<endl;
	     return -2;
	  }
      int rtmp,dtmp;
      strcpy(ftypetmp,"");
      int stst;
      if (poolNo==0)
      {
         stst = sscanf(fnametmp,"Run%3d-GDD_%30[^ /,\n\t]",&rtmp,ftypetmp);
         dtmp=0;
      }
      else
         stst = sscanf(fnametmp,"Run%3d-Pool%d_%30[^ /,\n\t]",&rtmp,&dtmp,ftypetmp);
//       int stst = sscanf(fnametmp,"S%03d-%d-%d-%3f-%3f%30s",&rtmp,&vm,&vd,&Rd,&dgap,ftypetmp);
      cout <<"___________________________________________________\n\narguments read = "<<stst<<endl;
      filetype.assign(ftypetmp);
      cout<<"runNo = "<<rtmp<<" ("<<runNo<<")"<<endl;
      cout<<"poolNo = "<<dtmp<<endl;
      //cout<<BLUE<<"input file ==>"<<fnametmp<<"<=="<<endlr;

      //       cout<<BLUE<<"filetype -->"<<ftypetmp<<"<--"<<endlr;

      if (stst>2) 
          cout<<"filetype = "<<filetype<<endl;
      fclose(ftmp);
      sprintf(command,"rm %s\n",afile);
      cout <<"executing command -->"<<command<<"<--"<<endl;
      tst=system(command);
    }
    else
    {
      cout<<"Available runs: "<<MINRUN<<" - "<<MAXRUN<<endl;
      return (-2);
    }
      //cout <<"executed command -->"<<command<<"<--"<<endl;
    
    
  const char *cftype = filetype.c_str();  /// add here any directory supplement
  //cout<<MAGENTA<<"-->"<<cftype<<"<--"<<endlr;
  
  char ftype[1000];
  strcpy (ftype,cftype);
  char *pch;
  char tmptype[100];
  sprintf(tmptype,"%s_tree.root",RTYPE);
  pch = strstr(ftype,tmptype);
  //   sprintf(ofname,"%s_tree.root",pch);
  *pch = '\0';
 
  char poolid[20];
  sprintf(poolid,"Pool%d",poolNo);
  char runid[20];
  sprintf(runid,"Run%03d",runNo);
//   char distance[20];
//   sprintf(distance,"%3.1fm",0.);
  char rtype[1000];
  strcpy(rtype,ftype);

  sprintf(fname,"%s/%s",basedirname,fnametmp);
  
  replaceEOL(fname);
  //cout<<BLUE<<"Input filename =>"<<fname<<"<="<<endlr;

	TFile *ifile = new TFile(fname);

	// Calculate rms/baselines
  
  if (!ifile->IsOpen())
  {
    cout<<RED<<"Attention! File \n"<<fname<<"\ncorresponding to the run "<<runid<<" Does not exist!!!"<<endlr;
    return (-3);
  }

/// read the oscilloscope info tree and store the parameters in "oscsetup" structure
  cout<<GREEN<<"Reading info tree..."<<endlr;
  TTree *infotree;
  
  sprintf(treename,"OsciloscopeSetup");
  sprintf(treetitle,"Osciloscope Setup");
  infotree = (TTree*) ifile->Get(treename);
   
  if (infotree==NULL)
  {
      cout<<RED<<treetitle<<" tree does not exist in file "<<fname<<endl;
      cout<<"Exiting..."<<endlr;
      return (-4);
  }
  
  char bnametmp[100];
  int srsNo;
  sprintf(bnametmp,"srsCh");
  TBranch* branch = infotree->GetBranch(bnametmp);
  if (branch==NULL) cout<<"Branch "<<bnametmp<<" not found"<< endl;  
  branch->SetAddress(&srsNo);
  float V1, V2, ZZ;
  sprintf(bnametmp,"V1");
  branch = infotree->GetBranch(bnametmp);
  if (branch==NULL) cout<<"Branch "<<bnametmp<<" not found"<< endl;  
  branch->SetAddress(&V1);
  branch = infotree->GetBranch("V2");
  branch->SetAddress(&V2);
  branch = infotree->GetBranch("Z");
  branch->SetAddress(&ZZ);
  cout<<GREEN<<"Info tree OK..."<<endlr;

  char DetName[10];
  branch = infotree->GetBranch("DetName");
  branch->SetAddress(DetName);
  char Photocathode[12];
  branch = infotree->GetBranch("Photocathode");
  branch->SetAddress(Photocathode);
  char Amplifier[20];
  branch = infotree->GetBranch("Amplifier");
  branch->SetAddress(Amplifier);
  int amplifierNo;
  branch = infotree->GetBranch("AmplifierNo");
  branch->SetAddress(&amplifierNo);

  bool MCP[4]={0,0,0,0};
  bool MM[4]={0,0,0,0};
  bool SIGMOIDFIT[4]={0,0,0,0};
  bool DOUBLESIGMOIDFIT[4]={0,0,0,0};
  // double sig_tshift[4]={1.0,1.0,1.0,1.0};
  double sig_tshift[4]={0.0,0.0,0.0,0.0};

  OSCSETUP* oscsetup = new OSCSETUP;
  for (int i=0;i<4;i++)
  {
      infotree->GetEntry(i);
      oscsetup->srsCh = srsNo;
      oscsetup->V1[i]=V1;
      oscsetup->V2[i]=V2;
      oscsetup->Z[i]=ZZ;
      strcpy(oscsetup->DetName[i],DetName);
      strcpy(oscsetup->Photocathode[i],Photocathode);
      strcpy(oscsetup->Amplifier[i],Amplifier);
      oscsetup->AmplifierNo[i]=amplifierNo;
      
      cout<<GREEN<<"Detector = "<<DetName<<endl;  

      if(strncmp(DetName,"MCP",3)==0)
      {
        MCP[i]=kTRUE;
        sig_tshift[i]=0.1;  //ns
        SIGMOIDFIT[i]=kTRUE;
      	DOUBLESIGMOIDFIT[i]=false;
      }
      if(strncmp(DetName,"MM",2)==0)
      {
        MM[i]=kTRUE;
        sig_tshift[i]=0.1;  //ns
        SIGMOIDFIT[i]=kTRUE;
      	DOUBLESIGMOIDFIT[i]=true;
      }
   	  if (amplifierNo==1 || amplifierNo == 3)
      {
         SIGMOIDFIT[i]=kTRUE;
   	  	 DOUBLESIGMOIDFIT[i]=true;
         sig_tshift[i]=0.5; //ns
      }

      cout<<GREEN<<"photocathode = "<<Photocathode<<endl;    
      cout<<"amplifier = "<<Amplifier<<endl;
      vm[i]=V1;
      vd[i]=V2;
      dv[i]=V2+V1;
      cout<<"Vm = "<<vm[i]<<" , ";
      cout<<"Vd = "<<vd[i]<<endlr;
  }

  bool callibrationRun=kFALSE;

  cout<<BLUE<<"SRS channel = "<<oscsetup->srsCh <<endlr;
  if (oscsetup->srsCh<1 || oscsetup->srsCh>4)
  {
     cout<<RED<<"No or invalid SRS channel. Treat as callibration run."<< endlr;
     callibrationRun = kTRUE;
  }
  
  
  cout<<GREEN<<"Info tree OK..."<<endlr;
    
  RUNPAR *runpar = new RUNPAR;
  RUNPAR4 *runpar4 = new RUNPAR4;
  runpar4->poolNo = poolNo;
  runpar4->runNo = runNo;
  
  for(int i=0;i<4;i++)
  {
     runpar4->preamNo[i]=oscsetup->AmplifierNo[i]; 
  }
  
//   char mesh[1000];
//   sprintf(mesh,"%d",vm);
//   char drift[1000];
//   sprintf(drift,"%d",vd);
  
  
  cout <<"********************\n\n\n Starting processing input tree now...\n\n"<< endl;  
  cout<<"\n\n\nWorking directory: "<<workdirname<<endl<<endl;
  
  
  ///OPEN the data tree here 

  TTree *tree;
  sprintf(treename,"RawDataTree");
  sprintf(treetitle,"Raw data tree");
  tree = (TTree*)ifile->Get(treename);

  if (tree==NULL) 
  {
      cout<<RED<<"Data tree from file "<<fname<<" is NULL!\nExiting..."<<endlr;
      return -5;
  }
  ///reading osciloscope parameters
  
  
  double dt = 0.1; /// 0.1 ns sampling step by defauld. The actual value will be obtained from the tree

  char cthreshold[20];
  sprintf(cthreshold,"%.1fmV",fabs(threshold*1000.));
  char cpethreshold[20];
  sprintf(cpethreshold,"%.1fmV",fabs(peTh*1000.));
  
  TString rtypes(runid);
  rtypes+="_";
  rtypes+=poolid;
//   rtypes+="_";
//   rtypes+=mesh;
//   rtypes+="_";
//   rtypes+=drift;
//   rtypes+="_";
//   rtypes+=distance;
//   rtypes+="_";
//   rtypes+=cthreshold;  
  
  char catt[30];
  
  if (runb)
    rtypes+="_b";

  TString ctypes(runid);
  ctypes+="_";
  ctypes+=poolid;
//   ctypes+="_";
//   ctypes+=mesh;
//   ctypes+="_";
//   ctypes+=drift;
//   ctypes+="_";
//   ctypes+=distance;  
//   ctypes+="_";
//   ctypes+=cthreshold;
//   //ctypes+="-";
//   //ctypes+=ftype;
  
  /// open root tree file
  //  sprintf(fname,"%s%s-Vm%s-Vd%s.root",dirname,rtype,mesh,drift);
  
  const int ARRAYSIZE = tree->GetMaximum("sumpoints")+10;
  cout <<"Array size ="<< ARRAYSIZE<<endl;
  
//   int spoints = 0;
//   int spoints2=0,nsegments=0;
  int maxpoints = 0;
  int evNo;

   
///reading the trakcing information stored in the datatree
  int trackOK, eventTracks; 
  double chi2track[200]; 
  double disttonextcluster[200][6];
  double totchargenextcluster[200][6];

  double Dt=0.;	
  double t0[4]={0.,0.,0.,0.};
  unsigned long long int epoch;
  unsigned long long int nn;

   
  int itrigger;//time of the trigger
  double ttrig; // time of the trigger in seconds 
  char date[4][50];

/// the 4 waveforms
  double *amplC[4];
  for (int i=0; i<4; i++)
    amplC[i] = new double[ARRAYSIZE]; 

///  a smoothed waveform  
  double *samplC;
//   for (int i=0; i<2; i++)
    samplC = new double[ARRAYSIZE]; 
///  a waveform after integration + differentiation    
  double *idamplC;
//for (int i=0; i<4; i++)
    idamplC = new double[ARRAYSIZE]; 
///////////////////////////////////////////////////
//For this error : I added the new iamplC
 //   ././AnalyseTreePicosec.C:1193:48: error: use of undeclared identifier 'iamplC'
   //     iwaveform = new TGraph(maxpoints,ptime,iamplC);
///////////////////////////////////////////////////    
      double *iamplC;
  //for (int i=0; i<4; i++)
    iamplC = new double[ARRAYSIZE]; 
  
  int fitstatus1[4], fitstatus2[4];
  double frmax[4];
  double frmin[4];

  double refP1[200][3];
  double refP2[200][3];

  double bslC[4];
  double rmsC[4];

  double hitX_C[4][200];
  double hitY_C[4][200];
  
cout<<"________________++_________________" << endl;

  double *ptime;
  ptime = new double[ARRAYSIZE];

  double *dampl[4];
  for(int i=0; i<4; i++)
    dampl[i] = new double[ARRAYSIZE];   /// differentiation

/// temporary arrays for the pulse processing
  double *dsampl, *idampl, *iampl, *sampl;
  sampl = new double[ARRAYSIZE];   /// smoothed signal
  dsampl = new double[ARRAYSIZE];  /// differentiation of smoothed array
  idampl = new double[ARRAYSIZE];  /// integration+differentiation
  iampl = new double[ARRAYSIZE];   /// integration with a fixed time constant

 cout<<"__________________________________" << endl;

  ///reading osciloscope parameters
  
  branch = tree->GetBranch("sumpoints");
  branch->SetAddress(&maxpoints);
  branch = tree->GetBranch("eventNo");
  branch->SetAddress(&evNo);
  branch = tree->GetBranch("srsNo");
  branch->SetAddress(&srsNo);
  branch = tree->GetBranch("dt");
  branch->SetAddress(&dt);
  branch = tree->GetBranch("epoch");
  branch->SetAddress(&epoch);
  branch = tree->GetBranch("nn");
  branch->SetAddress(&nn);
  branch = tree->GetBranch("t0");
  branch->SetAddress(t0);

  branch = tree->GetBranch("ttrig");
  branch->SetAddress(&ttrig);
  branch = tree->GetBranch("itrigger");
  branch->SetAddress(&itrigger); 

  branch = tree->GetBranch("fitstatus1");
  branch->SetAddress(fitstatus1);
  branch = tree->GetBranch("fitstatus2");
  branch->SetAddress(fitstatus2);

  branch = tree->GetBranch("rmax");
  branch->SetAddress(frmax);
  branch = tree->GetBranch("rmin");
  branch->SetAddress(frmin);
  ///reading the tracking parameters
  branch = tree->GetBranch("trackOK");
  branch->SetAddress(&trackOK);
  branch = tree->GetBranch("eventTracks");
  branch->SetAddress(&eventTracks);
  branch = tree->GetBranch("chi2track");
  branch->SetAddress(chi2track);

  branch = tree->GetBranch("refP1");
  branch->SetAddress(&refP1[0][0]);
  
  branch = tree->GetBranch("refP2");
  branch->SetAddress(&refP2[0][0]);

  branch = tree->GetBranch("disttonextcluster");
  branch->SetAddress(&disttonextcluster[0][0]);
 

  branch = tree->GetBranch("bslC");
  branch->SetAddress(bslC);
  branch = tree->GetBranch("rmsC");
  branch->SetAddress(rmsC);
  
  int active[4]={0,0,0,0};
  int actch=0;

  for(int i=0; i<4; i++)
  {
      TString channel = TString::Itoa(i+1,10);
      TString bname = "hitX_C"+channel;
      TString btype = bname+"[eventTracks]/D";
      branch = tree->GetBranch(bname);
      if (branch==NULL) 
          continue;
      
      active[i]=1;
      actch++;
      branch->SetAddress(hitX_C[i]);
      
      bname = "hitY_C"+channel;
      btype = bname+"[eventTracks]/D";
      branch = tree->GetBranch(bname);
      branch->SetAddress(hitY_C[i]);

      bname = "amplC"+channel;
      btype = bname+"[sumpoints]/D";
      branch = tree->GetBranch(bname);
      branch->SetAddress(&amplC[i][0]);

  }
  //branch = tree->GetBranch("hitX_C")
  //branch = tree->GetBranch("amplSum");
  //branch->SetAddress(amplSum);
  cout<<RED<<"ACTIVE CHANNELS " <<actch<<endlr;
  for(int i=0; i<4; i++) cout<<(active[i]==1?GREEN:RED)<<"Channel "<<i<<" is " << active[i]<<endlr;

	map<int, RMSBaselineCalculator> rmsBaselineCalculators; //calculators per channel
	for (int ci=0;ci<4;ci++)
	{
		if (!active[ci]) continue;
		RMSBaselineCalculator calc_ci(fname, "RawDataTree", ci, 80, 2, 20, INTEGRATION_TIME_TRIG);
		calc_ci.Process();
		rmsBaselineCalculators[ci] = calc_ci;
	}

  
  int nevents = tree->GetEntries();

	// for (int event_i=0; event_i < nevents; event_i++){
	// 	tree->GetEntry(event_i);
	// 	cout << "Epoch: " << epoch << endl;
	// 	for (int ci=0; ci < 4; ci++) {
	// 		cout << "Channel " << ci + 1 << " get_epoch_rms: " << rmsBaselineCalculators[ci].get_epoch_rms(epoch)<<endl;
	// 	}
	// }
	// return 11;
//   if (runNo==224) nevents-=100;
  cout <<"Found "<<nevents<<" events for the tree "<<endl;  

  double bsl=0.; 
  
  char tmpdir[500];
  sprintf(tmpdir,"%s",gSystem->pwd());
  cout<<"Actual directory: "<<tmpdir<<endl; 

  char plotdirname[500];
  sprintf(plotdirname,"%s/Run%03d/Pool%d",PLOTDIR,runNo,poolNo);//,abs(threshold));
  gSystem->mkdir(plotdirname,kTRUE);
  gSystem->ChangeDirectory(plotdirname); 
  gSystem->ChangeDirectory(tmpdir);


  char allplotdirname[500];
  sprintf(allplotdirname,"%s/moreplots/",plotdirname);//,abs(threshold));
  cout<<MAGENTA<<"more plots at: "<<allplotdirname<<endlr;
  gSystem->mkdir(allplotdirname,kTRUE);
  gSystem->ChangeDirectory(allplotdirname); 
  
  ///  prepare time array 
  ptime[0]=0;
  tree->GetEntry(1);

  for (int i=1;i<ARRAYSIZE;i++)
  {
    ptime[i]=ptime[i-1]+dt;//dt in ns--> ptime in ns...
  }

  cout<<" dt = "<<dt<<endl;

  const double ctoaf = 789.3;   /// emperical value to convert the value of the threshold in amplitude to threshold in charge.
  double peThCh = fabs (ctoaf * peTh);
  peThCh = fabs (ctoaf * peTh);
  cout<<"Analysis thresholds for \"good peak\":   Amplitude = "<< peTh*1000. <<" mV ,  Charge = "<<peThCh <<" pC \n"<<endl;

/// follows a simple event display  - activate with draw = 1 in the input command
/// graphs to draw ecent waveforms and derivatives / integrals
  TGraph* waveform;
//  TGraph* waveform2;
  TGraph* derivative;
  TGraph* derivative2; 
  TGraph* integralh;
  // TGraph* integralh2;
  TGraph* sig_waveform;
  TGraph* sig_d_waveform;
/// canvades for the event display
  char cname[100]; 
  TCanvas *evdcanv[4];
  TCanvas *fitcanv[4];

  if (draw)
  {
      for (int i=0; i<4; i++)
      {
        TString channel = TString::Itoa(i+1,10);        
        if (active[i]==1 && activeDraw[i])
        {
            evdcanv[i] = new TCanvas("EventDisplayC"+channel,"Event display C"+channel,800*2,600*3);
            //evdcanv[i]->SetGrid(1,1);
            evdcanv[i]->Divide(1,3);
            fitcanv[i] = new TCanvas("FitDisplayC"+channel,"Fit display C"+channel,800*2,600*3);
            //fitcanv[i]->SetGrid(1,1);
            fitcanv[i]->Divide(1,3);
        }

      }
//     ecanv = new TCanvas("EventDisplay","Event display");
//     dcanv = new TCanvas("DerivativeDisplay","Derivative display");
//     icanv = new TCanvas("IntegralDisplay","Integral display");
  }
  int eventNo=0; 
  
  int clr[40]={kRed,kBlue,kGreen,kMagenta};
  int col0=1;
  for (int i=4;i<40;i+=4)
  {
    clr[i]=kRed+1*col0;
    clr[i+1]=kBlue+1*col0;
    clr[i+2]=kGreen+1*col0;
    clr[i+3]=kMagenta+1*col0;
    col0++;
  }

/// create output tree    

const int MAXTRIG=100; //maximum number of triggers per channel, i.e. npeaks
//   PICPARAM *ppar, *ppars[4], *spar;
//  PEAKPARAM *ppar, *ppars[4], *spar[4];
  PEAKPARAM *ppar, *spar[4];
  TClonesArray *sparArr[4];
  
   
  for (int i=0; i<4; i++)
    {
      spar[i] = new PEAKPARAM[MAXTRIG];
      sparArr[i] = new TClonesArray("PEAKPARAM",MAXTRIG);
    }
  ppar = new PEAKPARAM;

///parameters for quick monitor of the analysed run
  double T10;// = new double[20000];
  double T90;// = new double[20000];
  double TB10;// = new double[20000];
//   double Tstart;// = new double[20000];
//   double Tend;// = new double[20000];
  double Ampl;// = new double[20000];
  double Charge;// = new double[20000];
  double Width;// = new double[20000];
  double TOT;// = new double[20000];
  double BSLch;// = new double[20000];


//   sprintf(ofname,"%s/Run%03d-Pool%d-%s_treeParam_%.1fmV.root",anadirname,rtype,fabs(threshold*1000.));
  sprintf(ofname,"%s/Run%03d-Pool%d_%streeParam.root",anadirname,runNo,poolNo,rtype);
  cout<<RED<<"Output filename: "<<ofname<<endlr;

  TFile *ofile;
  if (draw)
  {
    sprintf(ofname,"%s/Run%03d-Pool%d_%stest_treeParam.root",anadirname,runNo,poolNo,rtype);
    ofile = new TFile(ofname,"RECREATE");
  }
  else
    ofile = new TFile(ofname,"RECREATE");
  
  sprintf(treename,"ParameterTree");
  sprintf(treetitle,"Pulse Parameter tree");
  
  TTree *otree = new TTree(treename,treetitle);
  
  long double tnow = 0., tlast[]={0.,0.,0.,0.}, tlastgoodpeaks[]={0.,0.,0.,0.};
  long double dtlast[]={0.,0.,0.,0.}, dtlastgoodpeaks[]={0.,0.,0.,0.};
  int npeaks[4] = {0, 0, 0, 0}, ngoodPeaks[4] = {0, 0, 0, 0};
  int ntrigsTot[4]={0,0,0,0};
  int ngoodTrigsTot[4]={0,0,0,0}; 

  otree->Branch("eventNo", &eventNo, "eventNo/I");
  otree->Branch("evtime", &tnow, "evtime/l");

  for (int i=0;i<4;i++)
  {
      if (!active[i]) continue;

      TString channel1 = TString::Itoa(i+1,10);
      TString bname1 = "npeaks_C"+channel1;
      TString btype1 = bname1+"/I";
      otree->Branch(bname1, &npeaks[i], btype1);

//       for(int j=0;j<MAXTRIG;j++)
	  {
		  TString channel = TString::Itoa(i+1,10);
//           TString trigger_num = TString::Itoa(1,10);
		  TString bname = "peakparam_C"+channel;
          //TString btype = bname+"["+bname1+"]";
		  otree->Branch(bname, &sparArr[i]);
          // This creates branches named peakparam_C1_trigger_num_0, peakparam_C1_trigger_num_1, etc.
          // The branch is of type PEAKPARAM, which is a struct defined in MyFunctions.C
	  }
//       TString channel = TString::Itoa(i+1,10);
//       TString bname = "peakparam_C"+channel;
//       TString btype = bname+"["+bname1+"]/D";
//       otree->Branch(bname, spar[i], btype);
  }


  const double microsec = 1e-3;  ///values are in ns by default. Multiply any number by this to convert to microsec
  const double msec = 1e-6;
  
  int nbins = 128;
  double amplMax[4]={0.16};
  double amplMin[4]={0.};
  double chMax[4] = {50.0};
  double Thresholds[4];

  for(int i=0;i<4;i++)
  {
    Thresholds[i]=threshold;
    double rstep = (frmax[i]-frmin[i])/256.;
    cout<<GREEN<<"Rmin["<<i<<"] = "<<frmin[i]<<"  , Rmax = "<<frmax[i]<<"  step = "<<rstep*mV<<" mV"<<endl;
	if (fabs(threshold)<3*rstep)
    {
        cout << "setting new threshold from "<< Thresholds[i] *mV <<" mV   to    ";
        Thresholds[i] = -floor(10000 * 3* rstep)/10000.;
        cout<<Thresholds[i]*mV<<" mV"<<endl;
//         Thresholds[i]=threshold;
    }
	amplMax[i]=(rstep*mV*256);
  }

  double pwMax[4] = {40.,40.,40.,40.};
  double rtMax[4] = {12.,12.,12.,12.};
  double totMax[4]= {40.,40,40,40};

  double rbins = 20;
  double rmax = 20.;
  int dtbins = 800;
  double tmax = 0.002;  //  [ms]
  double single_point_bkg_rejection_probability = pow(total_bkg_rejection_probability, 1.0 / maxpoints );
  double single_point_bkg_rejection_sigmas = TMath::NormQuantile(1 - single_point_bkg_rejection_probability);
#ifdef DEBUGMSG
	cout <<RED << "single_point_bkg_rejection_probability = " << single_point_bkg_rejection_probability << endlr;
#endif
/// vvv This for loop makes threshold +0.02!!!
  for(int i=0;i<4;i++)
  {
    if (oscsetup->AmplifierNo[i]==1 || oscsetup->AmplifierNo[i]==4)
    {
        rtMax[i]=0.6;
        pwMax[i]=20;
        totMax[i]=350;
        chMax[i]=amplMax[i]/mV*(12.5 * ConvFactorCERN) ;
        Thresholds[i] = 20/mV;
    }
    else if (oscsetup->AmplifierNo[i]==3)
    {
        rtMax[i]=180*3;
        pwMax[i]=600*3;
        totMax[i]=600*3;
        chMax[i]=120*3;
        Thresholds[i] = 3/mV;
    }
    else if (oscsetup->AmplifierNo[i]==2)
    {
        rtMax[i]=180*6;
        pwMax[i]=600*6;
        totMax[i]=600*6;
        chMax[i]=120*20;
        Thresholds[i] = 20/mV;
    }
//     else if (oscsetup->AmplifierNo[i]==4)
//     {
//         rtMax[i]=5;
//         pwMax[i]=600*6;
//         totMax[i]=600*6;
//         chMax[i]=120*20;
//         Thresholds[i] = 20/mV;
//     }
    else
      Thresholds[i] = 20/mV;
    
    if (Thresholds[i]>0)
  		Thresholds[i]*=-1;

    if(strncmp(oscsetup->DetName[i],"MCP",3)==0)
    {
        Thresholds[i] = 5/mV;
        rtMax[i]=1.;
        pwMax[i]=20.;
        totMax[i]=10.;
        chMax[i]=amplMax[i]/mV;
    }
    
  }
  /// ^^^ This for loop makes threshold +0.02!!!

  for (int i=0; i<4; i++)
	Thresholds[i] = rmsC[i] * single_point_bkg_rejection_sigmas; // so that its always negative

  for(int i=0;i<4;i++)
  {

    double rstep = (frmax[i]-frmin[i])/256.;
    cout<<GREEN<<"Rmin["<<i<<"] = "<<frmin[i]<<"  , Rmax = "<<frmax[i]<<"  step = "<<rstep*mV<<" mV"<<endl;
    cout<<"Treshold = "<<Thresholds[i]*mV<<" mV"<<endl;
    if (fabs(Thresholds[i])<3*rstep)
    {
        cout << "setting new threshold from "<< Thresholds[i] *mV <<" mV   to    ";
        Thresholds[i] = -floor(10000 * 3* rstep)/10000.;
        cout<<Thresholds[i]*mV<<" mV"<<endl;
//         Thresholds[i]=threshold;
    }
  }

/// calculation of duration 
  tree->GetEntry(0);
  long double epochS = 1.*epoch;
  double framesize = maxpoints * dt * microsec;  //microsec
  tree->GetEntry(nevents-1);
  // tree->GetEntry(nevents-1);
  long double epochF = (1.*epoch + nn*1e-9)+0.004;
  cout<< "Epoch S = "<< epochS<<endl;
  cout<< "Epoch F = "<< epochF<<endl;

  double exprate = nevents/(epochF-epochS) * 75.;  // 75 is the number of peaks per waveform;


  if (!isCalibration)
  {
     exprate = 800 / 8. ; /// 800 events per 8 sec spill
  }
  
    
  int expdt=1;
  if (maxpoints < 5000)
  {
    expdt = ((int) (1./exprate))/1 + 1;   /// Normalizing to seconds for the slow detector 
    cout<<MAGENTA<<"Expected average rate = "<< exprate<<" / sec ,  average DT = "<<1./exprate<<" , implemented DT = "<<expdt<<endlr;
  }


 
  int ntim = 10;                         /// in case needing the rates per 10 or 20 or ... sec make it 10, 20, ...
  double timebinwidth = ntim * 1;   /// multiples of 10 (or 20...) in sec!!!
  
  if (!isCalibration)
  {
//    timebinwidth = 1;  ///1 sec. MAke it 8 to have ~events per spill, or 48 sec to we have events per supercycle
    timebinwidth =  1 * ( ((int)(epochF-epochS)) / 3600 +1);  
 
  }
 

  int tbins = (int) ((epochF-epochS) / (timebinwidth) );   
  
  tbins+=1;
  epochF = epochS + (tbins) * timebinwidth ;
  
  
   double period = timebinwidth / 1.0;
//   int tbins2 = (int)((epochF-epochS)/(period));//

   cout<<RED<<"timebinwidth = "<<timebinwidth << endlr;
   cout<<"__________________________________" << endl;
   


  tree->GetEntry(nevents-1);
  int longpulse = 1; 

  
  if (nevents<5000)
    nbins=64;
  if (nevents>40000)
    nbins=256;
  
  
  if (maxpoints < 11000)
  {
    tmax = (1./exprate ) * 6.;
    dtbins = 100; 
    longpulse = 0;
    rmax = (int) (5*exprate*period);
    if (runNo==219) rmax*=1;
    rbins = rmax; 
     if (rmax>100) rbins = 100;
     if (rmax>1000) rbins = 500;
//     cout<<"rmax = "<<rmax <<"    tmax = "<<tmax<<endl; return 1;
  }
  
  if (nevents<10000)
      dtbins=200;
  if (nevents<5000)
      dtbins=100;
  if (nevents<1000)
      dtbins=50;
  if (nevents<500)
      dtbins=25;
  
  if (!isCalibration)
  {
    period = timebinwidth / 1.0;
    rmax = (200*timebinwidth);  /// suposing maximum rate for the DAQ 200 events / sec.
    int itimebinwidth = int (timebinwidth);
    rbins = 200 / 4 ;
   
  }

  cout<<RED<<"Run duration = "<< epochF - epochS<<" sec.   time bin = "<<timebinwidth <<"  ,  nbins = "<<tbins<<"  ,  rmax = "<<rmax <<endlr;
  
  cout<<RED<<"Expected average rate = "<< exprate<<" / sec ,  average DT = "<<1./exprate<<" , implemented DT = "<<tmax/dtbins<<"  Tmax = "<<tmax<<endlr;
  
  
  ///////// Histogram definitions here /////////////////////////
  
  TH1D *hAMPL[4];
  TH1D *htotCH[4];
//   TH1D *hsCH[4];
  TH1D *hRT[4];
  TH1D *hPW[4];
  TH1D *hPD[4];
  TH1D *hDt[4];
  TH1D *hCutsAMPL[4];
  TH1D *hCutsCH[4];
  TH1D *hCutsDt[4];

  TH1D *hDAMPL[4];
  TH1D *heCH[4];
  TH1D *heCHfit[4];
  TH1D *heCHfixed[4];
  TH1D *htotCHfixed[4];
  TH1D *hionCH[4];


  //TH1D *hCutsRT[4];
  //TH1D *hCutsPW[4];
  //TH1D *hCutsTOT[4];
  TH1D *hRateCuts[4];
  TH1D *htotCHovAmpl[4];
  TH1D *htotCHovAmplCuts[4];
  TH1D *hAmplovDampl[4];
  TH1D *heCHovAmpl[4];
  TH1D *heCHovtotCh[4];
  
  TH1D *hRate[4];
  TH1D *hRateStructure[4];
  TH1D *hSRStriggerStruct[4];
  TH1D *hRateEvolution[4];
  TH1D *hRateEvolutionCheck[4];
  //TH1D *hCutsRateStructure[4];
  TH1D *hRateEvolutionCuts[4];
  //TH1D *hCutsRateStructTrigger[4];//time corrected by trigger position
  //TH1D *hBkgGoodPeaksRateEvolution[4];
//   TH1D *hSparkEvolution[4];
  //   TH1D *hPD;
  
  TH2D *h2deCHvsiCH[4];
  TH2D *h2dtotCHvsAmpl[4];
  TH2D *h2dAmplvsPD[4];
  TH2D *h2dAmplvsRT[4];
  TH2D *h2dAmplvsDampl[4];
  TH2D *h2dCutsCHevolution[4];
  TH2D *h2dCutsAmplEvolution[4];
  TH2D *h2deCHvsAmpl[4];
  TH2D *h2deCHvstotCh[4];
  
  TH2D *h2dXY[4];
  TH2D *h2dXYcuts[4];


//   TH1D* hSparkEvolution;
//cin.get();

  char hname[200], htitle[200], axtitle[200];
  
  for (int i=0;i<4; i++)
  {
    if (!active[i]) continue;
    cout<<"Defining histograms for channel "<<i+1<<endl;
///________________________________________________________________
/// the following plots concern the neutron rates (per pulse)

// Nphotons per pulse if > 0, as found by peak detection
    sprintf(hname,"Run%03d_pool%d_C%d_Rate_%s",runNo,poolNo,i+1,ftype);   /// This is used only when long pulse????
    sprintf(htitle,"C%d rate",i+1);    
    hRate[i]=new TH1D(hname,htitle,rbins,0.,rmax);

// // Nphotons per pulse as found by rate evolution plot!!!
    sprintf(htitle,"Rate (events per %gs) from evolution plot",period);
    sprintf(hname,"Run%03d_pool%d_C%d_Rate_Cuts_per%3.1fseconds_%s",runNo,poolNo,i+1,period,ftype);
    hRateCuts[i]=new TH1D(hname,htitle,rbins,0.,rmax);
    sprintf(axtitle,"average peak rate [events / %gsec]",period);
    hRateCuts[i]->GetXaxis()->SetTitle(axtitle);
//     cout<<RED<<"Histo rate "<<i+1<<" ready"<<endlr;
///________________________________________________________________
/// the following plots concern the neutron (and gamma) rate evolution 
/// (per pulse periode multiple, as it has been calculated from the timebinwidth) 

// all particle average rate per pulse . Must exclude spark & recovery events    
    sprintf(hname,"Run%03d_pool%d_C%d_RateEvolutionAll_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"All events");
    hRateEvolution[i]=new TH1D(hname,htitle,tbins,epochS,epochF);
    hRateEvolution[i]->GetYaxis()->SetTitle("#LT events #GT / sec");
    hRateEvolution[i]->GetYaxis()->SetLabelSize(0.03);
    hRateEvolution[i]->GetYaxis()->SetTitleSize(0.04);
    hRateEvolution[i]->GetYaxis()->SetTitleOffset(1.3);
    hRateEvolution[i]->GetXaxis()->SetTimeDisplay(1);
    hRateEvolution[i]->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%m/%d/%y}");
    hRateEvolution[i]->GetXaxis()->SetLabelSize(0.03);
    hRateEvolution[i]->GetXaxis()->SetLabelOffset(0.02);
    hRateEvolution[i]->SetMinimum(0);
    
// neutron average rate per pulse . Exclude spark events && pulsedBeam (pulse uncorrelated) goodpeaks    
    sprintf(hname,"Run%03d_pool%d_C%d_RateEvolutiongoodPeaks[ci]_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Rate evolution (cuts)");
    hRateEvolutionCuts[i]=new TH1D(hname,htitle,tbins,epochS,epochF);
    hRateEvolutionCuts[i]->GetYaxis()->SetTitle("#LT events #GT / sec");
    hRateEvolutionCuts[i]->GetYaxis()->SetLabelSize(0.03);
    hRateEvolutionCuts[i]->GetYaxis()->SetTitleSize(0.04);
    hRateEvolutionCuts[i]->GetYaxis()->SetTitleOffset(1.3);
    hRateEvolutionCuts[i]->GetXaxis()->SetTimeDisplay(1);
    hRateEvolutionCuts[i]->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%m/%d/%y}");
    hRateEvolutionCuts[i]->GetXaxis()->SetLabelSize(0.03);
    hRateEvolutionCuts[i]->GetXaxis()->SetLabelOffset(0.02);
    hRateEvolutionCuts[i]->SetMinimum(0);

    timebinwidth = hRateEvolution[i]->GetBinWidth(1);
    cout <<"Time bin = "<< timebinwidth<<endl;

// // sparks or baseline recovery events. Useful for rate calculations also    
//     sprintf(hname,"Run%03d_pool%d_C%d_Spark_Evolution_%s",runNo,poolNo,i+1,ftype);
//     sprintf(htitle,"Sparks or recovery");
//     hSparkEvolution=new TH1D(hname,htitle,tbins,epochS,epochF);
    
 // auxiliary plot for the normalization of the rate per pulse. 
    sprintf(hname,"Run%03d_pool%d_C%d_RateEvolution_Check_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Rate normalization plot C%d",i+1);
    hRateEvolutionCheck[i]=new TH1D(hname,htitle,tbins,epochS,epochF);

///________________________________________________________________
/// the following plots concern the rate structures WITHOUT correction for the beam trigger

    
     //sprintf(hname,"Run%03d_pool%d_C%d_RateStructure_GE_%s",runNo,poolNo,i+1,ftype);
     //sprintf(htitle,"Single p.e.");
     //hCutsRateStructure[i]=new TH1D(hname,htitle,250,0.,framesize);
     //hCutsRateStructure[i]->GetXaxis()->SetTitle("t [#mus]");
     //hCutsRateStructure[i]->SetMinimum(0);

    sprintf(hname,"Run%03d_pool%d_C%d_RateStructure_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Single p.e.");
    hRateStructure[i]=new TH1D(hname,htitle,250,0.,framesize);
    hRateStructure[i]->GetXaxis()->SetTitle("t [#mus]");
//     hRateStructure->SetMinimum(0);
    
    sprintf(hname,"Run%03d_pool%d_C%d_TriggerStructure_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"SRS trigger ");
    hSRStriggerStruct[i]=new TH1D(hname,htitle,250,0.,framesize);
    hSRStriggerStruct[i]->SetLineColor(kGreen+2);

///________________________________________________________________
///  DeltaT plots 

    sprintf(hname,"Run%03d_pool%d_C%d_consecutivePulses_#DeltaT_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Concecutive pulses #DeltaT");
    hDt[i] = new TH1D(hname,htitle,dtbins,0.,tmax);
    hDt[i]->GetXaxis()->SetTitle("#DeltaT [sec]");

    sprintf(hname,"Run%03d_pool%d_C%d_consecutivePulses_#DeltaT_Cuts_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"#DeltaT goodpeaks");
    hCutsDt[i] = new TH1D(hname,htitle,dtbins,0.,tmax);
    hCutsDt[i]->GetXaxis()->SetTitle("#DeltaT [sec]");

///________________________________________________________________
///  Amplitude plots 

    
    sprintf(hname,"Run%03d_pool%d_C%d_Amplitude_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Pulse amplitude C%d",i+1);
    hAMPL[i]=new TH1D(hname,htitle,nbins,0.,amplMax[i]);
    hAMPL[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");

    sprintf(htitle,"Pulse amplitude (after cuts)");
    sprintf(hname,"Run%03d_pool%d_C%d_Amplitude_Cuts_%s",runNo,poolNo,i+1,ftype);
    hCutsAMPL[i]=new TH1D(hname,htitle,nbins,0.,amplMax[i]);
    hCutsAMPL[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");
    
    sprintf(hname,"Run%03d_pool%d_C%d_DataAmplitude_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Pulse amplitude from data C%d",i+1);
    hDAMPL[i]=new TH1D(hname,htitle,nbins,0.,amplMax[i]);
    hDAMPL[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");

    sprintf(htitle,"Pulse amplitude evolution C%d",i+1);
    sprintf(hname,"Run%03d_pool%d_C%d_AmplitudeEvolution_%s",runNo,poolNo,i+1,ftype);
    h2dCutsAmplEvolution[i] = new TH2D(hname,htitle,100,epochS,epochF,nbins/2,0.,amplMax[i]);
    h2dCutsAmplEvolution[i]->GetXaxis()->SetTimeDisplay(1);
    h2dCutsAmplEvolution[i]->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%m/%d/%y}");
    h2dCutsAmplEvolution[i]->GetXaxis()->SetLabelSize(0.03);
    h2dCutsAmplEvolution[i]->GetXaxis()->SetLabelOffset(0.02);
    h2dCutsAmplEvolution[i]->GetYaxis()->SetTitle("Pulse amplitude [mV]");

///________________________________________________________________
///  Charge plots 

    sprintf(hname,"Run%03d_pool%d_C%d_Charge_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Total charge C%d",i+1);
    htotCH[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    htotCH[i]->GetXaxis()->SetTitle("charge [pC]");
    
    sprintf(htitle,"Total charge (after cuts)");
    sprintf(hname,"Run%03d_pool%d_C%d_Charge_Cuts_%s",runNo,poolNo,i+1,ftype);
    hCutsCH[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    hCutsCH[i]->GetXaxis()->SetTitle("charge [pC]");

    sprintf(hname,"Run%03d_pool%d_C%d_Charge_Fixed_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Total charge (fixed duration)C%d",i+1);
    htotCHfixed[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    htotCHfixed[i]->GetXaxis()->SetTitle("charge [pC]");

    sprintf(hname,"Run%03d_pool%d_C%d_eCharge_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"e-peak charge C%d",i+1);
    heCH[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    heCH[i]->GetXaxis()->SetTitle("charge [pC]");

    sprintf(hname,"Run%03d_pool%d_C%d_eChargeFit_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"e-peak charge from fit C%d",i+1);
    heCHfit[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    heCHfit[i]->GetXaxis()->SetTitle("charge [pC]");

    sprintf(hname,"Run%03d_pool%d_C%d_eChargeFixed_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"e-peak charge (fixed duration) C%d",i+1);
    heCHfixed[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    heCHfixed[i]->GetXaxis()->SetTitle("charge [pC]");

    sprintf(hname,"Run%03d_pool%d_C%d_ionCharge_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"Ion tail charge C%d",i+1);
    hionCH[i]=new TH1D(hname,htitle,nbins*1,0.,chMax[i]);
    hionCH[i]->GetXaxis()->SetTitle("charge [pC]");

    sprintf(htitle,"Pulse charge evolution C%d",i+1);
    sprintf(hname,"Run%03d_pool%d_C%d_ChargeEvolution_%s",runNo,poolNo,i+1,ftype);
    h2dCutsCHevolution[i] = new TH2D(hname,htitle,100,epochS,epochF,nbins/2,0.,chMax[i]);
    h2dCutsCHevolution[i]->GetXaxis()->SetTimeDisplay(1);
    h2dCutsCHevolution[i]->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%m/%d/%y}");
    h2dCutsCHevolution[i]->GetXaxis()->SetLabelSize(0.03);
    h2dCutsCHevolution[i]->GetXaxis()->SetLabelOffset(0.02);
    h2dCutsCHevolution[i]->GetYaxis()->SetTitle("charge [pC]");
    
    
///________________________________________________________________    
///  pulse time properties plots     
    
    sprintf(htitle,"Rise Time C%d",i+1);
    sprintf(fname2,"Run%03d_pool%d_C%d_Risetime_%s",runNo,poolNo,i+1,ftype);
    hRT[i]=new TH1D(fname2,htitle,int(50*rtMax[i]/dt),0.,rtMax[i]);
    hRT[i]->GetXaxis()->SetTitle("Risetime [ns]");  
    
    sprintf(htitle,"Pulse Width C%d",i+1);
    sprintf(fname2,"Run%03d_pool%d_C%d_Pulse_Width_%s",runNo,poolNo,i+1,ftype);
    hPW[i]=new TH1D(fname2,htitle,(int)floor(pwMax[i]/dt),0.,pwMax[i]);
    hPW[i]->GetXaxis()->SetTitle("e-Peak Width [ns]");

    sprintf(htitle,"Pulse Duration C%d",i+1);
    sprintf(fname2,"Run%03d_pool%d_C%d_Pulse_Duration_%s",runNo,poolNo,i+1,ftype);
    hPD[i]=new TH1D(fname2,htitle,(int)floor(totMax[i]/dt),0.,totMax[i]);
    hPD[i]->GetXaxis()->SetTitle("Pulse Duration [ns]");
///________________________________________________________________
///  Correlation plots    
    double chovamax = 500.;
    double amplovdampl = 500.;
    if (dv[i]<=50 || vm[i]>490) chovamax = 800.;
    
    if (oscsetup->AmplifierNo[i]==1 || oscsetup->AmplifierNo[i]==4)
    {
       chovamax =20;   
       amplovdampl =10.;
    }
    if (oscsetup->AmplifierNo[i]==3)
    {
       chovamax =4000;   
       amplovdampl =20.;
    }
    if (oscsetup->AmplifierNo[i]==2)
    {
       chovamax =30000;   
       amplovdampl =20.;
    }
    else
    {
       chovamax =30000;   
       amplovdampl =20.;
    }

    chovamax = (htotCH[i]->GetXaxis()->GetXmax())/(hAMPL[i]->GetXaxis()->GetXmax()) *2.;
    
    sprintf(htitle,"Total Charge over Amplitude");
    sprintf(hname,"Run%03d_pool%d_C%d_Charge_ov_Amplitude_%s",runNo,poolNo,i+1,ftype);
    htotCHovAmpl[i]=new TH1D(hname,htitle,nbins*2,0.,2*chovamax);
    htotCHovAmpl[i]->GetXaxis()->SetTitle("ratio [pC/mV]");

    sprintf(htitle,"Total Charge over Amplitude (after cuts)");
    sprintf(hname,"Run%03d_pool%d_C%d_selected_Charge_ov_Amplitude_%s",runNo,poolNo,i+1,ftype);
    htotCHovAmplCuts[i]=new TH1D(hname,htitle,nbins*2,0.,2*chovamax);
    htotCHovAmplCuts[i]->SetLineColor(kRed);
    htotCHovAmplCuts[i]->GetXaxis()->SetTitle("ratio [pC/mV]");

  	sprintf(htitle,"e-peak Charge over Amplitude");
  	sprintf(hname,"Run%03d_pool%d_C%d_eCharge_ov_Amplitude_%s",runNo,poolNo,i+1,ftype);
  	heCHovAmpl[i]=new TH1D(hname,htitle,nbins*2,0.,chovamax);
  	heCHovAmpl[i]->GetXaxis()->SetTitle("ratio [pC/mV]");

  	sprintf(htitle,"e-peak Amplitude (max data) over e-epeak Amplitude Fit");
  	sprintf(hname,"Run%03d_pool%d_C%d_selected_Amplitude_ov_Damplitude_%s",runNo,poolNo,i+1,ftype);
  	// hAmplovDampl[i]=new TH1D(hname,htitle,nbins*2,0.,amplMax[i]);
  	hAmplovDampl[i]=new TH1D(hname,htitle,nbins*2,0.,2.);
  	hAmplovDampl[i]->SetLineColor(kBlack);
  	hAmplovDampl[i]->GetXaxis()->SetTitle("ratio [V/V]");

    sprintf(hname,"Run%03d_pool%d_Ch%d_eCh_ov_iCh_%s",runNo,poolNo,i+1,ftype);
    sprintf(htitle,"e-peak charge / total charge C%d",i+1);
  	heCHovtotCh[i]=new TH1D(hname,htitle,nbins*2,0.,1.);
  	heCHovtotCh[i]->SetLineColor(kBlue);
  	heCHovtotCh[i]->GetXaxis()->SetTitle("ratio [pC/pC]");


    sprintf(fname2,"Run%03d_pool%d_C%d_Amplitude_vs_dataAmplitude_C_%s",runNo,poolNo,i+1,ftype);
    h2dAmplvsDampl[i]=new TH2D(fname2,fname2,nbins/2,0.,amplMax[i],nbins/2,0.,amplMax[i]);
    h2dAmplvsDampl[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");
    h2dAmplvsDampl[i]->GetYaxis()->SetTitle("Pulse amplitude [mV]");
    h2dAmplvsDampl[i]->SetStats(0);

    sprintf(fname2,"Run%03d_pool%d_C%d_Charge_vs_Amplitude_C_%s",runNo,poolNo,i+1,ftype);
    h2dtotCHvsAmpl[i]=new TH2D(fname2,fname2,nbins,0.,amplMax[i],nbins,0.,chMax[i]);
    h2dtotCHvsAmpl[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");
    h2dtotCHvsAmpl[i]->GetYaxis()->SetTitle("Total charge [pC]");
    h2dtotCHvsAmpl[i]->SetStats(0);

    sprintf(fname2,"Run%03d_pool%d_C%d_eCharge_vs_Amplitude_C_%s",runNo,poolNo,i+1,ftype);
    h2deCHvsAmpl[i]=new TH2D(fname2,fname2,nbins,0.,amplMax[i],nbins,0.,chMax[i]);
    h2deCHvsAmpl[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");
    h2deCHvsAmpl[i]->GetYaxis()->SetTitle("e-peak charge [pC]");
    h2deCHvsAmpl[i]->SetStats(0);

    sprintf(fname2,"Run%03d_pool%d_C%d_eCharge_vs_totCharge_C_%s",runNo,poolNo,i+1,ftype);
    h2deCHvstotCh[i]=new TH2D(fname2,fname2,nbins,0.,chMax[i],nbins,0.,chMax[i]);
    h2deCHvstotCh[i]->GetXaxis()->SetTitle("Total charge [pC]");
    h2deCHvstotCh[i]->GetYaxis()->SetTitle("e-peak charge [pC]");
    h2deCHvstotCh[i]->SetStats(0);

    sprintf(fname2,"Run%03d_pool%d_C%d_Ampl_vs_Rise_Time_C_%s",runNo,poolNo,i+1,ftype);
    h2dAmplvsRT[i]=new TH2D(fname2,fname2,nbins,0.,amplMax[i],(int)(50*rtMax[i]),0.,rtMax[i]);
    h2dAmplvsRT[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");
    h2dAmplvsRT[i]->GetYaxis()->SetTitle("Risetime [ns]");
    h2dAmplvsRT[i]->SetStats(0);
    
    sprintf(fname2,"Run%03d_pool%d_C%d_Amplitude_vs_PD_C_%s",runNo,poolNo,i+1,ftype);
    h2dAmplvsPD[i]=new TH2D(fname2,fname2,nbins,0.,amplMax[i],int(totMax[i]/dt),0.,totMax[i]);
    h2dAmplvsPD[i]->GetXaxis()->SetTitle("Pulse amplitude [mV]");
    h2dAmplvsPD[i]->GetYaxis()->SetTitle("Pulse Duration [ns]");
    h2dAmplvsPD[i]->SetStats(0);

    sprintf(fname2,"Run%03d_pool%d_C%d_eCharge_vs_iCharge_C_%s",runNo,poolNo,i+1,ftype);
    h2deCHvsiCH[i]=new TH2D(fname2,fname2,nbins,0.,chMax[i],nbins,0.,chMax[i]);
    h2deCHvsiCH[i]->GetXaxis()->SetTitle("e-Charge [pC]");
    h2deCHvsiCH[i]->GetYaxis()->SetTitle("i-Charge [pC]");
    h2deCHvsiCH[i]->SetStats(0);

    sprintf(htitle,"Hit map C%d",i+1);
    sprintf(fname2,"Run%03d_pool%d_C%d_XYmap_C_%s",runNo,poolNo,i+1,ftype);
    h2dXY[i]=new TH2D(fname2,htitle,200,-100.,100.,200,-100.,100.);
    h2dXY[i]->GetXaxis()->SetTitle("X [mm]");
    h2dXY[i]->GetYaxis()->SetTitle("Y [mm]");
    h2dXY[i]->SetStats(0);

    sprintf(htitle,"Hit map C%d",i+1);
    sprintf(fname2,"Run%03d_pool%d_C%d_XYmap_Cuts_C_%s",runNo,poolNo,i+1,ftype);
    h2dXYcuts[i]=new TH2D(fname2,htitle,200,-100.,100.,200,-100.,100.);
    h2dXYcuts[i]->GetXaxis()->SetTitle("X [mm]");
    h2dXYcuts[i]->GetYaxis()->SetTitle("Y [mm]");
    h2dXYcuts[i]->SetStats(0);
  }
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// END the 4-loop of histograms here!
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int npt = 2;
  double DT = 0.5;  /// default derivation time in [ns]
  double DTI = 6.0;  /// default integration time in [ns]
  int evpm =1000;
  if (nevents >100000) evpm = 10000;
  if (nevents >50000) evpm = 5000;
  if (maxpoints >50000) evpm = 10;

  eventNo = 0;
  int totNtrigs=0;
  int totNsparks = 0;
  long double drawdt=0.;
//   if (oscsetup->AmplifierNo[i]==5) eventNo=2000; ///skip first events because of the time change in the oscilloscope.
//   nevents=10000;
  cout <<"Start processing the "<<nevents<<" events"<<endl;  
  //

  int successfulFits_sigmoid = 0;
  int totalFits_sigmoid = 0;
  int successfulFits_double_sigmoid = 0;
  int totalFits_double_sigmoid = 0;

  int total_secondary_count = 0;
  int total_thin_count = 0;
	while (eventNo<nevents)
  {
//	if (eventNo<14255 || eventNo>14257) { eventNo++; continue;};
  	//cout << "Event Number: " << eventNo << endl;
  	// if (eventNo < 1200) { eventNo++; continue; }
  	//if (eventNo!=47) { eventNo++; continue;}
#ifdef DEBUGMSG
    cout << "Event number: " << eventNo << endl;
#endif

    if (draw)
	{
	  cout<<endl<<"______________________________________________________________________"<<endl<<endl;;
#ifdef DEBUGMSG
    	//cout<<endl<<"Entering 1st draw_____________________________________________________"<<endl<<endl;;
#endif
    	cout<<KRESET<<"Event to draw : ";
	  cin>>eventNo;
	}

      if (eventNo<0 || eventNo>=nevents) break;
      ///    if (eventNo<25 || eventNo>=150) {eventNo++; continue;}

      tree->GetEntry(eventNo);
  	  //Here we are setting new rms and bsl calculation from the RMSBaselineCalculator class
  	  for (int ci=0; ci<4; ci++) {
  	  	if (!active[ci]) continue;
	  	  rmsC[ci] = rmsBaselineCalculators[ci].get_epoch_rms(epoch);
  	  	//  bslC[ci] = rmsBaselineCalculators[ci].get_epoch_baseline(epoch);  // Imperfect waveform-by-waveform baseline correction done in previous code. Will adjust later with adjust_basline
  	  }

      double maxc[]={0.,0.,0.,0.};
      double minc[]={0.,0.,0.,0.};
      double maxd[]={0.,0.,0.,0.};
      double mind[]={0.,0.,0.,0.};
      long double evtime =  1.*epoch+1.*nn*1e-9;
//       long double evtime = 1.*epoch+1.*nn*1e-9 + t0;
//       printf("t0 = %8.8lf ,  evtime = %8.8LF\n",t0,evtime);

      int ntrigs = 0;
      int ntrigsCuts =0;
      int detspark = 0;
      double totcharge = 0.;

//       cout<<"frmax = "<<frmax<<"   "<<frmin<<"    "<<frmin - frmax<<endl;

/// reset time array in case dt was modified during data taking
  ptime[0]=0;
  for (int i=1;i<maxpoints;i++)
        ptime[i]=ptime[i-1]+dt;


      /////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// START the 4-channel loop here!
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////

  for (int ci = 0; ci < 4; ++ci)
  {
      if (!active[ci]) continue;
      
      if (draw)
        if (!activeDraw[ci]) continue;
        
      peTh = Thresholds[ci];

      ntrigs = 0;
      sparArr[ci]->Clear();
      ntrigsCuts =0;
      detspark = 0;
      totcharge = 0.;
      //cout<< "maxpoints = "<< maxpoints << endl;
      //cout<<"sampling of = dt = "<<dt<<endl;
      for (int i=0;i<maxpoints;i++)
        totcharge+= amplC[ci][i];
//         totcharge+= bslC[ci] - amplC[ci][i];

      totcharge*= -1;   ///pulses are negative


      totcharge /= maxpoints;
      totcharge *= DT*1e3;   /// make it fC

      if (fabs(totcharge)>200. || rmsC[ci]> 1.0015)   ///to be fine-tuned - values from nBLM
      {
        detspark = 1;
        totNsparks++;
      }

      npt = TMath::FloorNint(DT/dt)+1;   /// to be fine-tuned. Amplifier depended!
//  return 0;

    if (draw)
    {
#ifdef DEBUGMSG
        //cout<<endl<<"Entering 2nd if(draw) Fine tune for SmoothArray, DerivateArray and IntegratePulse that will be used for the analysis _________________________"<<endl<<endl;;
#endif
		//cout<<GREEN<<"Baseline rms per channel "<<ci+1<<" = "<<rmsC[ci]<<endl;
    	//cin.get();
        //cout<<endl<<"Event "<< eventNo<<" Channel "<<ci+1<<"\t fit1 "<<fitstatus1[ci]<<" fit2 "<<fitstatus2[ci]<< " bsl "<<bslC[ci]<<" rms "<<rmsC[ci]<< " totcharge "<<totcharge<< endl;
        //cout<<"Pulse length = "<<maxpoints<<endl;
        long double epochX = (1.*epoch + nn*1e-9);
        TTimeStamp *tstamp = new TTimeStamp(epochX+(rootConv-unixConv),0);
        printf("Event time : %s = %10.9Lf  DT = %10.9LF s = %10.7LF ms\n",tstamp->AsString("l"),epochX,epochX-drawdt,1000.*(epochX-drawdt));
        drawdt = epochX;
        cout <<"Default derivation for "<<DT<<" ns , with sampling of "<<dt<<" ns ==> derivation points npt = "<<npt<< endl;


        evdcanv[ci]->cd(1);//gPad->SetGrid(1,1);
        waveform = new TGraph(maxpoints,ptime,amplC[ci]);
        //waveform2 = new TGraph(maxpoints,ptime,amplC[ci]);
        maxc[ci]=TMath::MaxElement(maxpoints,amplC[ci]);
        minc[ci]=TMath::MinElement(maxpoints,amplC[ci]); //the peak amplitude of the pulse
	    sprintf(cname,"Event %d Waveform C%d\n",evNo,ci+1);
        waveform->SetTitle(cname);
        waveform->SetLineColor(clr[ci+1]);
        waveform->SetMarkerColor(clr[ci]);
        waveform->SetFillColor(0);

        double fmin = frmin[ci]-frmax[ci];
        waveform->GetHistogram()->SetMinimum(fmin);
        cout<<"Setting minimum at "<<frmin[ci]-frmax[ci]<<endl;

       waveform->Draw("apl");
       //waveform->Draw("pl");

        waveform->GetHistogram()->GetYaxis()->SetRangeUser(fmin,-fmin/8.);
        evdcanv[ci]->Modified();
        evdcanv[ci]->Update();


///  Smoothing array when no "bit filter" !!!!

      evdcanv[ci]->cd(2); gPad->SetGrid(1,1);
      double maxelement_dampl[4]= {0,0,0,0};
      double minelement_dampl[4] = {0,0,0,0};
      int nsmooth = 3;
    	// cout<<BLUE<< "Smooting starts here for: "<< npt <<" points"<<endlr;
    	// cout<<BLUE<<"                           "<<endlr;
    	// cout<<BLUE<<" data to Sampl"<<endlr;
    	// cout<<YELLOW<<"                           "<<endlr;
    	// cout<<YELLOW<<" Derivation of to Sampl to dampl"<<endlr;
      SmoothArray(amplC[ci], samplC, maxpoints, nsmooth, 1);
      DerivateArray(samplC,dampl[ci],maxpoints,dt,npt,1); ///with the number of points

      maxelement_dampl[ci] = TMath::MaxElement(maxpoints, dampl[ci]);
      minelement_dampl[ci] = TMath::MinElement(maxpoints, dampl[ci]);
      derivative = new TGraph(maxpoints,ptime,dampl[ci]); // the dapl arr has changed
      derivative2 = new TGraph(maxpoints,ptime,dampl[ci]);
      sprintf(cname,"Derivative C%d\n",ci+1);
      derivative->SetMarkerColor(clr[ci+8]);
      derivative->SetLineColor(clr[ci+8]);
      derivative->SetFillColor(0);
      derivative->SetTitle(cname);
      derivative->Draw("apl");
      derivative->SetMaximum(maxelement_dampl[ci]+0.1);
      derivative->SetMinimum(minelement_dampl[ci]-0.1);
      evdcanv[ci]->Modified();
      evdcanv[ci]->Update();
      //continue;
///  Derivate smoothed signals for analysis  (may not be used...)
    //if (oscsetup->AmplifierNo[ci]==1) nsmooth = 3;
// //       SmoothArray(amplC[ci], samplC, maxpoints, nsmooth, 1);
      //cout<<RED<<"Ready to derivate smoothed array"<<endlr;
      // cout<<RED<< "Smooting here for: "<< npt <<" points"<<endlr;
      // cout << RED << "                           " << endlr;
      // cout << RED << " data to Sampl" << endlr;
      // cout << YELLOW << "                           " << endlr;
      // cout << YELLOW << " Derivation of to Sampl to dsampl now" << endlr;;
      DerivateArray(samplC, dsampl,maxpoints,dt,npt,1);
      //cout<<YELLOW<<"Ready to plot Smoothed derivative"<<endlr;
      TGraph *graph22 = new TGraph(maxpoints,ptime,dsampl);
      sprintf(cname,"Smoothed Derivative C%d\n",ci+1);
      graph22->SetMarkerColor(clr[ci+4]);
      graph22->SetLineColor(clr[ci+4]);
      graph22->SetLineWidth(2);
      graph22->SetFillColor(0);
      graph22->SetTitle(cname);
      //graph22->Draw("pl");
      //graph22->Draw("same");
      evdcanv[ci]->Modified();
      evdcanv[ci]->Update();

      evdcanv[ci]->cd(3); gPad->SetGrid(1,1);
      double intgrA = IntegrateA(maxpoints,dsampl,idamplC,dt);
      TGraph *iwaveform = new TGraph(maxpoints,ptime,idamplC);
      sprintf(cname,"Event %d Int-Waveform C%d\n",evNo,ci+1);
      iwaveform->SetTitle(cname);
      iwaveform->SetLineColor(clr[ci+4]);
      iwaveform->SetMarkerColor(clr[ci+4]);
      iwaveform->SetFillColor(0);
      iwaveform->Draw("apl");
      //evdcanv[ci]->Modified();
      //evdcanv[ci]->Update();

  ///Integration of Pulse here

      int nint = 20; //default to change you need to see MyFunctions.h
      nint = N_INTEGRATION_POINTS;
      //double DTI2 = 2.;  ///default
      double DTI2 = 2.; //time integration for 2ns
      nint = TMath::FloorNint(DTI2/dt)+1;
      // cout<<"Integration points = "<<nint<<endl;
      double intgr = IntegratePulse(maxpoints,idamplC,iampl,dt,nint*dt);
     // continue;
      maxd[ci]=TMath::MaxElement(maxpoints,iampl);
      mind[ci]=TMath::MinElement(maxpoints,iampl);
      integralh = new TGraph(maxpoints,ptime,iampl);
      // integralh2 = new TGraph(maxpoints,ptime,iampl);
      sprintf(cname,"Event %d - Integral %g of C%d\n",evNo, nint*dt,ci+1);
      integralh->SetMarkerColor(clr[ci+4]);
      integralh->SetLineColor(clr[ci+4]);
      integralh->SetFillColor(0);
      integralh->SetTitle(cname);
      integralh->SetMaximum(maxelement_dampl[ci]+0.1);
      integralh->SetMinimum(mind[ci]-0.1);
      integralh->Draw("apl");
      evdcanv[ci]->Modified();
      evdcanv[ci]->Update();

      DTI2 = 1.5;
      nint = TMath::FloorNint(DTI2/dt)+1;;
      intgr = IntegratePulse(maxpoints,idamplC,iampl,dt,nint*dt);
      TGraph *graph23 = new TGraph(maxpoints,ptime,iampl);
      sprintf(cname,"Event %d - Integral %g of C%d\n",evNo, nint*dt,ci+1);
      graph23->SetMarkerColor(clr[ci+10]);
      graph23->SetLineColor(clr[ci+10]);
      graph23->SetFillColor(0);
      graph23->SetTitle(cname);
      //graph23->Draw("pl");

      DTI2 = 20.;
      nint = TMath::FloorNint(DTI2/dt)+1;
      intgr = IntegratePulse(maxpoints,idamplC,iampl,dt,nint*dt);
      graph23 = new TGraph(maxpoints,ptime,iampl);
      sprintf(cname,"Event %d - Integral %g of C%d\n",evNo,nint*dt,ci+1);
      graph23->SetMarkerColor(clr[ci+15]);
      graph23->SetLineColor(clr[ci+15]);
      graph23->SetFillColor(0);
      graph23->SetTitle(cname);
      //graph23->Draw("pl");
      evdcanv[ci]->Modified();
      evdcanv[ci]->Update();

      //cout<<RED<<"OLA KALA"<<endlr;
      //continue;

  }///end (if(draw))

    /// smooth sumn array for analysis
    int nsmooth = 1;
  	//int nsmooth = 3;
    if (oscsetup->AmplifierNo[ci]==1 || oscsetup->AmplifierNo[ci] == 3) // smoothing with 3 points
    {
        DT = 2.;
        nsmooth = 3;
        npt = TMath::FloorNint(DT/dt)+1; //2ns derivation window
    }
//     if (oscsetup->AmplifierNo[ci]==1) nsmooth = 15;
#ifdef DEBUGMSG
  	cout<<BLUE<< "Smooting for analysis here : "<< npt <<" points"<<endlr;
  	cout<<BLUE<<"                           "<<endlr;
  	cout<<BLUE<<" data to Sampl"<<endlr;
  	cout<<YELLOW<<"                           "<<endlr;
  	cout<<YELLOW<<" Derivation of to ampl to dsampl"<<endlr;
#endif

    SmoothArray(amplC[ci], sampl, maxpoints, nsmooth, 1);

    /// derivate data array
    DerivateArray(amplC[ci],dsampl,maxpoints,dt,npt,1);
    /// derivate smoothed data array
	//DerivateArray(sampl,dsampl,maxpoints,dt,npt,1);
    //continue;
    //return 11;

    int itrig =0;
    ntrigs = 0;
    ntrigsCuts=0;
    int pulsedBeam = 0;
    itrig = itrigger;
    if (itrig>0)
        pulsedBeam = 1;
        //pulsedBeam = callibrationRun;
/// Having the following "if (detspark)" here means that an event with a spark, of a "baseline recovery event will not be analysed!!!
    if (detspark)
    {
//        hSparkEvolution->Fill(evtime);
       cout<<" Spark ("<< (totcharge>10.) <<") or recovery ("<< (rmsC[ci]> 0.0015) <<") at event "<<eventNo<<endl;   /// correct the values
       eventNo++;
       continue;
    }
  	vector<pair<double, double>> trigger_windows;
  	if (strncmp(oscsetup->DetName[ci], "MM", 2) == 0) {
		adjust_baseline(maxpoints, ptime, sampl);
  		double trigger_threshold = rmsBaselineCalculators[ci].get_epoch_integral_rms(epoch) * single_point_bkg_rejection_sigmas;
		TriggerResult trigger_results = GetTriggerWindows(ptime, maxpoints, sampl, dt, trigger_threshold);
  	  	//cout<< "Event number: " << eventNo << " Channel number: " << ci << endl;
  		//cin.get();
  		//print event number that had secondary pulses
  		if (trigger_results.secondary_count > 0) {
  			// cout<<"Event number: "<<eventNo<<" Secondary pulses detected: "<<trigger_results.secondary_count<<endl;
  		}
  		trigger_windows = trigger_results.pulse_bounds_filtered;
  		total_secondary_count += trigger_results.secondary_count;
  		total_thin_count += trigger_results.thin_count;

  		// cout << "Secondary pulses rejected: " << secondary_counts << endl;
		// cout << "Thin pulses rejected: " << thin_counts << endl;
	}
  	auto trigger_windows_iterator = trigger_windows.begin();

    int ti = 0;
    bool isMM = kFALSE;
    
    while (ti < maxpoints-50 && ntrigs<MAXTRIG)
	{
       ppar->Reset();
	  //cout<<"check: "<<ntrigs+1<<endl;
	  ppar->rms=rmsC[ci];
	  ppar->bsl=bslC[ci];
///*************************************************************

      if(strncmp(oscsetup->DetName[ci], "MCP", 3) == 0 && SIGMOIDFIT[ci] && DOUBLESIGMOIDFIT[ci]==false)  // Checks if "MCP" is at the start of the string
      {
#ifdef DEBUGMSG
        cout<<BLUE<<"Channel "<<ci+1<<" is MCP"<<endlr;

      	cout<<RED<<"Threshold "<<Thresholds[ci]<<endlr;
#endif
        ti = AnalyseLongPulseMCP(maxpoints,evNo,sampl,dt,dsampl,ppar,Thresholds[ci],sig_tshift[ci], ti);

  		//ti = AnalyseLongPulseMCP(maxpoints,evNo,sampl,dt,dsampl,ppar,-0.040,sig_tshift[ci], ti);

    	if (ti<0) break;
		if (ti < maxpoints-50) {
    		successfulFits_sigmoid += ppar->SigmoidfitSuccess;
    		totalFits_sigmoid++;
    	}
      }
      else if (strncmp(oscsetup->DetName[ci], "MM", 2) == 0)
      {
         isMM = kTRUE;
      	// break;
      	//cout<<BLUE<<"Channel "<<ci+1<<" uses the fit. Threshold = "<<Thresholds[ci]*mV<<" mV"<<endlr;
	      //ti = AnalyseLongPulseCiv(maxpoints,evNo,sampl,dsampl,ppar,threshold, dt, ti);    /// all the analysis is done here!!!!
#ifdef DEBUGMSG
      	// Open a file to write the data
      	if(ci == 1 || ci == 3) {
      		double trigger_threshold = rmsBaselineCalculators[ci].get_epoch_integral_rms(epoch) * single_point_bkg_rejection_sigmas;
      		string filename = "waveform_data_" + to_string(ci) + ".txt";
      		ofstream outFile(filename);
      		if (!outFile) {
      			cerr << "Error: Could not open file for writing!" << endl;
      			return 1;
      		}

      		// Write event metadata
      		outFile << "Event data:\n";
      		outFile << "const int points = " << maxpoints << ";\n";
      		outFile << "dt: " << dt << "\n";
      		outFile << "RMS: " << ppar->rms << "\n";
      		outFile << "BSL: " << ppar->bsl << "\n";
      		outFile << "Threshold: " << Thresholds[ci] << "\n";
      		outFile<< "trigger_threshold: " << trigger_threshold << "\n";
      		outFile << "Sig_tshift: " << sig_tshift[ci] << "\n";
      		outFile << "Event number: " << evNo << "\n";
			outFile << "INTEGRATION_TIME_TRIGGER: " << INTEGRATION_TIME_TRIG << "\n";
      		outFile << " total_bkg_rejection_probability: " << total_bkg_rejection_probability << "\n";
      		outFile << " ion_tail_end_point_threshold_fraction: " << ion_tail_end_point_threshold_fraction << "\n";
			outFile << " CIVIDEC_PULSE_DURATION: " << CIVIDEC_PULSE_DURATION << "\n";
      		outFile << " CIVIDEC_PEAK_DURATION: " << CIVIDEC_PEAK_DURATION << "\n";
      		// Write the sampl array
      		outFile << "double data[" << maxpoints << "] = {";
      		for (int i = 0; i < maxpoints; ++i) {
      			outFile << fixed << setprecision(6) << sampl[i];
      			if (i < maxpoints - 1) outFile << ", ";
      		}
      		outFile << "};\n";
			//Write time array
      		outFile << "double time[" << maxpoints << "] = {";
      		for (int i = 0; i < maxpoints; ++i) {
	  			outFile << fixed << setprecision(6) << ptime[i];
	  			if (i < maxpoints - 1) outFile << ", ";
	  		}
      		outFile << "};\n";
      		// Write the dsampl array
      		outFile << "double drv[" << maxpoints << "] = {";
      		for (int i = 0; i < maxpoints; ++i) {
      			outFile << fixed << setprecision(6) << dsampl[i];
      			if (i < maxpoints - 1) outFile << ", ";
      		}
      		outFile << "};\n";
      		outFile.close();
      		cout << "Waveform data written to " << filename << endl;
      	}
#endif

#ifdef DEBUGMSG
      cout<<RED<<"Threshold "<<Thresholds[ci]<<endlr;
	  //print smoothing and integration points
      	// cout<<" Smoothing points = "<<nsmooth<<endl;
      	// cout<<" Integration points = "<<nint<<endl;
#endif

      	// Give trigger window to the function
      	if (trigger_windows_iterator == trigger_windows.end()) {
      		break;
      	}

      	pair<double, double> trigger_window = *trigger_windows_iterator;
      	// cout << "Trigger window: (" << trigger_window.first << ", " << trigger_window.second << ")" << endl;
      	int i_start = convert_x_to_index(ptime, maxpoints, trigger_window.first);
      	// cout << "i_start: " << i_start << endl;
      	int i_end = convert_x_to_index(ptime, maxpoints, trigger_window.second);
      	// cout << "i_end: " << i_end << endl;

      	if (i_end < maxpoints - 50) {
      		// cout << "Analyzing trigger window..." << endl;
      		AnalysePicosecBounds(maxpoints, evNo, sampl, dt, i_start, i_end, ppar);    /// all the analysis is done here!!!!
      		successfulFits_sigmoid += ppar->SigmoidfitSuccess;
      		totalFits_sigmoid++;
      		successfulFits_double_sigmoid += ppar->doubleSigmoidfitSuccess;
      		totalFits_double_sigmoid++;
      	}

      	++trigger_windows_iterator;
      }
      else
      {
          ti = AnalyseLongPulse(maxpoints,sampl,dsampl,ppar,Thresholds[ci], dt, ti);
          if (ti<0) break;
      }

      new ((*sparArr[ci])[ntrigs]) PEAKPARAM(*ppar);// only for storing the peak parameters to the output tree

   	  AddPar(ppar,&spar[ci][ntrigs]); // this is for the output visualization plots and histograms
///*************************************************************

      double risetime = spar[ci][ntrigs].t90-spar[ci][ntrigs].t10;
      if (isMM) risetime = spar[ci][ntrigs].risetime;
      
      double pulseDuration = (spar[ci][ntrigs].ftime_pos - spar[ci][ntrigs].stime_pos)*dt;

      
	  hAMPL[ci]->Fill(spar[ci][ntrigs].ampl*mV);

// 	  hsCH[ci]->Fill(spar[ci][ntrigs].charge-spar[ci][ntrigs].bslch);
	  hRT[ci]->Fill(risetime);
      hPD[ci]->Fill(pulseDuration);

      if (!pulsedBeam)
	  {
	    hRateStructure[ci]->Fill(spar[ci][ntrigs].t10*microsec);//t10 is ns, so to pass to us
	  }

	  //cout << "rate structure = " << spar[ci][ntrigs].t10/1000. << endl;
	  h2dAmplvsRT[ci]->Fill(spar[ci][ntrigs].ampl*mV,risetime);
	  h2dAmplvsPD[ci]->Fill(spar[ci][ntrigs].ampl*mV,pulseDuration);

      if (isMM)  /// Micromegas specific plots
      {
         hDAMPL[ci]->Fill(spar[ci][ntrigs].dampl*mV);
         htotCH[ci]->Fill(spar[ci][ntrigs].totcharge);
         htotCHfixed[ci]->Fill(spar[ci][ntrigs].totchargefixed);
         heCHfixed[ci]->Fill(spar[ci][ntrigs].echargefixed);
         heCHfit[ci]->Fill(spar[ci][ntrigs].echargefit);
         heCH[ci]->Fill(spar[ci][ntrigs].echarge);
         hionCH[ci]->Fill(spar[ci][ntrigs].ioncharge);

         htotCHovAmpl[ci]->Fill((spar[ci][ntrigs].totcharge)/(spar[ci][ntrigs].ampl*mV) );
         hAmplovDampl[ci]->Fill((spar[ci][ntrigs].ampl*mV)/(spar[ci][ntrigs].dampl*mV) );
         heCHovAmpl[ci]->Fill((spar[ci][ntrigs].echarge)/(spar[ci][ntrigs].ampl*mV) );
         heCHovtotCh[ci]->Fill((spar[ci][ntrigs].echarge)/(spar[ci][ntrigs].totcharge) );
         // cout<<YELLOW<<"spar[ci][ntrigs].ioncharge = "<<spar[ci][ntrigs].ioncharge<<endl;
         h2deCHvsiCH[ci]->Fill((spar[ci][ntrigs].ioncharge),(spar[ci][ntrigs].echarge) );
	     // cout<<YELLOW<<"spar[ci][ntrigs].ioncharge = "<<spar[ci][ntrigs].ioncharge<<endl;
      	//cin.get();
		 h2dtotCHvsAmpl[ci]->Fill(spar[ci][ntrigs].ampl*mV,spar[ci][ntrigs].totcharge);
         h2dAmplvsDampl[ci]->Fill(spar[ci][ntrigs].dampl*mV,spar[ci][ntrigs].ampl*mV);
         h2deCHvsAmpl[ci]->Fill(spar[ci][ntrigs].ampl*mV,spar[ci][ntrigs].echarge);
      	 h2deCHvstotCh[ci]->Fill(spar[ci][ntrigs].echarge,spar[ci][ntrigs].totcharge);
#ifdef DEBUGMSG         
         cout<<BLUE<<"spar[ci][ntrigs].totcharge = "<<spar[ci][ntrigs].totcharge<<endl;
         cout<<BLUE<<"spar[ci][ntrigs].totchargefixed = "<<spar[ci][ntrigs].totchargefixed<<endl;
         cout<<BLUE<<"spar[ci][ntrigs].echarge = "<<spar[ci][ntrigs].echarge<<endl;
         cout<<BLUE<<"spar[ci][ntrigs].echargefit = "<<spar[ci][ntrigs].echargefit<<endl;
         cout<<BLUE<<"spar[ci][ntrigs].echargefixed = "<<spar[ci][ntrigs].echargefixed<<endl;
         cout<<BLUE<<"risetime = "<<risetime<<endl;
#endif
        
      }
      else  /// other detector plots
      {
        htotCH[ci]->Fill(spar[ci][ntrigs].charge);
        h2dtotCHvsAmpl[ci]->Fill(spar[ci][ntrigs].ampl*mV,spar[ci][ntrigs].charge);
        htotCHovAmpl[ci]->Fill((spar[ci][ntrigs].charge)/(spar[ci][ntrigs].ampl*mV));
      }
      

      
      int cut1 = spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl>250.;
	  cut1 =1;
	  if (oscsetup->AmplifierNo[ci]==1)
      {
          cut1 = spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl>1.5 && spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl<6.;
          cut1 *= risetime>0.25 && risetime<2.0;
      }
	  if (oscsetup->AmplifierNo[ci]==3)
      {
          cut1 = spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl>400 && spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl<850;
          cut1 *= (TOT>80);
          cut1 *= (risetime>45);
          cut1 *= (Width>170);
//           cut1=!cut1;
      }
	  if (oscsetup->AmplifierNo[ci]==2)
      {
          cut1 = spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl>400 && spar[ci][ntrigs].charge/spar[ci][ntrigs].ampl<850;
          cut1 *= (TOT>80);
          cut1 *= (risetime>45);
          cut1 *= (Width>170);
//           cut1=!cut1;
      }
      cut1=1;

      hPW[ci]->Fill(spar[ci][ntrigs].width);
	  if (ntrigs>0)  /// here we take the dt between concecutive peaks on the same waveform!
	  {
	    hDt[ci]->Fill((spar[ci][ntrigs].t10-spar[ci][ntrigs-1].t10)*1e-9);//in seconds
	  }
// 	  else
          tnow = 1.*epoch + nn*1e-9 + spar[ci][0].t10*1e-9;   /// normally we arrive here only in the first peak per pulse!!!!

	  // //doubling histos for nthreshold
	  if(1 || ((spar[ci][ntrigs].ampl) > (peTh) && (spar[ci][ntrigs].charge) > (peThCh) && cut1))//(spar[ci][ntrigs].tot>60) && (spar[ci][ntrigs].tot<100)(spar[ci][ntrigs].width>60) && (spar[ci][ntrigs].width<100))
	  {
	      if (isMM)
          {
            hCutsCH[ci]->Fill(spar[ci][ntrigs].totcharge);
            htotCHovAmplCuts[ci]->Fill((spar[ci][ntrigs].totcharge)/(spar[ci][ntrigs].ampl*mV));
            h2dCutsCHevolution[ci]->Fill(tnow,spar[ci][ntrigs].totcharge);
          }
          else
          {
            hCutsCH[ci]->Fill(spar[ci][ntrigs].charge);
            htotCHovAmplCuts[ci]->Fill((spar[ci][ntrigs].charge)/spar[ci][ntrigs].ampl*mV);
            h2dCutsCHevolution[ci]->Fill(tnow,spar[ci][ntrigs].charge);
          }
          hCutsAMPL[ci]->Fill(spar[ci][ntrigs].ampl*mV);
	      h2dCutsAmplEvolution[ci]->Fill(tnow,spar[ci][ntrigs].ampl*mV);

          if (ntrigsCuts>0)  /// here we take the dt between concecutive peaks on the same waveform!
          {
            hCutsDt[ci]->Fill((spar[ci][ntrigsCuts].t10-spar[ci][ntrigsCuts-1].t10)*1e-9);
          }

           ntrigsCuts++;
	    }
	  //------------------------

	  ntrigs++;
	}///end loop ti   ==> scan along one waveform to find multiple peaks

      npeaks[ci] = ntrigs;
      ngoodPeaks[ci] = ntrigsCuts;
#ifdef DEBUGMSG
      cout<<'C' << ci+1 <<" :   Found "<<ntrigsCuts<<" goodpeaks"<<endl;
#endif
      //       if (ngoodPeaks[ci]>1) cout <<"event "<<eventNo<<" n = "<<ngoodPeaks[ci]<<"  "<< epoch - epochS<< "  "<< evtime - epochS<< endl;

     
      ntrigsTot[ci] += npeaks[ci];
      ngoodTrigsTot[ci] += ngoodPeaks[ci];
      //cout << "npeaks[ci]" << npeaks[ci] << endl;

/// fill in the XY plots. Single good peak is required for the cuts
      for (int i=0;i<eventTracks;i++)
      {
        h2dXY[ci]->Fill(hitX_C[ci][i],hitY_C[ci][i]);
      }
      
      if (ngoodPeaks[ci]==1)
      {
          h2dXYcuts[ci]->Fill(hitX_C[ci][0],hitY_C[ci][0]);
      }

/// fill in rate evolution plots      
      double ts = ptime[0];//ns
      double tf = ptime[maxpoints-1];//ns
      double tlive = (tf-ts) * 1e-9;//to make it sec
      double peakspersec = 1.* npeaks[ci] / tlive;
      double goodpeakspersec = 1.*ntrigsCuts / tlive;
///       double tnow = epoch + nn*1e-9 + spar[npeaks[ci]-1].t10*1e-9;  /// time of last peak in event or of the unique event if short frame//sec

      long double depoch = 1.* (long double) epoch;
      long double tnowgoodpeaks = depoch + nn*1e-9 + spar[ci][ntrigsCuts-1].t10*1e-9 ;  /// used for dt calculation for goodpeaks


      hRateEvolution[ci]->Fill(evtime,npeaks[ci]);
      hRateEvolutionCheck[ci]->Fill(evtime);

      if (callibrationRun)
      {
         //hBkgGoodPeaksRateEvolution[ci]->Fill(evtime,ngoodPeaks[ci]);
         if (ngoodPeaks[ci]>0)
         {
           hRate[ci]->Fill(ngoodPeaks[ci]);
           hRateEvolutionCuts[ci]->Fill(evtime,ngoodPeaks[ci]);
         }
      }
      else
      {

        hSRStriggerStruct[ci]->Fill(ttrig*microsec);//t10 is ns, so to pass to us
        
        hRateEvolutionCuts[ci]->Fill(evtime,ngoodPeaks[ci]);
      }

      if (eventNo>0 && npeaks[ci]>0)  /// here we take the dt between the first peaks of concecutive waveforms
      {
          //intantaneous flux ("peak flux")
          //cout<<"tnow = "<<tnow<<"  nn = "<<nn*1e-9<<endl;
          dtlast[ci] = tnow - tlast[ci];
//     	    cout<<GREEN<<'C'<<ci+1<<" :   tlast = "<<tlast[ci]<<" + "<< + nn*1e-9+spar[ci][npeaks[ci]-1].t10*1e-9<<"  dtlast = "<<dtlast[ci]<<endlr;
          if (tlast[ci]>0)
              hDt[ci]->Fill(dtlast[ci]);
         tlast[ci] = tnow ;
      }


      if (eventNo>0 && ngoodPeaks[ci]>0)
      {
          dtlastgoodpeaks[ci] = tnowgoodpeaks - tlastgoodpeaks[ci];
          hCutsDt[ci]->Fill(dtlastgoodpeaks[ci]);
          tlastgoodpeaks[ci] = tnowgoodpeaks ;
      }

 /// 3rd IF DRAW
      //----------------------
    if (draw)
	{
    //cout<<endl<<"Entering 3rd if (draw) after the analysis________________________"<<endl<<endl;;
	  cout<<"Found "<<ntrigs<<" pulses "<<endl;

	  evdcanv[ci]->cd(1);
	TGraph *graphSum = new TGraph(maxpoints,ptime,amplC[ci]);
	graphSum->SetMarkerColor(clr[ci+2]);
    graphSum->SetLineColor(clr[ci+2]);
    graphSum->SetFillColor(0);
    graphSum->SetLineWidth(2);
    graphSum->SetLineStyle(7);
    sprintf(cname,"Raw signal channel %d", evNo, ci+1);
    graphSum->SetTitle(cname);
    graphSum->GetXaxis()->SetTitle("Time [ns]");
    graphSum->GetYaxis()->SetTitle("Amplitude [V]");
	// graphSum->Draw("apl");

	TGraph *sgraphSum = new TGraph(maxpoints,ptime,sampl);
	sgraphSum->SetLineColor(7);
	sgraphSum->SetFillColor(0);
	sgraphSum->SetLineWidth(2);
	sprintf(cname,"Smoothed signal 1 event %d",evNo);
	sgraphSum->SetTitle(cname);
    sgraphSum->GetXaxis()->SetTitle("Time [ns]");
    sgraphSum->GetYaxis()->SetTitle("Amplitude [V]");
	// sgraphSum->Draw("pl");

	  evdcanv[ci]->cd(2);
    //TGraph* graph22 = new TGraph(maxpoints,ptime,dsampl);
	  TGraph* graph22 = new TGraph(maxpoints,ptime, dampl[ci]);
    double maxelement_dampl[4]= {0,0,0,0};
    maxelement_dampl[ci] = TMath::MaxElement(maxpoints, dampl[ci]);
	  graph22->SetMarkerColor(1);
	  graph22->SetLineColor(1);
	  graph22->SetFillColor(0);
	  graph22->SetLineWidth(2);
	  graph22->SetLineStyle(7);
	  sprintf(cname,"Derivative\n");
	  graph22->SetTitle(cname);
	  //graph22->Draw("pl");
	  graphSum->SetMarkerColor(clr[ci+2]);
	  graphSum->SetLineColor(clr[ci+2]);
	  graphSum->SetFillColor(0);
	  graphSum->SetLineWidth(2);
	  graphSum->SetLineStyle(7);
	  sprintf(cname,"Raw signal channel %d", evNo, ci+1);
	  graphSum->SetTitle(cname);
	  graphSum->GetXaxis()->SetTitle("Time [ns]");
	  graphSum->GetYaxis()->SetTitle("Amplitude [V]");
	  graphSum->Draw("pl");


    evdcanv[ci]->cd(3);

    TMultiGraph *mg = new TMultiGraph();
    int nint = 20; //default
    nint = N_INTEGRATION_POINTS;
    //cout<<"N integration points = "<<N_INTEGRATION_POINTS<<endl;
    //double DTI2 = 2.;  ///default
    double DTI2 = 150.; //time window for integration in ns
    nint = TMath::FloorNint(DTI2/dt)+1; // make it int at >= of the setted time window
    //cout<<"Integration points second time with DTI2 = "<<nint<<endl;
    double intgr;
    intgr = IntegratePulse(maxpoints,idamplC,iampl,dt,nint*dt);
    // continue;
    maxd[ci]=TMath::MaxElement(maxpoints,iampl);
    mind[ci]=TMath::MinElement(maxpoints,iampl);
    integralh = new TGraph(maxpoints,ptime,iampl);
    sprintf(cname,"Event %d - Integral %g of C%d\n",evNo, nint*dt,ci+1);
    integralh->SetMarkerColor(clr[ci+4]);
    integralh->SetLineColor(clr[ci+4]);
    integralh->SetFillColor(0);
    integralh->SetLineWidth(2);
    integralh->SetLineStyle(10);
    integralh->SetTitle(cname);
    //integralh->SetMaximum(maxelement_dampl[ci]+0.1);
    //integralh->SetMinimum(mind[ci]-0.1);
    mg->Add(integralh);
    //integralh->Draw("apl");
    derivative2->SetMarkerColor(clr[ci+11]);
    derivative2->SetLineColor(clr[ci+11]);
    derivative2->SetFillColor(0);
    derivative2->SetLineWidth(2);
    derivative2->SetLineStyle(10);
    sprintf(cname,"Derivative C%d\n",ci+1);
    derivative2->SetTitle(cname);
    //derivative2->Draw("pl");
	  mg->Add(derivative2);
    graphSum->SetMarkerColor(clr[ci+2]);
    graphSum->SetLineColor(clr[ci+2]);
    graphSum->SetFillColor(0);
    graphSum->SetLineWidth(2);
    graphSum->SetLineStyle(10);
    sprintf(cname,"Raw signal channel %d", evNo, ci+1);
    graphSum->SetTitle(cname);
    mg->Add(graphSum);

    DTI2 = 150.; //time window for integration in ns
    nint = TMath::FloorNint(DTI2/dt)+1; // make it int at >= of the setted time window
    intgr = IntegratePulse(maxpoints,dampl[ci],iampl,dt,nint*dt);
    integralh = new TGraph(maxpoints,ptime,iampl);
    sprintf(cname,"Event %d - Integral Derivative %g of C%d\n",evNo, nint*dt,ci+1);
    integralh->SetMarkerColor(clr[ci+1]);
    integralh->SetMarkerStyle(7);
    integralh->SetLineColor(clr[ci+1]);
    integralh->SetFillColor(0);
    integralh->SetLineWidth(2);
    integralh->SetLineStyle(10);
    integralh->SetTitle(cname);
    //integralh->SetMaximum(maxelement_dampl[ci]+0.1);
    //integralh->SetMinimum(mind[ci]-0.1);
    mg->Add(integralh);
    mg->GetXaxis()->SetTitle("Time [ns]"); //this is to the integral graph only
    mg->GetYaxis()->SetTitle("Amplitude [V]");

    mg->Draw("ap");


	  if (waveform!=0)
	  {
      TH1F* h1 = (TH1F*) waveform->GetHistogram();
      TH1F* h2 = (TH1F*) derivative->GetHistogram();
      TH1F* h3 = (TH1F*) integralh->GetHistogram();


      h1->GetXaxis()->SetTitle("Time [ns]");
      h1->GetYaxis()->SetTitle("Amplitude [V]");
      h2->GetXaxis()->SetTitle("Time [ns]");
      h2->GetYaxis()->SetTitle("Amplitude [V]");


	    evdcanv[ci]->Modified();
	    evdcanv[ci]->cd(2);
	    evdcanv[ci]->Update();
	  }

// 	  for (int i=0;i<ntrigs;i++)
	  for (int i=0;i<npeaks[ci];i++)
	  {
	     evdcanv[ci]->cd(1);
	      TLine *line1 = new TLine(spar[ci][i].t10,-spar[ci][i].ampl,spar[ci][i].tb10,-spar[ci][i].ampl);
//         TLatex *label1 = new TLatex(0.5, 0.5, "Peak maximum position");
//         label1->SetLineColor(4);
        line1->SetLineColor(4);
//         label1->Draw();
        line1->Draw();
	      TLine *line2 = new TLine(spar[ci][i].stime_pos*dt,spar[ci][i].bsl,spar[ci][i].ftime_pos*dt,spar[ci][i].bsl);
// 	      TLatex *label2 = new TLatex(0.5, 0.5, "Starting/Ending points of waveform");
        line2->SetLineColor(2);
//         label2->SetLineColor(2);
//         label2->Draw();
        line2->Draw();
	      TLine *line3 = new TLine(spar[ci][i].ttrig*dt,Thresholds[ci],spar[ci][i].ttrig*dt+spar[ci][i].tot[0],Thresholds[ci]);
// 	      TLatex *label3 = new TLatex(0.5, 0.5, "Starting/Ending points of E-peak");
          //cout<<"(x1,y1) = "<<spar[ci][i].ttrig<<" , "<<Thresholds[ci] << "  (x2,y2) = "<<spar[ci][i].ttrig+spar[ci][i].tot[0]<<" , "<<Thresholds[ci]<<endl;
        line3->SetLineColor(3);
//         label3->SetLineColor(3);
//         label3->Draw();
        line3->Draw();
        gPad->BuildLegend(0.75,0.8,0.99,0.99);

        // 	    //cout<<i<<" \t "<<spar[i].ampl<<endl;
        evdcanv[ci]->cd(2);
        TLine *line4 = new TLine(spar[ci][i].t10,-spar[ci][i].ampl,spar[ci][i].tb10,-spar[ci][i].ampl);
        line4->SetLineColor(4);
//         TLatex *label4 = new TLatex(0.5, 0.5, "Peak maximum position");
//         label4->SetLineColor(4);
//         label4->Draw();
        line4->Draw();

        TLine *line5 = new TLine(spar[ci][i].stime_pos*dt,spar[ci][i].bsl,spar[ci][i].ftime_pos*dt,spar[ci][i].bsl);
        line5->SetLineColor(2);
//         TLatex *label5 = new TLatex(0.5, 0.5, "Starting/Ending points of waveform");
//         label5->SetLineColor(2);
//         label5->Draw();
        line5->Draw();

        TLine *line6 = new TLine(spar[ci][i].ttrig,Thresholds[ci],spar[ci][i].ttrig+spar[ci][i].tot[0],Thresholds[ci]);
        line6->SetLineColor(3);
//         TLatex *label6 = new TLatex(0.5, 0.5, "Starting/Ending points of E-peak");
//         label6->SetLineColor(3);
//         label6->Draw();
        line6->Draw();

        TLine *line7 = new TLine(spar[ci][i].maxtime,-spar[ci][i].bsl, spar[ci][i].maxtime, -spar[ci][i].ampl);
        line7->SetLineColor(7);
//         TLatex *label7 = new TLatex(0.5, 0.5, "Vertical region of peak max");
//         label7->SetLineColor(8);
//         label7->Draw();
        line7->Draw();

        gPad->BuildLegend(0.75,0.8,0.99,0.99);

        evdcanv[ci]->Modified();
        evdcanv[ci]->Update();

        /// Sigmoid Fit Draw
	  	// cout<<RED<<"Entering Sigmoid if"<<endlr;
	  	// cout<<YELLOW<<ci<<endlr;
	  	// cout<<"npeaks = "<<npeaks[ci]<<endl;

       if (SIGMOIDFIT[ci])
       {
         // cout<<RED<<"Sigmoidfit flag executted for channel "<<ci+1<<endlr;
         fitcanv[ci]->cd(1);
         // cout<<RED<<"spar[channel][trig] "<<spar[ci][i].tot_sig_end_pos<<endlr;
         TimeSigmoidDraw(maxpoints, amplC[ci], ptime, &spar[ci][i], evNo, fitcanv[ci]);
       	// cout<<RED<<"Exiting TimeSigmoid DRaw()"<<endlr;
         gStyle->SetOptFit(1111);
         fitcanv[ci]->Modified();
         fitcanv[ci]->BuildLegend(0.75,0.8,0.99,0.99);
         fitcanv[ci]->Update();

         fitcanv[ci]->cd(2);
       	// cout<<RED<<"Entering FullSigmoidDraw()"<<endlr;
         FullSigmoidDraw(maxpoints, amplC[ci], ptime, dt, &spar[ci][i], evNo, fitcanv[ci]);
         gStyle->SetOptFit(1111);
         fitcanv[ci]->Modified();
         fitcanv[ci]->BuildLegend(0.75,0.8,0.99,0.99);
         fitcanv[ci]->Update();

       }
	  	// cout<<RED<<"Entering Doublesigmoid if"<<endlr;

	  	if (DOUBLESIGMOIDFIT[ci]==false)
	  	{
	  		//cout<<RED<<"Sigmoidfit flag executted for channel "<<ci+1<<endlr;
	  		fitcanv[ci]->cd(1);
	  		TimeSigmoidDraw(maxpoints, amplC[ci], ptime, &spar[ci][i], evNo, fitcanv[ci]);
	  		gStyle->SetOptFit(1111);
	  		fitcanv[ci]->Modified();
	  		fitcanv[ci]->BuildLegend(0.75,0.8,0.99,0.99);
	  		fitcanv[ci]->Update();
	  	}

    }
    	// cout<<RED<<"npeaks = "<<npeaks[ci]<<endlr;

	  evdcanv[ci]->Modified();
	  evdcanv[ci]->Update();

	  evdcanv[ci]->cd(2);
      gPad->BuildLegend(0.75,0.8,0.99,0.99);
      gPad->Update();
	  //    continue;

	  if (integralh!=NULL)
	  {
        evdcanv[ci]->cd(3);
        TH1F* h2 = (TH1F*) integralh->GetHistogram();
	    h2->SetMaximum(TMath::MaxElement(4,maxd)+0.1*TMath::Abs(TMath::MaxElement(4,maxd)) );
	    h2->SetMinimum(TMath::MinElement(4,mind)-0.1*TMath::Abs(TMath::MinElement(4,mind)) );
	    gPad->Modified();
	    gPad->BuildLegend(0.75,0.8,0.99,0.99);
	    gPad->Update();
	  }



	}
      ///---------------end 2nd if(draw)--------------------------------

    if (ntrigs>20000)
	{
	  cout<<ntrigs<<" pulses in event "<< eventNo <<" channel "<<ci+1 << endl;
	  return (ntrigs);
	}


    
  }
       /////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 /// END of the channel loop (ci) inside events tree entry processing
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////
	if (eventNo % (evpm) == 0)
    {
        cout<<BLUE<<"Event "<<setw(6)<< eventNo<<endl;
        for (int ci = 0; ci < 4; ++ci) 
		{
			cout<<"Channel C"<<ci+1<<" , found ntrigs ="<<npeaks[ci]<<" pulses , "<<ntrigsTot[ci]<<" total ,  after cuts ="<<ngoodPeaks[ci]<<" pulses , "<<ntrigsTot[ci]<<" total "<<endlr;
		}
	}
    //  for (int ci=0;ci<4; ci++) cout<<RED<<"Size of Array = "<<sparArr[ci]->GetEntriesFast()<<endlr;

    otree->Fill();
  
    eventNo++;
    //if (eventNo > 1250) break;

      
  } /// end of the tree while (eventNo < nevents)

	for (int ci = 0; ci < 4; ++ci) {
		cout<<BLUE<<"Total number of pulses in channel "<<ci+1<<" = "<<ntrigsTot[ci]<<endl;
	}
  	cout<<endl<<"Succesful Fits = "<<successfulFits_sigmoid<<endl;
  	cout<<"Total Fits = "<<totalFits_sigmoid<<endl;
  	cout<<"Percentage of successful fits = "<<100.*successfulFits_sigmoid/totalFits_sigmoid<<"%"<<endl;

  	cout<<endl<<"Succesful Fits Double Sigmoid= "<<successfulFits_double_sigmoid<<endl;
  	cout<<"Total Fits Double Sigmoid = "<<totalFits_double_sigmoid<<endl;
  	cout<<"Percentage of successful fits Double Sigmoid= "<<100.*successfulFits_double_sigmoid/totalFits_double_sigmoid<<"%"<<endl;

	cout<<endl<<"Total secondary pulses rejected = "<<total_secondary_count<<endl;
	cout<<endl<<"Total thin pulses rejected = "<<total_thin_count<<endl;
	cout<<endl<<"Percentage of secondary pulses rejected in all the pulses = "<<100.*total_secondary_count/(ntrigsTot[0]+ntrigsTot[1]+ntrigsTot[2]+ntrigsTot[3])<<"%"<<endl;
	cout<<endl<<"Percentage of thin pulses rejected in all the pulses = "<<100.*total_thin_count/(ntrigsTot[0]+ntrigsTot[1]+ntrigsTot[2]+ntrigsTot[3])<<"%"<<endl;

  if (draw) 
  {
    cout<<endl<<"Entering the 4rth if(draw)______________________________________________"<<endl<<endlr;;
    gSystem->ChangeDirectory(tmpdir);
    return 0;  /// finished the event display, so exit without drawing the histograms (empty)
  }
/// fraph drawing follows.
// return 0;  
  
  bool firstactive=true;
  

  //For CASE 1 (pulsed beam)
  //hRateEvolution->Scale(1./timebinwidth);
  int abins = 100;
  double *err0;
  for (int ci = 0; ci < 4; ++ci)
  {
    if (!active[ci]) continue;
      
    bool isMM = kFALSE;
    if (strncmp(oscsetup->DetName[ci], "MM", 2) == 0) isMM = kTRUE;  

    cout<<BLUE<<"Preparing rate evolution for channel "<<MAGENTA<<ci+1<<endlr;
    if (firstactive)
    {
      abins = hRateEvolution[ci]->GetNbinsX();
      err0 = new double[abins+1];
      for(int i=0;i<=abins;i++)
          err0[i]=0;
      
      firstactive=false;
    }
    
    if (callibrationRun)
    {
        hRateEvolutionCheck[ci]->SetError(err0);
        cout<<GREEN<<"done"<<endlr;

        cout<<BLUE<<"Normalizing callibration rate evolution for channel "<<MAGENTA<<ci+1<<endlr;
        DivideH(hRateEvolution[ci],hRateEvolutionCheck[ci]);
        DivideH(hRateEvolutionCuts[ci],hRateEvolutionCheck[ci]);
//         cout<<BLUE<<"Scaling for timebinwidth = "<<MAGENTA<<timebinwidth<<endlr;
//         ScaleHistoErr(hSparkEvolution[ci],1.0/timebinwidth);
        cout<<GREEN<<"done"<<endlr;
    }
    

    /// divide with binwidth to make rates per second
    for (int i = 0; i<hRateEvolutionCuts[ci]->GetNbinsX();i++)
    {
        
        double y = hRateEvolutionCuts[ci]->GetBinContent(i);
      
        y /= timebinwidth; /// convert it to events per time periode width
        if (y>0)
           hRateCuts[ci]->Fill(y);
    }

    
    TString channel = TString::Itoa(ci+1,10);
    TString chname = ctypes+"_C"+channel;

    cout<<BLUE<<"Entering plots directory "<<MAGENTA<< plotdirname <<BLUE<<" for channel "<<MAGENTA<<ci+1<<endlr;
    gSystem->ChangeDirectory(plotdirname);
    cout<<GREEN<<"done"<<endlr;
   
/// ///////////////////////////////////
///  Amplitude and charge     ///
/// ///////////////////////////////////

    TCanvas *campl = new TCanvas("Amplitude"+chname,"Amplitude "+chname);
    campl->SetLogy();
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TF1* f2 = new TF1("fexp","[0]*exp(-[1]*x)",0.8*Thresholds[ci],peTh*2.);
    f2->SetTitle("noise fit");
    f2->SetLineColor(kGreen+2);
    f2->SetParameter(1,100.);
    f2->SetParLimits(1,0.01,1e7);

    hCutsAMPL[ci]->GetXaxis()->SetRangeUser(peTh,amplMax[ci]*mV*0.75);
    double maxampl = hCutsAMPL[ci]->GetMaximum();
    int maxabin = hCutsAMPL[ci]->GetMaximumBin();
    hCutsAMPL[ci]->GetXaxis()->UnZoom();
    double maxax = hCutsAMPL[ci]->GetBinLowEdge(maxabin);
    double minax = peTh*mV;
    for (int i=maxabin;i>1;i--)
    {
        double y = hCutsAMPL[ci]->GetBinContent(i);
        if (y<0.5*maxampl)
        {
            minax = hCutsAMPL[ci]->GetBinLowEdge(i);
            break;
        }
    }
    
    minax = peTh*mV+0.05;

    TF1 *fpolyaAmpl = new TF1 ("PolyaAmpl"+channel,"([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])",minax,amplMax[ci]);
    fpolyaAmpl->SetTitle("Polya fit (ampl.)");
    fpolyaAmpl->SetParName(0,"c");
    fpolyaAmpl->SetParName(1,"#bar{Q}");
    fpolyaAmpl->SetParName(2,"#theta");
    fpolyaAmpl->SetLineColor(kGreen+2);
    fpolyaAmpl->SetRange(minax,amplMax[ci]);
    fpolyaAmpl->SetNpx(1000);

    fpolyaAmpl->SetParameter(0,1000); // amplitude
    fpolyaAmpl->SetParameter(1,30); //Charge
    fpolyaAmpl->SetParameter(2,30); //theta
    fpolyaAmpl->SetParLimits(2,0.001,5000); //theta

    hCutsAMPL[ci]->Sumw2(); 
    hCutsAMPL[ci]->Fit(fpolyaAmpl,"QR0","",minax,amplMax[ci]);
    hCutsAMPL[ci]->Fit(fpolyaAmpl,"QR0","",minax,amplMax[ci]);
    hCutsAMPL[ci]->Fit(fpolyaAmpl,"QB0","",minax,amplMax[ci]);
//     hCutsAMPL->Fit(fpolyaA,"BM+");

    cout<<"\n\n==> minax = "<<minax<<endl;
// hAMPL->SetStats(0);
    gStyle->SetOptStat(100);
    gStyle->SetOptFit(111);
    cout<<"Fit "<<channel<<"  "<<1<<endl;
    hCutsAMPL[ci]->Fit(fpolyaAmpl,"B0","",minax,amplMax[ci]);
    cout<<"Fit "<<channel<<"  "<<2<<endl;
    hCutsAMPL[ci]->Fit(fpolyaAmpl,"RM","",minax,amplMax[ci]);

    hCutsAMPL[ci]->SetLineColor(2);
    hCutsAMPL[ci]->SetMaximum(hAMPL[ci]->GetMaximum()*2.);
    hCutsAMPL[ci]->Draw("");
    hAMPL[ci]->Draw("same");
    
    if (isMM)
    {
       hDAMPL[ci]->SetLineColor(kBlack);
       hDAMPL[ci]->Draw("same");
    }

    
    fpolyaAmpl->SetRange(minax,amplMax[ci]);
    fpolyaAmpl->Draw("same");
    // f2->Draw("same");
    campl->Update();
    campl->BuildLegend(0.6,0.5,0.9,0.65);
  
//  TPad* subpad1=(TPad*)cch->GetPrimitive(Form("%s_%d",rcanv->GetName(),4));
    campl->Update();  /// VERY IMPORTANT!!!!! crashes without update!!!
    TPaveStats* st11 = (TPaveStats*)campl->GetPrimitive("stats");
    st11->SetX1NDC(0.6); //new x start position
    st11->SetY1NDC(0.65); //new x end position 
    st11->SetX2NDC(0.9); //new x start position
    st11->SetY2NDC(0.9); //new x end position 
    st11->SetTextSize(0.031);
    campl->Modified();

    campl->SaveAs(".png");
    campl->SaveAs(".pdf");

 
    TCanvas *cch = new TCanvas("Charge"+chname,"Charge "+chname);

    double minchx = peThCh;

    minchx = peThCh + +1.5*ctoaf/1000.;
    minchx = minax*ctoaf/1000.*10;
    minchx = 150.;
    cout<<"\n\n==> minchx = "<<minchx<<endl;
    // fpolyaAmpl->SetRange(minchx,chMax[ci]);

    TF1 *fpolyaCH = new TF1 ("PolyaCharge"+channel,"([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])",minchx,chMax[ci]);
    fpolyaCH->SetTitle("Polya fit (charge)");
    fpolyaCH->SetParName(0,"c");
    fpolyaCH->SetParName(1,"#bar{Q}");
    fpolyaCH->SetParName(2,"#theta");
    fpolyaCH->SetLineColor(kGreen+2);
    fpolyaCH->SetRange(minax,amplMax[ci]);
    fpolyaCH->SetNpx(1000);

    fpolyaCH->SetParameter(0,1000); // amplitude
    fpolyaCH->SetParameter(1,500); //Charge
    fpolyaCH->SetParameter(2,30); //theta
    fpolyaCH->SetParLimits(2,0.001,5000); //theta


    //hCutsCH[ci]->SetStats(1);
    cout<<"Fit "<<channel<<"  "<<6<<endl;
    //hCutsCH[ci]->Fit(fpolyaCH,"B0","",minchx,chMax[ci]);
    cout<<"Fit "<<channel<<"  "<<7<<endl;
    //hCutsCH[ci]->Fit(fpolyaCH,"BM","",minchx,chMax[ci]);

    hCutsCH[ci]->SetMaximum(2.*htotCH[ci]->GetMaximum());
    cch->SetLogy(); 
    hCutsCH[ci]->SetLineColor(2);
    hCutsCH[ci]->Draw();
//   htotCH->Sumw2();
//   hsCH->Draw("same");
//  hsCH[ci]->SetLineColor(kMagenta);
    htotCH[ci]->Draw("same");
    fpolyaCH->SetRange(minchx,chMax[ci]);
    fpolyaCH->Draw("same");

    if (isMM)
    {
       heCH[ci]->SetLineColor(kMagenta);
       heCHfit[ci]->SetLineColor(kBlack);
       heCHfixed[ci]->SetLineColor(kBlue-2);
       hionCH[ci]->SetLineColor(kGreen+1);
       heCH[ci]->Draw("same");
       heCHfit[ci]->Draw("same");
       heCHfixed[ci]->Draw("same");
       hionCH[ci]->Draw("same");
       hCutsCH[ci]->SetMaximum(hCutsCH[ci]->GetMaximum()*5);
    }
    
    cch->Update();
    cch->BuildLegend(0.6,0.5,0.9,0.65);
    
   
  
//  TPad* subpad1=(TPad*)cch->GetPrimitive(Form("%s_%d",rcanv->GetName(),4));
    cch->Update();  /// VERY IMPORTANT!!!!! crashes without update!!!
    TPaveStats* st1 = (TPaveStats*)cch->GetPrimitive("stats");
    // if (st==NULL) return -3243;
    st1->SetX1NDC(0.6); //new x start position
    st1->SetY1NDC(0.65); //new x end position 
    st1->SetX2NDC(0.9); //new x start position
    st1->SetY2NDC(0.9); //new x end position 
    st1->SetTextSize(0.031);
    cch->Modified();
    
    cch->SaveAs(".png");
    cch->SaveAs(".pdf");
///________________________________________________________________
  
  
/// ///////////////////////////////////
///  		Rates
/// ///////////////////////////////////

    gSystem->ChangeDirectory(plotdirname);
    
    TCanvas *rcanv = new TCanvas("Rates_"+chname,"Rates "+chname,1600,1100);
    rcanv->Divide(3,2);
    rcanv->cd(1);
    gStyle->SetOptFit(111);
    gStyle->SetOptStat(1110);


    double mean = hRateCuts[ci]->GetMean();
    int imax = hRateCuts[ci]->GetMaximumBin();
    mean = hRateCuts[ci]->GetBinLowEdge(imax);
    double rms  = hRateCuts[ci]->GetRMS();
    double maxrate = hRateCuts[ci]->GetMaximum();
    double nentries  = hRateCuts[ci]->GetEntries();
    int fgaus = 0;
    TF1 *frate;
//   mean=3.3;
//   rms=1.0;
//   if (runNo>=219 && exprate>10) {mean = 500.; rms=7;} 
  
    if (mean>6.5)
    {
        frate = new TF1("gaus","gaus",0.,2000.);
        frate->SetParameter(1,1.1*mean);
        frate->SetParameter(2,rms);
        frate->SetParLimits(0,0.2*maxrate,3.*maxrate);
        frate->SetParLimits(1,0.5*mean,2.*mean);
        frate->SetParLimits(2,sqrt(mean)/50.,2.*sqrt(mean));
        fgaus=1;
    }
    
    else 
    {
        frate = new TF1("pois","[0]*TMath::Poisson(x,[1])",0,20); 
        frate->SetParameter(0,1);
        frate->SetParameter(1,mean);
    }
    frate->SetNpx(2000);
    hRateCuts[ci]->Sumw2();
    cout<<"Fit "<<channel<<"  "<<8<<endl;
    hRateCuts[ci]->Fit(frate,"WBN+","",mean-1.5*rms,mean+1.5*rms);
    cout<<"Fit "<<channel<<"  "<<9<<endl;
    hRateCuts[ci]->Fit(frate,"MEB+","",mean-1.5*rms,mean+1.5*rms);
    hRateCuts[ci]->SetLineColor(2);
    hRateCuts[ci]->Draw("");
//   frate->Draw("same");

    
    char txt2[200];

    TF1 *fbsl2;
    double mean2 = hRate[ci]->GetMean();
    double rms2  = hRate[ci]->GetRMS();
    int fgaus2 = 0;
// //   mean2 = 3.;
    if (longpulse)   /// Only when long pulse the hRate histo can be used
    {  
        hRate[ci]->Draw("");
        //rms=1.1;
        if (mean2>6.5)
        {
            fbsl2 = new TF1("gaus2","gaus",0.,2000.);
            fbsl2->SetParameter(1,mean2);
            fbsl2->SetParameter(2,rms2);
            fgaus2=1;
        }
        else 
        {
            fbsl2 = new TF1("pois2","[0]*TMath::Poisson(x,[1])",0,20); 
            fbsl2->SetParameter(0,1);
            fbsl2->SetParameter(1,mean2);
        }
        hRate[ci]->Sumw2();
        cout<<"Fit "<<10<<endl;
        hRate[ci]->Fit(fbsl2,"W+","same",1.6,20.);
        cout<<"Fit "<<11<<endl;
        hRate[ci]->Fit(fbsl2,"RE+","same",1.6,20.);
        fbsl2->Draw("same");  // no point to draw these
    }

  
    double crate =frate->GetParameter(1);

    TPad *subpad=(TPad*)rcanv->GetPrimitive(Form("%s_%d",rcanv->GetName(),1));
    subpad->Update();  /// VERY IMPORTANT!!!!! crashes without update!!!
    TPaveStats *st = (TPaveStats*)subpad->GetPrimitive("stats");
    // if (st==NULL) return -3243;
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetTextSize(0.03);
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetX2NDC(0.9); //new x start position
    st->SetY2NDC(0.9); //new x end position 
    st->SetTextSize(0.031);
    subpad->Modified();

    
    //   rcanv->Update();
    rcanv->cd(2);
    //   gStyle->SetOptStat(0);
    hRateEvolution[ci]->SetStats(0);  
    //hRateEvolutionCuts[ci]->SetStats(0);  
    double maxy = hRateEvolution[ci]->GetMaximum()*1.6;
     hRateEvolution[ci]->SetMaximum(maxy);
      hRateEvolution[ci]->Draw();
    hRateEvolutionCuts[ci]->SetLineColor(2);
    hRateEvolutionCuts[ci]->SetLineWidth(2);
    hRateEvolutionCuts[ci]->Draw("same"); 
  

    rcanv->cd(3);
    gPad->SetLogy();
    gStyle->SetOptFit(111);
    gStyle->SetOptStat(1110);
  
  //  hAMPL->SetStats(0);
  //  hAMPL->Draw();
    hAMPL[ci]->Draw("");
    hCutsAMPL[ci]->SetLineColor(2);
    hCutsAMPL[ci]->Draw("same");
    
    subpad=(TPad*)rcanv->GetPrimitive(Form("%s_%d",rcanv->GetName(),3));
    subpad->Update();  /// VERY IMPORTANT!!!!! crashes without update!!!
    st = (TPaveStats*)subpad->GetPrimitive("stats");
    // if (st==NULL) return -3243;
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetTextSize(0.03);
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetX2NDC(0.9); //new x start position
    st->SetY2NDC(0.9); //new x end position 
    st->SetTextSize(0.031);
    subpad->Modified();

//   rcanv->Update();
    rcanv->cd(6);
    gPad->SetLogy();
    //hCutsCH[ci]->Draw("");
    htotCH[ci]->Draw("");
    //hCutsCH[ci]->SetLineColor(2);
    
    subpad=(TPad*)rcanv->GetPrimitive(Form("%s_%d",rcanv->GetName(),6));
    subpad->Update();  /// VERY IMPORTANT!!!!! crashes without update!!!
    st = (TPaveStats*)subpad->GetPrimitive("stats");
    // if (st==NULL) return -3243;
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetTextSize(0.03);
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetX2NDC(0.9); //new x start position
    st->SetY2NDC(0.9); //new x end position 
    st->SetTextSize(0.031);
    
    subpad->Modified();

    rcanv->cd(4);   
    double ttot = tmax-0.;
    //exponential fit
    TF1* f1 = new TF1("f1","[0]*exp(-[1]*x)",0.,200.);

    f1->SetParLimits(1,1e-6,1e6);
    f1->SetParLimits(0,1e-6,1e6);

    f1->SetParameter(0,nevents/5.);
    f1->SetParameter(1,exprate);
    f1->SetLineColor(4);
    
    hDt[ci]->Sumw2(1);
    hCutsDt[ci]->Sumw2(1);

    
    double tbeam = 5e-6;

    cout<<"Fit "<<channel<<"  "<<12<<endl;
    hDt[ci]->Fit(f1,"+QB0","",3*tmax/100.,tmax);//crate[4]);
    
    double rates=1.*f1->GetParameter(1);
    char title[100];
    sprintf(title,"f(x) = %g e^{-%g x}",f1->GetParameter(0),f1->GetParameter(1));
    f1->SetTitle(title);
    
    gStyle->SetOptFit(111);
    gStyle->SetOptStat(110);
    cout<<"total rate = "<<crate<<" c/s \t calculated = "<<rates<<" c/s "<<endl;
    hDt[ci]->SetLineColor(4);
    hDt[ci]->SetMarkerColor(4);
    hDt[ci]->SetMarkerStyle(2);
    hDt[ci]->Draw("E0");
    hDt[ci]->SetMinimum(0.1);
    hDt[ci]->SetMaximum(hDt[ci]->GetMaximum()*2.);

    hCutsDt[ci]->Draw("same");

    TF1* f11 = new TF1("f11","[0]*exp(-[1]*x)",0.,tmax*2.);


    f1->Draw("same"); 
    //  f11->Draw("same"); 
    
    hCutsDt[ci]->SetLineColor(kRed+1);
    TF1* f12 = new TF1("f12","[0]*exp(-[1]*x)",0.,tmax*2);
    f12->SetParLimits(1,1e-6,1e6);
    f12->SetParameter(0,120.);
    f12->SetParameter(1,exprate);
    f12->SetParLimits(0,0.1e-6,1e6);
    f12->SetLineColor(kRed+2);
    cout<<"Fit "<<channel<<"  "<<13<<endl;
    hCutsDt[ci]->Fit(f12,"+QB0","",3*tmax/100.,tmax);//crate[4]);
    f12->Draw("same"); 
    
    TPad* pad = (TPad*) gPad;
    pad->SetLogy();
  
    subpad=(TPad*)rcanv->GetPrimitive(Form("%s_%d",rcanv->GetName(),4));
    subpad->Update();  /// VERY IMPORTANT!!!!! crashes without update!!!
    st = (TPaveStats*)subpad->GetPrimitive("stats");
    // if (st==NULL) return -3243;
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetTextSize(0.03);
    st->SetX1NDC(0.55); //new x start position
    st->SetY1NDC(0.65); //new x end position 
    st->SetX2NDC(0.9); //new x start position
    st->SetY2NDC(0.9); //new x end position 
    st->SetTextSize(0.031);
    
    subpad->Modified();

  
    cout << "ntrigsTotal = " << ntrigsTot[ci] << endl;
    double SF =  ntrigsTot[ci] / (f1->GetParameter(1));
    cout << "SF = " << SF << endl;
  
  
    //   rcanv->Update();
    rcanv->cd(5); 
    sprintf(txt2," #it{run:} #color[4]{#bf{%03d}} C%d \t #it{amplifier:} #color[4]{#bf{%s}} \t #it{R = } #color[4]{#bf{%s m}}",runNo,ci+1,oscsetup->Amplifier[ci],oscsetup->DetName[ci]); 
    TLatex *tex = new TLatex(0.02,0.95,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.065);
    tex->SetLineWidth(2);
    tex->SetTextAlign(12);
    tex->Draw();

//     sprintf(txt2," V_{m} = #color[4]{-%d}V, V_{d} = #color[4]{-%d}V, peak = #color[4]{%g}mV, n = #color[4]{%g}mV",vm,vd,fabs(Thresholds[ci]*1000.),fabs(peTh*1000.)); 
//     sprintf(txt2," V_{m} = #bf{-%d V}, V_{d} = #bf{-%d V}, peak = #bf{%g mV}, n = #bf{%g mV}",vm,vd,fabs(Thresholds[ci]*1000.),fabs(peTh*1000.)); 
    sprintf(txt2," V_{1} = #bf{-%f V},  V_{2} = #bf{-%f V}",oscsetup->V1[ci],oscsetup->V2[ci]); 
    tex = new TLatex(0.5,0.88,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
    sprintf(txt2," #it{Thresholds:}  peak = #bf{%g mV},  n = #bf{%g mV}",fabs(Thresholds[ci]*1000.),fabs(peTh*1000.)); 
    tex = new TLatex(0.5,0.80,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
    double calcr = frate->GetParameter(1);
    double calcrerr = frate->GetParError(1);
    if (fgaus)
    {
        calcr = frate->GetParameter(1);
        calcrerr = frate->GetParameter(2);
    }
    sprintf(txt2,"#LTr_{n}#GT = #bf{%3.0f #pm %3.1f} #frac{n}{#it{%gs}} = #bf{%3.0f #pm %3.1f} s^{-1}",calcr,calcrerr,timebinwidth,calcr/timebinwidth,calcrerr/timebinwidth);
    //  if (fgaus)
    //    sprintf(txt2,"#LTr_{n(#it{#Deltat=%gs})}#GT = #bf{%3.2f #pm %3.3f} #frac{n}{pulse}",timebinwidth,frate->GetParameter(1),frate->GetParameter(2));
    tex = new TLatex(0.5,0.7,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();

    
    sprintf(txt2,"#DeltaT fit, all events:");
    tex = new TLatex(0.5,0.6,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
    
    sprintf(txt2,"#LTr_{all events}#GT = #bf{%3.1f #pm %3.1f} #frac{counts}{sec}",f1->GetParameter(1),f1->GetParError(1));
    if (exprate<1)
        sprintf(txt2,"#LTr_{all events}#GT = #(){#bf{%3.1f #pm %3.2f}}#times10^{-3} #frac{counts}{sec}",f1->GetParameter(1)*1e3,f1->GetParError(1)*1e3);
    tex = new TLatex(0.5,0.55,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
    // runpar->grate=f1->GetParameter(1);
    // runpar->sgrate=f1->GetParError(1);
    
    sprintf(txt2,"#LTr_{n(#it{thr=%gmV})}#GT = #bf{%3.1f #pm %3.1f} #frac{counts}{sec}",peTh*mV,f12->GetParameter(1),f12->GetParError(1));
    if (exprate<1)
        sprintf(txt2,"#LTr_{n(#it{thr=%gmV})}#GT = #(){#bf{%3.1f #pm %3.1f}}#times10^{-3} #frac{counts}{sec}",peTh*mV,f12->GetParameter(1)*1000.,f12->GetParError(1)*1000.);
    tex = new TLatex(0.5,0.46,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
    // runpar->rate=f12->GetParameter(1);
    // runpar->srate=f12->GetParError(1);
    
    long double duration = (epochF-epochS)/3600. ;
    sprintf(txt2,"#it{#bf{%d sparks} were observed within #bf{%3.1LF} hours}",totNsparks,duration);
    if (totNsparks==0)
        sprintf(txt2,"#it{#bf{No spark} was observed within #bf{%3.1LF} hours}",duration);
    else if (totNsparks==1)
        sprintf(txt2,"#it{#bf{%d spark} was observed within #bf{%3.1LF} hours}",totNsparks,duration);

    tex = new TLatex(0.5,0.335,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();

    sprintf(txt2,"#it{Polya fits}");
    tex = new TLatex(0.5,0.25,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();

    sprintf(txt2,"#it{Ampl}: %s = #bf{%3.1f #pm %3.2f} mV , %s = #bf{%3.3g} ",fpolyaAmpl->GetParName(1),fpolyaAmpl->GetParameter(1)*1000./mV,fpolyaAmpl->GetParError(1)*1000./mV,fpolyaAmpl->GetParName(2),fpolyaAmpl->GetParameter(2));
    tex = new TLatex(0.5,0.18,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
    // runpar->ampl = fpolyaAmpl->GetParameter(1);
    // runpar->sampl = fpolyaAmpl->GetParameter(2);
    
    sprintf(txt2,"#it{Charge}: %s = #bf{%3.1f #pm %3.2f} , %s = #bf{%3.3g} ",fpolyaCH->GetParName(1),fpolyaCH->GetParameter(1),fpolyaCH->GetParError(1),fpolyaCH->GetParName(2),fpolyaCH->GetParameter(2));
    tex = new TLatex(0.5,0.1,txt2);
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->SetTextAlign(23);
    tex->Draw();
  
    // runpar->charge = fpolyaCH->GetParameter(1);
    // runpar->scharge = fpolyaCH->GetParameter(2);

    rcanv->Update();
  
    gSystem->ChangeDirectory(allplotdirname);
    rcanv->SaveAs(".png");
    rcanv->SaveAs(".pdf");

 
    TCanvas *canvrtev = new TCanvas("RateEvolution_"+chname,"Rate Evolution "+chname);
//     hRateEvolutionCuts[ci]->Draw(""); 
//     hSparkEvolution[ci]->Draw("sameE0");
      hRateEvolution[ci]->Draw("E0");
    hRateEvolutionCuts[ci]->Draw("SAMEE0"); 
    canvrtev->BuildLegend();
    canvrtev->Modified();
    canvrtev->Update();
    
    gSystem->ChangeDirectory(allplotdirname);
    canvrtev->SaveAs(".png");
    canvrtev->SaveAs(".pdf");
  
/// end of rate plots

// // // // 
 
    TCanvas *camplch = new TCanvas("CHovAmplitude"+chname,"Charge over Amplitude "+chname);
    gStyle->SetOptFit(111);
    gStyle->SetOptStat(1100);
    htotCHovAmplCuts[ci]->Sumw2();
    htotCHovAmplCuts[ci]->Draw();
    htotCHovAmpl[ci]->Draw("same");


    int maxbin = htotCHovAmpl[ci]->GetMaximumBin();
    double maxbinval = htotCHovAmpl[ci]->GetBinLowEdge(maxbin);
    TF1 *fgaus3 = new TF1("fgaus4","gaus",maxbin-100.,maxbin+100.);
    fgaus3->SetParameter(2,sqrt(maxbinval));
    fgaus3->SetParameter(0,htotCHovAmpl[ci]->GetMaximum());
    fgaus3->SetParameter(1,maxbinval);
    cout<<"Fit "<<channel<<"  "<<14<<endl;
    htotCHovAmpl[ci]->Fit(fgaus3,"M+","",maxbinval-250.0,maxbinval+250.0);
    // runpar->chovampl = fgaus3->GetParameter(1);
    // runpar->schovampl = fgaus3->GetParameter(2);


    gStyle->SetOptFit(111);

    TCanvas *sccanv = new TCanvas("SignalCuts"+chname,"Signal Cuts "+chname,1600,1100);
    sccanv->Divide(4,2);
 
    sccanv->cd(1);
    gStyle->SetOptStat(1100);
    htotCHovAmpl[ci]->Draw();
    htotCHovAmplCuts[ci]->Draw("same");

    sccanv->cd(2);
    hRateStructure[ci]->SetStats(0);
    hRateStructure[ci]->GetXaxis()->SetLabelSize(0.05);
    hRateStructure[ci]->GetXaxis()->SetTitleSize(0.05);
    hRateStructure[ci]->GetYaxis()->SetLabelSize(0.05);
    maxy = hRateStructure[ci]->GetMaximum()*1.2;
    hRateStructure[ci]->SetMaximum(maxy);
    hRateStructure[ci]->Draw();
    //hCutsRateStructure[ci]->SetLineColor(2);
    //hCutsRateStructure[ci]->Draw("same");
    hSRStriggerStruct[ci]->Draw("same");

    TCanvas * canvstruct = new TCanvas("RateStructures_"+chname,"Rate Structures "+chname) ;
    hRateStructure[ci]->Draw();
    //hCutsRateStructure[ci]->Draw("same");
    canvstruct->BuildLegend();
    canvstruct->Modified();
    canvstruct->Update();
    gSystem->ChangeDirectory(allplotdirname);
    canvstruct->SaveAs(".png");
    canvstruct->SaveAs(".pdf");

    
    sccanv->cd(6);
    gStyle->SetOptStat(1100);
    //hCutsRT[ci]->Draw("");
    hRT[ci]->Draw("");

 
    sccanv->cd(4);

    hPW[ci]->Draw("");

    sccanv->cd(3);
    //hCutsTOT[ci]->Draw("");
    hPD[ci]->Draw("");

    sccanv->cd(7);
    h2dAmplvsRT[ci]->Draw("colz");
    sccanv->cd(5);
    h2dtotCHvsAmpl[ci]->Draw("colz");
    sccanv->cd(8);
    h2dAmplvsPD[ci]->Draw("colz");

    sccanv->Modified();
    sccanv->Update();
    sccanv->SaveAs(".png");
    sccanv->SaveAs(".pdf");
 
    
    TCanvas *evcanv = new TCanvas("gainEvolution_"+chname,"Pulse size evolution "+chname,1600,1800);
    evcanv->Divide(1,2);
    evcanv->cd(1);
    h2dCutsAmplEvolution[ci]->Draw("colz");
    evcanv->cd(2);
    h2dCutsCHevolution[ci]->Draw("colz");
    evcanv->Modified();
    evcanv->Update();
    evcanv->SaveAs(".png");
    evcanv->SaveAs(".pdf");
   
    

    if (isMM)
    {
      TCanvas *champlcanv = new TCanvas("ChAmplCorrelations"+chname,"Charge - Amplitude Correlations "+chname,1600,1100);
      champlcanv->Divide(4,2);
  
      champlcanv->cd(1);
      gStyle->SetOptStat(1100);
  
      htotCHovAmpl[ci]->Draw();
      htotCHovAmplCuts[ci]->SetLineColor(2);
      htotCHovAmplCuts[ci]->Draw("same");
  
      champlcanv->cd(4+1);
      h2dtotCHvsAmpl[ci]->Draw();
      
      champlcanv->cd(2);
      gStyle->SetOptStat(1100);
  
      heCHovAmpl[ci]->Draw();

      champlcanv->cd(4+2);
      h2deCHvsAmpl[ci]->Draw();
      
      champlcanv->cd(3);
      gStyle->SetOptStat(1100);
  
      heCHovtotCh[ci]->Draw();

      champlcanv->cd(4+3);
      h2deCHvstotCh[ci]->Draw();
      
      champlcanv->cd(4);
      gStyle->SetOptStat(1100);
  
      hAmplovDampl[ci]->Draw();

      champlcanv->cd(4+4);
      h2dAmplvsDampl[ci]->Draw();
      

      champlcanv->Modified();
      champlcanv->Update();
      champlcanv->SaveAs(".png");
      champlcanv->SaveAs(".pdf");
      
      ///--------
      TCanvas *echichcanv = new TCanvas("eCHiCH"+chname,"e-Charge vs i-Charge "+chname);
      h2deCHvsiCH[ci]->Draw("colz");
      
      echichcanv->Modified();
      echichcanv->Update();
      echichcanv->SaveAs(".png");
      echichcanv->SaveAs(".pdf");
      
      
    }  /// end of MM only canvases
    
  }    /// end of ci loop for canvases  

//   TString chname = ctypes;
  TCanvas *xycanv = new TCanvas("XYmaps","XY Maps "+ctypes,1600,1100);
  xycanv->Divide(actch,2);
  int ci_shift=0;
  
  for (int ci=0; ci<4; ci++)
  {
     if (active[ci]==0) continue;
     

    cout<<BLUE<<"Entering plots directory "<<MAGENTA<< plotdirname <<BLUE<<" for channel "<<MAGENTA<<ci+1<<endlr;
    gSystem->ChangeDirectory(plotdirname);
    cout<<GREEN<<"done"<<endlr;

    xycanv->cd(ci_shift+1);
    h2dXY[ci]->Draw();
    xycanv->cd(ci_shift+actch+1);
    h2dXYcuts[ci]->Draw();
     
     ci_shift++;
     
     xycanv->Modified();
     xycanv->Update();
     xycanv->SaveAs(".png");
     xycanv->SaveAs(".pdf");
     
  }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// end of ci loop
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 ofile->Write("",TObject::kOverwrite);
 gSystem->ChangeDirectory(tmpdir);
 cout<<BLUE<<"End of script!"<<endl;
 return 0; 
//// turn the warnings back on
}
#if defined __SUNPRO_CC
#   pragma enable_warn
#elif defined _MSC_VER
#   pragma warning(pop)
#endif
