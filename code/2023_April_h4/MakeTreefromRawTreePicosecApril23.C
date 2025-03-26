#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <vector>
#include <iomanip>
// #include "MyFunctions.h"
#include "MyFunctions_2023_April.C"
// #include "tracks.C"

// void searcheandreplace(char *fnames)
// {
//     char tmp[1000];
//     int nchar=0;
//     do {
//       tmp[nchar]=fnames[nchar];
//       nchar++;
//     } while (fnames[nchar-1]!='\0'); 
//     
//     for (int i=0; i<nchar-1; i++)
//     {
//        if (fnames[i]==' ')
//        {
//           fnames[i]='\\';
//           for (int j=i+1;j<nchar+1;j++)
// 	    fnames[j]=tmp[j-1];
// 	cout<<fnames<<endl;
// 	  searcheandreplace(&fnames[i+2]);
//        } 
// 
//     }  
// }  
// // void replaceEOL(char *fnames)
// // {
// //     int nchar=0;
// //     while (fnames[nchar]!='\n')
// //     {
// //       nchar++;
// //     }
// //     fnames[nchar]='\0';
// //     
// // }  
// // void replaceEOLT(char *fnames)
// // {
// //     int nchar=0;
// //     while (fnames[nchar]!='\n')
// //     {
// //       nchar++;
// //     }
// //     fnames[nchar]='\t';
// // //     fnames[nchar+1]='\0';
// //     
// // }  

// double Summa(int npoints,double *X,double *Y,double mean, double sigma)
// {
//    double sum=0.;
//    for(int i=0;i<npoints;i++)
//    {
//        if (X[i]>mean-3*sigma)
//        {
// 	 sum+=Y[i];
// //          cout<<"["<<i<<" | "<<X[i]<<"] - "; 
//        }
//        if (X[i]>mean+4*sigma)
// 	 break;
//    }
//    return (sum);
// }


int MakeTreefromRawTreePicosecApril23(int runNo=1, int poolNo=1,  string filetype = "")  ///put trigger = ChannelNo in case of external trigger
{
  // gROOT->LoadMacro("MyFunctions.C");
  gROOT->LoadMacro("MyFunctions_2023_April.C");
  gROOT->LoadMacro("tracks.C");
  
  
  const double mV =1000.;
  
  
  int vm=275;
  int vd=500;

  float bthick,dgap;
  int angle=0;
  
  char command[5000];
  char fname[4][5000];
  char fname2[5000],ofilename[5000],ifilename[5000];
  char fnametmp[2000];
  vector<string> fnames[4];

  char rtype[1000];
  strcpy(rtype,RTYPE);

  char basedirname[1000];
  char datadirname[1000];
  char workdirname[1000];
  char outdirname[1000];
  sprintf(basedirname,"%s",BASEDIRNAME);
  sprintf(datadirname,"%s",DATADIRNAME);
  sprintf(workdirname,"%s",WORKDIR);
  sprintf(outdirname,"%s",OUTDIRNAME);


  OSCSETUP* oscsetup = new OSCSETUP;
  
  if (runNo>=MINRUN && runNo<=MAXRUN)
  {
    char afile[500]; 
    sprintf(afile,"%s/tmpfile.tmp",workdirname);
    FILE *ftmp=fopen(afile,"w");
   
    if (ftmp == NULL)
    {
      cout<<afile<<" can not be created. Probably the directory '"<<workdirname<<"' does not exit. Exiting..."<<endl;
      exit (-12);
    }
    fclose(ftmp);
  
    if (poolNo==0)
      sprintf(command,"cd %s\nls -d Run%03d-GDD*%s_%s*.root > %s 2>/dev/null",datadirname,runNo,filetype.c_str(),rtype,afile);
    else
      sprintf(command,"cd %s\nls -d Run%03d-Pool%d*%s_%s*.root > %s 2>/dev/null",datadirname,runNo,poolNo,filetype.c_str(),rtype,afile);

    int tst=system(command);
   cout<<command<<endl<<"returned: "<<tst<<endl;
    if (tst !=0)
    {
      //cout<<command<<endl<<"returned: "<<tst<<endl;
      cout<<"Probably the tree of the run "<<runNo<<" in Pool "<<poolNo<<" was not found in directory "<<datadirname<<endl<<"Exiting..."<<endl;
      return tst;
    }
    ftmp=fopen(afile,"r");
    if (fgets(fnametmp,200,ftmp) == NULL)
    {
      cout<<"Failed to read the filename for Pool "<<poolNo<<" run "<<runNo<<" at "<<datadirname<<endl<<"Exiting..."<<endl;
      return -2;
    }

    fclose(ftmp);
    sprintf(command,"rm %s\n",afile);
    tst=system(command);

    int rtmp = runNo, dtmp = 1;
    char ftypetmp[500];
    strcpy(ftypetmp,"");
    filetype.assign(ftypetmp);
//     int stst = sscanf(fnametmp,"%3d-%02d-%d-%d_%30[^ /,\n\t]",&rtmp,&dtmp,&vm,&vd,ftypetmp);
    cout <<BLUE<<"___________________________________________________\n\nReading osciloscope setup for Run "<<runNo<<" , Pool "<<poolNo<<endlr;

    int stst = GetOsciloscopeSetup(runNo, poolNo, basedirname, oscsetup);
    if (stst == 0)
    {
        cout<<RED<< "Run "<<runNo<<" , Pool "<<poolNo<<" was not found in OsciloscopeSetup.txt"<<endl;
        cout<<"Exiting..."<<endlr;
        return -10;
    }
    else if (stst == -1)
    {
        cout<<RED<<"OsciloscopeSetup.txt was not found in "<<basedirname<<endl;
        cout<<"Exiting..."<<endlr;
        return -1;
    }
    else if (stst == -11)
    {
        cout<<RED<<"Error while reading OsciloscopeSetup.txt. File may be corrupted or miss-written"<<endl;
        cout<<"Exiting..."<<endlr;
        return -11;
    }
//     cout<<"V12 = "<<oscsetup->V1[1]<<endl;
    
  }
  else
  {
    cout<<"Available runs: "<<MINRUN<<" - "<<MAXRUN<<endl;
    return (-2);
  }
  
  const int trigger = oscsetup->srsCh;
  cout<<"SRS trigger at channel "<<trigger<<endl<<endlr;
  if (trigger>4){
    cout<<RED<<"Wrong trigger channel : "<<trigger<<endlr;
    return -15;
  }
  
  
 // int amplNo=oscsetup->AmplifierNo[1];
  //char amplifier[20];

  // WhichAmplifier(amplNo, amplifier);

  
  
  // return -33;
  const char *ftype = filetype.c_str();  /// add here any directory supplement
  
  char photocathode[100];
  strcpy(photocathode,ftype);
  
  char ofname[1000];
  strcpy (ofname,fnametmp);
  char *pch;
  pch = strstr(ofname,"raw_tree.root");
//   sprintf(ofname,"%s_tree.root",pch);
  strcpy (pch,"_tree.root");
#ifdef DEBUG
  strcpy (pch,"_DEBUGtree.root");
#endif
  cout<<"output filename = "<<ofname<<endl<<endl;;

  
  char mesh[1000];
  sprintf(mesh,"%d",vm);
  char drift[1000];
  sprintf(drift,"%d",vd);

  int drft;
  sscanf(drift,"%d",&drft);

  const int MaxFiles = MAX_N_FILES;
  const int MAXSEG = 100000;
  const int FRAMESIZE = 2000000;
  
  double *ptime;
  ptime = new double[FRAMESIZE];


  cout <<RED<<"********************\n\n\n Now starting processing "<<fnametmp<<"..."<<endl<< endlr;  

  replaceEOL(fnametmp);
  sprintf(ifilename,"%s/%s",datadirname,fnametmp);
  TFile *infile = new TFile(ifilename);
  if (!infile->IsOpen()){
    cout<<"Failed to open \""<<ifilename<<"\"\ncorresponding to the run "<<runNo<<"\nExiting..."<<endl;
    return (-1);
  }
  
  int spoints[]={0,0,0,0};
  int evNo=0;
  double Dt[4]={0.,0.};	
  double t0[]={0.,0.,0.,0.};
  double gain[]={0.,0.,0.,0.};
  double offset[]={0.,0.,0.,0.};
  double rmax[]={0.,0.,0.,0.};
  double rmin[]={0.,0.,0.,0.};
  double range[]={0.,0.,0.,0.};
  double dt, odt;
  unsigned long long int epoch;
  unsigned long long int nn;
  int itrigger;//time of the trigger
  char date[4][50];
  double *amplC[4];
  int active[] = {0,0,0,0};
//   for (int i=0; i<2; i++)
  
  char treename[200];
  sprintf(treename,"TreeWithRawData");
  char treetitle[200];
  sprintf(treetitle,"Raw data tree");

  TTree *intree = (TTree*)infile->Get(treename);
  TBranch* branch;

  int inpoints;
  branch = intree->GetBranch("npoints");
  branch->SetAddress(&inpoints);
  int eventNo;
  branch = intree->GetBranch("eventNo");
  branch->SetAddress(&eventNo);
  
  branch = intree->GetBranch("t0");
  branch->SetAddress(t0);

  branch = intree->GetBranch("gain");
  branch->SetAddress(gain);
  branch = intree->GetBranch("offset");
  branch->SetAddress(offset);
  branch = intree->GetBranch("rmax");
  branch->SetAddress(rmax);
  branch = intree->GetBranch("rmin");
  branch->SetAddress(rmin);

  double idt;
  branch = intree->GetBranch("dt");
  branch->SetAddress(&idt);
  
  unsigned long long int iepoch=0;
  branch = intree->GetBranch("epoch");
  branch->SetAddress(&iepoch);
  unsigned long long int inn;
  branch = intree->GetBranch("nn");
  branch->SetAddress(&inn);
  char idate[50];
  branch = intree->GetBranch("date");
  branch->SetAddress(&idate);

  for (int i=0;i<4; i++)
  {
    TString channel = TString::Itoa(i+1,10);
    TString bname = "npoints"+channel;
    TString btype = bname+"/I";
    branch = intree->GetBranch(bname);
    if (branch==NULL) 
      continue;
    active[i]=1;
    branch->SetAddress(&spoints[i]);
    cout<<GREEN<<"Maximum npoints of channel "<<i+1<<" = "<<intree->GetMaximum(bname)<<endlr;
    bname = "amplC"+channel;
    cout <<MAGENTA<<"found tree branch "<< bname << endl;
    btype = bname+"[npoints"+channel+"]/D";
    branch = intree->GetBranch(bname);
    amplC[i] = new double[FRAMESIZE]; 
    branch->SetAddress(&amplC[i][0]);
  }

//   intree->SetCacheSize(200000000);
//   long int bufsize = intree->GetCacheSize();
//   cout<<"Cache size = "<<bufsize<<endl;
// //   intree->InitializeBranchLists(true);
//   intree->SetBasketSize("*",32000);
//   intree->AddBranchToCache("*");
//   intree->LoadBaskets();
//   bufsize = intree->GetCacheSize();
//   cout<<"Cache size = "<<bufsize<<endl;
  const int nevents = intree->GetEntries();

//   intree->AddBranchToCache("*");
  //for (int i=0;i<nevents;i++)

  
  const int nadiv = 8;//previous 2
  
  int actch=0;
  int framesize=0;
  double bslmax[] = {1.,1.,1.,1.};
  double bslmin[] = {-1.,-1.,-1.,-1.};
  double meanoffset[4];
  double vpb[4];
  for (int i=0;i<4; i++)
  {
    double max=0.;
    double min=0.;
    if ((active[i]) )
    {
      actch++;
      intree->GetEntry(nevents/2);
      
      framesize = spoints[i];
      cout<<CYAN<<"\nFramesize = "<<framesize<<" points"<<"\t";
      cout<<"Gain = "<<gain[i]<<" V"<<"\t";
      cout<<"Offset = "<<offset[i]<<" V"<<endl;
      cout<<"Range Max = "<< rmax[i]*mV<<" mV";
      cout<<"\t\tRange Min = "<< rmin[i]*mV<<"\t";
      range[i]=rmax[i]-rmin[i];
      cout<<"Full range = "<<range[i]*mV<<" mV  -->  " << range[i]/8 *mV <<" mV / div  --> "<<range[i]/256*mV<<" mV / bit"<<endl;
      bslmax[i] = rmax[i];
      bslmin[i] = rmin[i] + range[i]*((8.-nadiv)/8.);
//cout<<"Corresponding to "<<endl<<endl<<endl;
//       break;  /// only the first active channel is used!!!!!!
       meanoffset[i]= (rmax[i] + rmin[i])/2. ;
      cout<<"Mean offset C"<<i+1<<" = "<< meanoffset[i]*mV <<" mV"<<endl;
      vpb[i] = range[i]/256.;
      cout<<"vpd = "<<vpb[i] <<endl;

    }
  }
  dt = idt*1e9;
  if (actch==0) 
  {
    cout<<"No active channel found, or the only active channel is declared as external trigger !!! Exiting..."<<endl;
    return -12;
  }
  //   drft=651;
  int runb = 0;
  
  char rfname[200];
  char poolid[20];
  if (poolNo==0)
    sprintf(poolid,"GDD");
  else
    sprintf(poolid,"Pool%d",poolNo);
  char runid[20];
  sprintf(runid,"Run%03d",runNo);
  char ctrigger[20];
  sprintf(ctrigger,"%d",trigger);
  cout<<RED<<"Here I am"<<endlr;

  TString rtypes(runid);
  rtypes+="-";
  //rtypes+=amplifier;
  //rtypes+="-";
  rtypes+=rtype;
  
  //char catt[30];
  
  if (runb)
    rtypes+="_b";

//   TString ctypes(runid);
//   ctypes+="-";
//   ctypes+=poolid;
//   ctypes+="-";
// //   ctypes+=mesh;
// //   ctypes+="-";
// //   ctypes+=drift;
// //   ctypes+="-";
// //   ctypes+=ctrigger;
// //   ctypes+="-";
//   ctypes+=ftype;
// 
  
  int maxpoints=0;

  
  double XX[4][200];
  double YY[4][200];
  double refP1[200][3];
  double refP2[200][3];
  int eventTracks;

  double *amplCo[4];
  for (int i=0; i<4; i++)
  {
    if (trigger!=i+1 && active[i])
      amplCo[i] = new double[FRAMESIZE]; 
  }
  
///  creating output file
  sprintf(ofilename,"%s/%s",outdirname,ofname);
  TFile *ofile = new TFile(ofilename,"RECREATE");

/// Creating the osciloscope info tree
  
  TTree *infotree;
  
  sprintf(treename,"OsciloscopeSetup");
  sprintf(treetitle,"Osciloscope Setup");
  infotree = new TTree(treename,treetitle);
   
  int srsChNo=-1;

  infotree->Branch("srsCh", &srsChNo, "srsCh/I");
  float V1, V2, ZZ;
  infotree->Branch("V1", &V1, "V1/F");
  infotree->Branch("V2", &V2, "V2/F");
  infotree->Branch("Z", &ZZ, "Z/F");
  char DetName[20];
  infotree->Branch("DetName",DetName , "DetName[10]/C");
  char Photocathode[20];
  infotree->Branch("Photocathode",Photocathode , "Photocathode[12]/C");
  char Amplifier[20];
  infotree->Branch("Amplifier",Amplifier , "Amplifier[12]/C");
  int AmplifierNo;
  infotree->Branch("AmplifierNo", &AmplifierNo, "AmplifierNo/I");

  srsChNo=oscsetup->srsCh;
  for (int i=0;i<4;i++)
  {
      V1=oscsetup->V1[i];
      V2=oscsetup->V2[i];
      ZZ=oscsetup->Z[i];
      strcpy(DetName,oscsetup->DetName[i]);
      strcpy(Photocathode,oscsetup->Photocathode[i]);
      strcpy(Amplifier,oscsetup->Amplifier[i]);
      AmplifierNo=oscsetup->AmplifierNo[i];
      infotree->Fill();
  }
///___________________________________________________

/// creating output  tree
  ofile->cd();
  TTree *outtree;
  int srsNo=-1;
  sprintf(treename,"RawDataTree");
  sprintf(treetitle,"Data tree");
  outtree = new TTree(treename,treetitle);

  outtree->Branch("eventNo", &evNo, "eventNo/I");
  outtree->Branch("srsNo", &srsNo, "srsNo/I");
  outtree->Branch("dt", &odt, "dt/D");
  outtree->Branch("epoch", &epoch, "epoch/l");
  outtree->Branch("nn", &nn, "nn/l");
  // tree->Branch("date", &date, "date/C");
  outtree->Branch("t0",t0, "t0[4]/D");
  outtree->Branch("itrigger",&itrigger,"itrigger/I");
  double ttrig;
  outtree->Branch("ttrig",&ttrig,"ttrig/D");
  int fitstatus1[]={0,0,0,0}, fitstatus2[]={0,0,0,0};
  outtree->Branch("fitstatus1",fitstatus1,"fitstatus1[4]/I");
  outtree->Branch("fitstatus2",fitstatus2,"fitstatus2[4]/I");
  double frangemax[4],frangemin[4];
  outtree->Branch("rmax",frangemax,"rmax[4]/D");
  outtree->Branch("rmin",frangemin,"rmin[4]/D");

  int trackOK;
  outtree->Branch("trackOK",&trackOK,"trackOK/I");
  outtree->Branch("eventTracks",&eventTracks,"eventTracks/I");
  double chi2track[200];
  outtree->Branch("chi2track",chi2track,"chi2track[eventTracks]/D");
  outtree->Branch("refP1",&refP1[0][0],"refP1[eventTracks][3]/D");
  outtree->Branch("refP2",&refP2[0][0],"refP2[eventTracks][3]/D");
  double disttonextcluster[200][6];
  outtree->Branch("disttonextcluster",disttonextcluster,"disttonextcluster[eventTracks][6]/D");
  double totchargenextcluster[200][6];
  outtree->Branch("totchargenextcluster",totchargenextcluster,"totchargenextcluster[eventTracks][6]/D");
  
  outtree->Branch("sumpoints", &maxpoints, "sumpoints/I");
  for (int i=0; i<4; i++)
  {
    if (trigger!=i+1 && active[i])
    {
       TString channel = TString::Itoa(i+1,10);
       TString bname = "amplC"+channel;
       TString btype = bname+"[sumpoints]/D";
       outtree->Branch(bname, &amplCo[i][0], btype);

       TString bname2 = "hitX_C"+channel;
       TString btype2 = bname2+"[eventTracks]/D";
       outtree->Branch(bname2, &XX[i][0], btype2);

       TString bname3 = "hitY_C"+channel;
       TString btype3 = bname3+"[eventTracks]/D";
       outtree->Branch(bname3, &YY[i][0], btype3);
    }
  }
//   outtree->Branch("amplSum2", amplSum[1], "amplSum2[sumpoints]/D");
  double bslC[]={0,0,0,0};
  double rmsC[]={0.1,0.1,0.1,0.1};
  outtree->Branch("bslC", bslC, "bslC[4]/D");
  outtree->Branch("rmsC", rmsC, "rmsC[4]/D");
  

  cout<<GREEN<<"\n\n-->  Tree is ready!"<<endlr;

  double binscale =1.;
  if (framesize < 1000)
    binscale = 1.0;
  int nbins = 256 / 8 * nadiv -1;
  
  char htitle[100];

  TH1F *hbsl[4], *hbsl0[4];
  TH1D *hbslall[4], *hbsl0all[4];

  for (int ci=0;ci<4;ci++)
  {
    if (!active[ci]) continue;
//     if (trigger==ci+1) continue;
    sprintf(fname2,"Signal_bsl_%d",ci+1);
    hbsl[ci] = new TH1F(fname2,fname2,nbins,bslmin[ci],bslmax[ci]);
    hbsl[ci]->GetXaxis()->SetTitle("[V]");
//     cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;
    sprintf(fname2,"Signal_bsl_tot_%d",ci+1);
    sprintf(htitle,"Signal baseline from all events");
    hbslall[ci] = new TH1D(fname2,htitle,nbins,bslmin[ci],bslmax[ci]);
    hbslall[ci]->GetXaxis()->SetTitle("[V]");
//     cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;

    sprintf(fname2,"Subtracted_bsl_%d",ci+1);
    hbsl0[ci] = new TH1F(fname2,fname2,(int)nbins*binscale,-range[ci]/8.,range[ci]/8.);
    hbsl0[ci]->GetXaxis()->SetTitle("[V]");
//     cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;
    sprintf(fname2,"Subtracted_bsl_tot_%d",ci+1);
    hbsl0all[ci] = new TH1D(fname2,fname2,(int)nbins*binscale,-range[ci]/8.,range[ci]/8.);
    hbsl0all[ci]->GetXaxis()->SetTitle("[V]");
//     cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;
  }


//   double *amplSum; //in this case (lina4) is just ampl - baseline 
//     amplSum = new double[FRAMESIZE];
  
  double *tbase = new double[FRAMESIZE];
  for (int i=0;i<FRAMESIZE;i++)
    tbase[i]=i*dt;
    
  double *EVTIME[4], *BSL[4], *RMS[4], *sBSL[4], *sRMS[4];
  
  
  for (int i=0; i<4; i++)
  {
      if (active[i])
      {
        EVTIME[i] = new double[nevents];  
        BSL[i] = new double[nevents];  
        RMS[i] = new double[nevents];  
        sBSL[i] = new double[nevents];  
        sRMS[i] = new double[nevents];  
      }
  }
  
//! **************************************************************  
/// read the tracking tree here
  int trackingFileOK = 0;
  char trackfilename[1000], tracktreename[1000];
  sprintf(trackfilename,"%s/anaRun%04d.root",TRACKDIRNAME,runNo);
  TFile *trackfile = new TFile(trackfilename,"READ");
  sprintf(tracktreename,"tracks");
  sprintf(treetitle,"Tracks");

  TTree *tracktree = (TTree*)trackfile->Get(tracktreename);
  
  Int_t           srstriggerctr=-1;
  Int_t           srstimestamp=-1;
  Int_t           ntracks=-1;
  Int_t           tracknumber=-1;
  Double_t        trackchi2=-1;
  std::vector<std::vector<double> > *hits;
  std::vector<double>  *distnextcluster;
  std::vector<double>  *totchanextcluster;
  hits = 0;
  distnextcluster = 0;
  totchanextcluster = 0;
  int trackTreeEntries=0;
  if (tracktree==NULL)
  {
      cout<<RED<<"Failed to read tracking file "<<trackfilename<<".\n\nThere will be no tracking information for the analysis for the run "<<runNo<<endlr;
      trackingFileOK=-1;
  }
  else if (tracktree->GetEntries()<1)
  {
      cout<<RED<<"The tracking tree in file "<<trackfilename<<"\nhas no entries!\nPlease make sure that it has been created correctly! \n\nThere will be no tracking information for the analysis for the run "<<runNo<<endlr;
      trackingFileOK=-2;
  }
  else
  {
    cout<<BLUE<<"Reading tracking tree in file "<<trackfilename<<"\nwith "<<tracktree->GetEntries()<<" entries"<<endlr;
    TBranch* b_srstriggerctr;
    tracktree->SetBranchAddress("srstriggerctr", &srstriggerctr, &b_srstriggerctr);
    TBranch* b_srstimestamp;
    tracktree->SetBranchAddress("srstimestamp", &srstimestamp, &b_srstimestamp);
    TBranch* b_ntracks;
    tracktree->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
    TBranch* b_tracknumber;
    tracktree->SetBranchAddress("tracknumber", &tracknumber, &b_tracknumber);
    TBranch* b_trackchi2;
    tracktree->SetBranchAddress("trackchi2", &trackchi2, &b_trackchi2);
    TBranch* b_hits;
    tracktree->SetBranchAddress("hits", &hits, &b_hits);
    TBranch* b_distnextcluster;
    tracktree->SetBranchAddress("distnextcluster", &distnextcluster, &b_distnextcluster);
    TBranch* b_totchanextcluster;
// cout<<" Bafbvsdfbsdv sfd sfd"<<endl;
    tracktree->SetBranchAddress("totchanextcluster", &totchanextcluster, &b_totchanextcluster);
// cout<<" kljrhglfdshsdfljgdsdv sfd sfd"<<endl;
    trackingFileOK=1;
    trackTreeEntries=tracktree->GetEntries();
  }
  if(trackingFileOK>0)
    tracktree->GetEntry(tracktree->GetEntries()-1);
    const int maxTrackID = (srstriggerctr>0?srstriggerctr:-1);

  TRACKDATA *trackcollection;

  if(trackingFileOK>0)
    trackcollection = new TRACKDATA[maxTrackID+1];
  else
    trackcollection = new TRACKDATA[1];
//   TRACKDATA trackcollection[maxTrackID];
     cout <<MAGENTA<<"Allocated "<<sizeof(TRACKDATA)*(maxTrackID+1)/1024/1024<<" MB for the tracking data arrays"<<endlr;
 
//   cout <<"Allocated "<<sizeof(trackcollection) / 1024  <<" kB"<<endl;

//   return 4;
   int detPositions=0;  
  
   for (int i=0;i<maxTrackID+1; i++)
   {
       trackcollection[i].ntracks=0;
   }
   for (int i=0;i<trackTreeEntries; i++)
   {
       tracktree->GetEntry(i);
       int evID =  srstriggerctr;
       trackcollection[evID].ntracks = ntracks;
       trackcollection[evID].trackchi2[tracknumber]=trackchi2;
       trackcollection[evID].srstriggerctr[tracknumber]=srstriggerctr;
       trackcollection[evID].srstimestamp[tracknumber]=srstimestamp;
       for (int j=0;j<6;j++)
       {
           trackcollection[evID].distnextcluster[tracknumber][j] = distnextcluster->at(j);
           trackcollection[evID].totchanextcluster[tracknumber][j] = totchanextcluster->at(j);
       }
       detPositions= hits->size();
       for(int counter=0; counter<detPositions ; counter++)
       {
//                     cout<<"Z of position "<<counter<<" = "<<hits->at(counter).at(2)<<endl;
          trackcollection[evID].hits[tracknumber][counter][0] = hits->at(counter).at(0);   
          trackcollection[evID].hits[tracknumber][counter][1] = hits->at(counter).at(1);   
          trackcollection[evID].hits[tracknumber][counter][2] = hits->at(counter).at(2);   
          
       }

   }
  
   int cposition[4]={-1,-1,-1,-1};
   int agoodevent=0;
   for (int i=0;i<trackTreeEntries; i++)
   {
       if (trackcollection[i].ntracks>0)
       {
           agoodevent=i;
           break;
       }
   }
   int p1=0, p2=2;  
     
   for (int i=0;i<4;i++)
   {
       double z=oscsetup->Z[i];
       for (int j=0;j<detPositions;j++)
       {
//            cout<<"   "<<trackcollection[agoodevent].hits[0][j][2];
           if (z<0) p1 = j;
           if (z>0) p2 = j;
           if (z==trackcollection[agoodevent].hits[0][j][2])
           {
              cposition[i]=j;
              break;
           }
       }
       cout<<(cposition[i]>=0?GREEN:RED)<<"Channel C"<<i+1<<" at position "<<cposition[i] <<endlr;
   }
   cout<<BLUE<<"Chosen reference positions = "<<p1<<" & "<<p2<<endlr;
   
   if (trackingFileOK>0)
       trackfile->Close();

   for (int i=0;i<100 && i<maxTrackID+1; i+=100)
   {
       if (trackcollection[i].ntracks<1) continue;
       cout<<"Event "<<setw(6)<<i<<" positions:";
       for (int j=0;j<detPositions;j++)
           cout<<"  "<<trackcollection[i].hits[0][j][2];

       cout<<endl;
   }

  
//    for (int i=0;i<10&&i<maxTrackID+1; i++)
//    {
//       if (trackcollection[i].ntracks>1)
//       {
//           cout<<" event "<< i <<" ("<<trackcollection[i].srstriggerctr[0]<<" . "<<trackcollection[i].srstriggerctr[1]<<") has "<<trackcollection[i].ntracks<<" tracks"<<endl;
//           cout<<"chi2[0] = "<<trackcollection[i].trackchi2[0]<<" chi2[1] = "<<trackcollection[i].trackchi2[1]<<endl;
//       }
//       else if (trackcollection[i].ntracks==0)
//           cout<<RED<<" event "<< i <<" ("<<trackcollection[i].srstriggerctr[0]<<") has NO tracks"<<endlr;
// 
//    }
 
///___________________________________________________  end of tracking tree reading
   
    ofile->cd();
  

#ifdef DEBUG 
//    TCanvas *cbsl = new TCanvas("baseline"+rtypes,"baseline "+rtypes,1600*4/4,1200*3/4);
//     if (trigger>0) actch-=1;
    TCanvas *cbsl = new TCanvas("baseline"+rtypes,"baseline "+rtypes,1600*4/4,800*actch*2/4);
    cout<<BLUE<<"-->  Atchu = "<<actch<<endlr;
    cbsl->Divide(2,actch);
    int linep = 0;
#endif  

    
  double tshift[2]={0.,0.};

  //   int dd,mm,yy,hour,min,sec;
//   char month[20];
  int nseg;

  double bsltmp=0.;
  double rmstmp=0.;
  
  TF1* fbsl[4];
  TF1* fbslall[4];  
  for (int i=0; i<4; i++)
  {
      sprintf(fname2,"gausC%d",i);
      fbsl[i] = new TF1(fname2,"gaus",-1.,1.);
      fbsl[i]->SetNpx(2000);
      sprintf(fname2,"gausC%dtot",i);
      fbslall[i] = new TF1(fname2,"gaus",-1.,1.);
      fbslall[i]->SetNpx(2000);
  }
//   sprintf(fname2,"gaus0");
//   TF1* fbsl0 = new TF1(fname2,"gaus",-1.,1.);
//   fbsl0->SetNpx(2000);

//   sprintf(fname2,"gaus0tot");
//   TF1* fbsl0all = new TF1(fname2,"gaus",-1.,1.);
//   fbsl0all->SetNpx(2000);

  gStyle->SetOptFit(11);
  gStyle->SetOptStat(1100);

  int ngpoints[] = {0,0,0,0};
  int eventS = 0;
  int eventF = nevents;
// #ifdef DEBUG 
//     eventS = 16438;
//     cout<<"Give event No to debug: ";
//     cin>>eventS;
//     eventF = eventS;
// #endif
    

  intree->GetEntry(0);
  intree->GetEntry(nevents-1);

// #ifndef DEBUG  
   cout<<"Loading tree baskets in virtual memory"<<endl;
//   int nbaskets=intree->LoadBaskets(4000000000);
//   cout<<"Loaded "<<nbaskets<<" baskets "<<endl;
// #endif
  
  int overflow=0;
  int lastevent=-1;
  double tdt = -1.;  
  int tindex=0;
//   eventF =1000;
  for (int nEv=eventS; nEv<nevents && nEv<eventF+1; nEv+=1)
  { 
#ifdef DEBUG 
//     eventS = 16438;
    cout<<"Give event No to debug: ";
    cin>>nEv;
    if (nEv<0) 
    {
        cbsl->cd(0);
        return (-1);
    }
    linep = 0;
//     eventF = eventS;
#endif

    
    intree->GetEntry(nEv);
       //cout<<"Event"<<endl;

    odt=idt*1e+9;
    evNo = eventNo;

    if (tdt != idt)
    {
      for (int k=0; k<inpoints+10; k++ ) /// time array
        ptime[k]=k*idt;// tmpx*1e9;
      dt = idt*1e9;  //ns
      tdt=idt;
    }
    epoch=iepoch;
    nn=inn;
    double  epochF = (1.*iepoch + inn*1e-9);
  
//! check how many tracks are reccorded for this event
    
    bool trigExist=false;
    bool validSRS=true;
    for (int ci=0; ci<4;ci++)
    {
       //cout<<"Searching for trigger channel, c"<<ci+1<<endl;
      if (!active[ci]) continue;
      
      itrigger = -100;
      

      if (trigger == ci+1)   /// calculate the SRS event number here
      {

        for (int k=0; k<spoints[ci] - (12.5/dt); k++ )  ///processing trigger (ci=0)
        {
            if(amplC[ci][k]<-0.4)    /// in Oct 2022 the max amplitude was -1V 
            {
                k+= (int) (12.5/dt);
                itrigger =k;
                break;
            }
        }
        ttrig = itrigger * dt;
        trigExist=true;

  //cout<<"itrigger = "<<itrigger<<" , ttrig = "<<ttrig<<" ns"<<endl;

         srsNo=0;
         int step =(int) (25./dt);
         if (step<=0) return -44;
           if (itrigger + 16*step>=spoints[ci]) validSRS = false;
         int bit=1;
         for (int j=16; j>0; j--) 
         {
             int k = itrigger + j*step;
             if (amplC[ci][k]<-0.4)
                 srsNo+=bit;
             bit*=2;
         }
         if (!validSRS){ 
           srsNo = -1;
//            cout<<MAGENTA<<"Frame size too short to fit in the SRS trigger for event "<< evNo<< " (nEv = "<<nEv<<"). Set to -1"<<endlr;
         }
         else {
           if (srsNo+(overflow)*65536 < lastevent) overflow++;
           srsNo += (overflow) * 65536;
           lastevent=srsNo;
         }
//          cout<<"SRS No = "<<srsNo<<" for event " <<evNo<<endl;    
      }

      frangemax[ci] = rmax[ci];
      frangemin[ci] = rmin[ci];


      hbsl[ci]->Reset("");

//       return 3;

      long double sumsig=0.;
      long double diff=0.;
    if (trigger!=ci+1)
    {
  
      for (int k=0; k<spoints[ci]; k++ ) /// read the segment data points (ci=1)
      {
        amplCo[ci][k]=amplC[ci][k];       /// storing original waveform
        hbsl[ci]->Fill(amplCo[ci][k]);    /// event baseline calculation
        sumsig+=(amplCo[ci][k]-meanoffset[ci]);
#ifdef DEBUGMSG 
    diff = (amplCo[ci][k]-meanoffset[ci]) + diff;
	cout<<k<<" = "<<setw(15)<<amplCo[ci][k]<<" --> "<<setprecision(10)<< setw(15)<<(amplCo[ci][k]-meanoffset[ci])  <<" mV  --> sum = "<<sumsig<<"  diff = "<<diff<<endl;
	diff = (amplCo[ci][k]-meanoffset[ci]);
#endif 
          
       }
    } 
#ifdef DEBUG 
#ifdef DEBUGMSG 
      cout<<"Sumsig = "<<sumsig  <<" --> "<<sumsig/502<<"   gain = "<<gain[ci]<<" --> gain/2 = "<<gain[ci]/2.<<endl ;
      cout<< " V/bit = "<< vpb[ci] << " --> V/bit/2 = " <<vpb[ci]/2. <<endl ;
#endif
      cbsl->cd(1+linep*2);
//       cout<<"changing to pad #"<<1+linep*2<<endl;
      TGraph* gr = new TGraph(spoints[ci],tbase,amplC[ci]);
      char graphname[100];
      sprintf(graphname,"channel %d",ci+1);
      gr->SetTitle(graphname);
      gr->Draw("a*l");

      cbsl->cd(2+linep*2);
      gStyle->SetOptStat(1111110);
      hbsl[ci]->Draw();
      cbsl->Modified();
      cbsl->Update();
      linep++;
#endif
      if (trigger==ci+1) continue;

     /// check if histo is empty. In that case either signal completely out of range or problematic event.
/// For this reason we check the sumsig. If it is 0 or very close to 0, reject completely the event!!!!
      
   int hentries= hbsl[ci]->Integral();
   if (hentries == 0)
   {
	 cout<<"_____________________________________________________________"<<endl;
	 cout<<"Event "<<evNo <<" has no point in the range of the baseline histo with "<<hbsl[ci]->GetEntries()<<" entries!"<<endl;
	 cout<<"There are "<<hbsl[ci]->GetBinContent(0)<<" underflow and "<<hbsl[ci]->GetBinContent(nbins+1)<<" points!"<<endl;
	 cout<<"Sumsig for " << spoints[ci] << " is " << sumsig*mV << " mV with a gain of " << vpb[ci]*mV<< " mV/bit"<<endl;
	 if (fabs(sumsig)<vpb[ci]/10. )
	 {
	    cout<<"The event is REJECTED!!!"<<endl;
	    cout<<"-------------------------------------------------------------"<<endl;
	    continue;    
     }
	 else 
	 {
	    int ovfl=hbsl[ci]->GetBinContent(nbins+1);
	    int unfl=hbsl[ci]->GetBinContent(0);
	    int hnp = spoints[ci]/2;
	    if (ovfl>=hnp-1 && unfl>=hnp-1 && fabs(ovfl-unfl)<=1)
	    {	      
	      cout<<"It seems that there might be one point within the range of the scope that makes the frame inegral non-zero."<<endl;
	      cout<<"However, this point is out of the bsl histo range, and the overrange / underange points are shared equally."<<endl; 
	      cout<<"The event is REJECTED!!!"<<endl;
	      cout<<"-------------------------------------------------------------"<<endl;
	      continue;
	    }
	    else
	    {
	      cout<<"The event is registered, thought the bsl histo integral is zero!"<<endl;
	      cout<<"-------------------------------------------------------------"<<endl;
	    }
	  }
    }
    else if (hentries <= 20)
    {	 
	   cout<<"_____________________________________________________________"<<endl;
	   cout<<"Event "<<evNo <<" has very few points point in the range of the baseline histo with "<<hbsl[ci]->GetEntries()<<" entries!"<<endl;
	   cout<<"There are "<<hbsl[ci]->GetBinContent(0)<<" underflow and "<<hbsl[ci]->GetBinContent(nbins+1)<<" points!"<<endl;
	   cout<<"Sumsig for " << spoints[ci] << " is " << sumsig*mV << " mV with a gain of " << vpb[ci]*mV<< " mV/bit"<<endl;
    }
      

      
      if (1 )// || (framesize<10000 && nEv<500) )
      {
        FilterHisto(hbsl[ci],0.03);

        int maxb      = hbsl[ci]->GetMaximumBin();
        double maxAmp =  hbsl[ci]->GetBinContent(maxb);
        double meant = hbsl[ci]->GetBinLowEdge(maxb);
        double rmst  = hbsl[ci]->GetRMS();
        double var = hbsl[ci]->GetBinWidth(maxb);

        if (rmst>0.005) rmst=0.005;
        
        fbsl[ci]->SetParLimits(0,0.5*maxAmp,2.*maxAmp);
        fbsl[ci]->SetParLimits(1,meant-var,meant+2*var);
        fbsl[ci]->SetParLimits(2,vpb[ci]/2.,10.*vpb[ci]);
        fbsl[ci]->SetParameter(0,maxAmp);
        fbsl[ci]->SetParameter(2,rmst);
        fbsl[ci]->SetParameter(1,meant+var/2.);
    // // 	cbsl->cd(1);

        sprintf(fname2,"channel %d event %d baseline",ci+1,evNo );
        hbsl[ci]->SetTitle(fname2);
        fitstatus1[ci] = hbsl[ci]->Fit(fbsl[ci],"BQN","",meant-3*rmst,meant+3*rmst);
        
        bsltmp=fbsl[ci]->GetParameter(1);
        bslC[ci]=fbsl[ci]->GetParameter(1);
        rmstmp=fbsl[ci]->GetParameter(2);
        rmsC[ci]=fbsl[ci]->GetParameter(2);
        
        EVTIME[ci][ngpoints[ci]]=epochF;
        BSL[ci][ngpoints[ci]]=fbsl[ci]->GetParameter(1)*1000.;
        RMS[ci][ngpoints[ci]]=fbsl[ci]->GetParameter(2)*1000.;
        sBSL[ci][ngpoints[ci]]=fbsl[ci]->GetParError(1)*1000.;
        sRMS[ci][ngpoints[ci]]=fbsl[ci]->GetParError(2)*1000.;
        ngpoints[ci]++;

    #ifdef DEBUG
        fbsl[ci]->Draw("same");
        cbsl->Modified();
        cbsl->Update();
    #endif
        hbslall[ci]->Add(hbsl[ci]);

      }

      maxpoints = inpoints;
      
      if (maxpoints <= 0) continue;


      int tst = SubtractBaseline(maxpoints,amplC[ci],amplCo[ci],bslC[ci]);//amplC corresponds to initial signal and is amplCo is saved in the output tree

// //       hbsl0->Reset("");  /// this is the subtracted event baseline !!!!
// // 
// //       for (int k=0; k<maxpoints; k++ ) /// read the segment data points
// //       {	     
// // 	if (fabs(amplSum[k])<0.005)
// // 	hbsl0->Fill(amplSum[k]);
// //       }
// // 	  
// //       /// calculate baseline
// //       int maxbt = hbsl0->GetMaximumBin();
// //       double baxby = hbsl0->GetBinContent(maxbt);
// //       double meantt = hbsl0->GetBinLowEdge(maxbt);
// //       double vart = hbsl0->GetBinWidth(maxbt);
// //       double rmstt = hbsl0->GetRMS();
// //       // 	  cout<<"meant C"<<ci+1<<" = "<<meant<<endl;
// // 
// //       if (rmstt<0.005) /// don't accumulate wierd events in total baseline histo;
// //       hbsl0all->Add(hbsl0);

/// fit and draw the basline-subtructed signal new baseline (must be ~0) 
    
// //       fbsl0->SetParameter(0,baxby+vart/2.);
// //       fbsl0->SetParameter(2,rmstt);
// // //       fbsl[ci]->SetParameter(2,0.002);
// //       fbsl0->SetParameter(1,meantt);
// //       
// //       fbsl0->SetParLimits(0,0.5*baxby,2.*baxby);
// //       fbsl0->SetParLimits(1,meantt-vart,meantt+2*vart);
// //       fbsl0->SetParLimits(2,vpb/2.,10.*vpb);
// // 
// //       cbsl->cd(2);
// //       sprintf(fname2,"bsl-corrected event %d baseline",evNo );
// //       hbsl0->SetTitle(fname2);
// //       fitstatus2 = hbsl0->Fit(fbsl0,"BQ","",meantt-0.0075,meantt+0.0075);
// //       bslsum=fbsl0->GetParameter(1);
// //       rmssum=fbsl0->GetParameter(2);
// //       hbsl0->Draw("E0");
      
    if (fitstatus1[ci] != 0 )//|| fitstatus2 || 0)
    {
        cout<<"evNo = "<<evNo<<"\t  f1 = "<<fitstatus1[ci]<<" channel = " << ci+1;
        cout<<" bsl = " <<bsltmp<<" rms = " << rmstmp;
// //         cout<<" bsl0 = " <<bslsum<<" rms0 = " << rmssum<<endl;
        cout<<" SumSig = " <<sumsig<<" histogram Integral = " << hentries<<endl;
// //         cbsl->Modified();
// //         cbsl->Update();
// 	return (fitstatus1);
    }
	
//       hbsl0->GetXaxis()->SetRangeUser(floor(bslsum*100.+0.5)/100-0.02,floor(bslsum*100.+0.5)/100+0.02);
#ifdef DEBUG 
    cbsl->Modified();
    cbsl->Update();
//     continue;
//     return -4;  
#endif
/// fit and draw the acumulated baseline of the signals
    if ((nEv==eventF-1) || (framesize<10000 && nEv%50==0))
    {
// //         cbsl->cd(3);
        double maxbt0 = hbslall[ci]->GetMaximumBin();
        double maxv = hbslall[ci]->GetBinContent(maxbt0);
        double meantt0 = hbslall[ci]->GetBinLowEdge(maxbt0);
        double rmstt0 = hbslall[ci]->GetRMS();
	
        fbslall[ci]->SetParameter(0,maxv);
        fbslall[ci]->SetParameter(1,meantt0);
        fbslall[ci]->SetParameter(2,rmstt0);
        fbslall[ci]->SetParLimits(0,0.5*maxv,1.5*maxv);
        fbslall[ci]->SetParLimits(1,bsltmp-0.01,bsltmp+0.01);
        fbslall[ci]->SetParLimits(2,0.0005,0.03);

        hbslall[ci]->Fit(fbslall[ci],"BQN","",bsltmp-0.01,bsltmp+0.01);
// // 	hbslall->Draw("E0");
        bsltmp=fbslall[ci]->GetParameter(1);
// // 	cbsl->Modified();
// // 	cbsl->Update();
      }

    }
/// end of the 4 channel loop      
     // if (trigExist)
        //cout<<"Found srs event ID = "<<srsNo<<endlr;
      //else
        //cout<<RED<<"There is no SRS trigger for this run"<<endlr;
    
    trackOK=0;
    

    if (trackingFileOK>0)
    {
        if (srsNo>maxTrackID+1 || srsNo<0)
            eventTracks=0;
        else
            eventTracks = trackcollection[srsNo].ntracks;
        
        trackOK=(eventTracks>0?1:0);
    }
    else eventTracks=0;

    for (int ntr=0;ntr<eventTracks;ntr++)
    {
        chi2track[ntr]=trackcollection[srsNo].trackchi2[ntr];
    }
    
///   Now assign the tracking data for the channel ci
    for (int ci=0;ci<4;ci++)
    {
      int cpos=cposition[ci];
      
      if (trackOK==1 && (ci+1 != trigger)&&trigExist)
      {
        for(int j=0;j<eventTracks;j++)
        {
            for(int l=0;l<3;l++)
            {
                refP1[j][l]= trackcollection[srsNo].hits[j][p1][l];
                refP2[j][l]= trackcollection[srsNo].hits[j][p2][l];
            }
            if (cpos>=0 && cpos<20)
            {
                XX[ci][j]=trackcollection[srsNo].hits[j][cpos][0];
                YY[ci][j]=trackcollection[srsNo].hits[j][cpos][1];
            }
            else 
            {
                double z = oscsetup->Z[ci];
                double t = (z-refP1[j][2]) / (refP2[j][2]-refP1[j][2]);
                XX[ci][j] = refP1[j][0] + t * (refP2[j][0] - refP1[j][0]);
                YY[ci][j] = refP1[j][1] + t * (refP2[j][1] - refP1[j][1]); 
            }

            for (int l=0;l<6;l++)
            {
                totchargenextcluster[j][l]=trackcollection[srsNo].totchanextcluster[j][l];
                disttonextcluster[j][l]=trackcollection[srsNo].distnextcluster[j][l];
            }
            
        }
      }
    }
#ifdef DEBUG 
    linep=0;
    for (int i=0;linep<actch && i<4;i++)
    {
        if (!active[i]) continue;
        cbsl->cd(2+linep*2);
//       cout<<"changing to pad #"<<1+linep*2<<endl;
        char txt2[1000];
        cout<<"Event "<<nEv <<" SRS "<<srsNo<<" has "<<eventTracks<<" tracks"<<endl;
        infotree->GetEntry(i);
        for (int j=0;j<eventTracks;j++)
        {
//             double mmax = hbsl[i]->GetMaximum();
            sprintf(txt2," #it{hit #%d} = #bf{(%3.1f,%3.1f,%3.1f)} ",j, XX[i][j],YY[i][j],ZZ); 
            TLatex* tex = new TLatex(bslmin[i]+range[i]/8.*0.2,4500-j*450.,txt2);
            tex->SetTextFont(132);
            tex->SetTextSize(0.045);
            tex->SetLineWidth(2);
            tex->SetTextAlign(12);
            if (trigger != i+1)tex->Draw();
        }
      linep++;
    }
    cbsl->Modified();
    cbsl->Update();
    continue;

#endif


// for (int i=0 ;i <77; i++ )
// {
//     tracktree->GetEntry(i);
//     int etmp=srstriggerctr;
//     cout<<"searching for evNo "<<srsNo <<" , found "<<etmp <<" at position "<<i<<"  "<< srstimestamp<<endl;
// }
// return 88;

/// Find the corresponding tracking event;

//     if (trackingFileOK>0)
//     {
//         trackOK=0;
//         int position=tindex;
//         while (tindex<tracktree->GetEntries())
//         {
//            tracktree->GetEntry(tindex);
//            int etmp=srstriggerctr;
//             cout<<"searching for evNo "<<srsNo <<" , found "<<etmp <<" at position "<<tindex<< endl;
//            if (etmp==srsNo) 
//            {
//                position=tindex;
//                tindex++;
//                trackOK=1;
//                break;
//            }
//            else if (etmp>srsNo)
//            {
//                tindex--;
// //               position = -1;
//                 cout<<MAGENTA<<"Event = "<<evNo <<" ,  srsNo = "<< srsNo <<" has no srstriggerctr "<<endlr;
//                trackOK=0;
//                break;
//            }
//            tindex++;
//         }
//         for (int i=1;i<2;i++)
//         {
//             if (trigger!=i+1 && trackOK==1)
//             {
//                 int nTracks = ntracks;
//                 int dsize = distnextcluster->size();
//                 int tsize = totchanextcluster->size();
// //                 if (totchanextcluster!=NULL) tsize=totchanextcluster->size();
// //                 xy[i][0] = distnextcluster->at(3);
// //                 xy[i][1] = distnextcluster->at(2);
// //               if (nTracks<2 ) continue;
//                 
//                 cout<<"srstriggerctr "<<srstriggerctr <<" srstimestamp "<<srstimestamp<<" ntracks "<< ntracks<<" tracknumber "<< tracknumber<<endl;
//                 cout<<"Event = "<<evNo <<" ,  srsNo = "<< srsNo <<" found at tindex = "<< position<<" of tracking tree. trackOK is "<<trackOK<<endl;
//                  cout<<distnextcluster->size() <<" elements, ";
//                  for (int j=0; j<dsize; j++) cout<<j<<" = "<<distnextcluster->at(j) <<" ";
//                  cout<<endl<<totchanextcluster->size() <<" elements, ";
//                  for (int j=0; j<tsize; j++) cout<<j<<" = "<<totchanextcluster->at(j) <<" ";             
//                  cout<<endl<<trackchi2<<"  ________________"<<endl;
//                 int counter=0;
//                 for(vector<vector<double> >::iterator it = hits->begin(); it != hits->end(); ++it)
//                 {
// //                     cout<<"Z of position "<<counter<<" = "<<hits->at(counter).at(2)<<endl;
//                     counter++;
//                 }
//                 
//             }
//             else
//             {
// //                 XX[i][0] = -1111;
// //                 xy[i][1] = -1111;
//             }
//         }
//             
// //         cout<<"Event = "<<evNo <<" ,  srsNo = "<< srsNo <<" found at tindex = "<< position<<" of tracking tree. trackOK is "<<trackOK<<endl;
//     }
//       
      
     outtree->Fill();

     if (nEv==20) 
     {
         outtree->OptimizeBaskets(100000000,1.1,"d" );
     }
      
     if (framesize>50000)
     {
        if ((nEv+1) % 10 ==0)
            cout<<"Processed "<<nEv+1<<" events (evNo = "<<evNo<<") out of "<<nevents<<endl;
     }
     else
     {
        if ((nEv) % 1000 ==0)
        {
            cout<<"Processed "<<nEv+1<<" events (evNo = "<<evNo<<") out of "<<nevents<<endl;
            TTimeStamp *tstamp = new TTimeStamp();
            tstamp->Set();
            cout<<"\t\tprocessing time : "<<tstamp->AsString("l")<<endl;
            cout<<"SRS No = "<<srsNo<<" for event " <<evNo<<endl;    
        }
     }

     if (nevents+1 % 10000 == 0)
	  ofile->Write("",TObject::kOverwrite);
      
  }   //// end of the tree events loop  
  
  
  
  
  TCanvas *ebsl[4];
  
  for (int ci=0;ci<4;ci++)
  {
    if (!active[ci] || trigger==ci+1) continue;
   
    TString channel = TString::Itoa(ci+1,10);
    TString cname(rtypes);
    cname+="-C";
    cname+=channel;
    cname+="-";
    cname+=oscsetup->Amplifier[ci];
    ebsl[ci] = new TCanvas("baselineEvolution"+cname,"Baseline Evolution "+cname);
    
    TGraphErrors* grBSL = new TGraphErrors(ngpoints[ci],EVTIME[ci],BSL[ci],0,sBSL[ci]);
    grBSL->SetTitle("BASELINE " + cname);
    grBSL->SetLineColor(4);
    grBSL->Draw("AP");
    TH1F *h1 = grBSL->GetHistogram();
    h1->SetBins(100000,EVTIME[ci][0]-0.05*(EVTIME[ci][ngpoints[ci]-1]-EVTIME[ci][0]),EVTIME[ci][ngpoints[ci]-1]+0.05*(EVTIME[ci][ngpoints[ci]-1]-EVTIME[ci][0]));
    grBSL->GetXaxis()->SetTimeDisplay(1);
    grBSL->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%d/%m/%y}");
    grBSL->GetXaxis()->SetLabelOffset(0.04);
    grBSL->GetYaxis()->SetTitle("[mV]");
    //   grBSL->SetMinimum(0);
    grBSL->Write("baselineEvolution"+channel);
  

    TCanvas *erms = new TCanvas("baselineRMSEvolution"+cname,"Baseline RMS Evolution "+cname);
    TGraphErrors* grRMS = new TGraphErrors(ngpoints[ci],EVTIME[ci],RMS[ci],0,sRMS[ci]);
    grRMS->SetTitle("RMS " + cname);
    grRMS->SetLineColor(2);
    grRMS->Draw("AP");
    TH1F *h2 = grRMS->GetHistogram();
    h2->SetBins(100000,EVTIME[ci][0]-0.05*(EVTIME[ci][ngpoints[ci]-1]-EVTIME[ci][0]),EVTIME[ci][ngpoints[ci]-1]+0.05*(EVTIME[ci][ngpoints[ci]-1]-EVTIME[ci][0]));
    grRMS->GetXaxis()->SetTimeDisplay(1);
    grRMS->GetXaxis()->SetLabelOffset(0.04);
    grRMS->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%d/%m/%y}");
    grRMS->GetYaxis()->SetTitle("[mV]");
  //   grRMS->SetMinimum(grRMS->GetMinimum()*0.7);
  //   grRMS->SetMaximum(grRMS->GetMaximum()*1.3);
    grRMS->Write("rslRMSEvolution"+channel);
    
    TCanvas *ebslrms = new TCanvas("BSL+RMS_Evolution"+cname,"Baseline & RMS Evolution "+cname);
    TGraphErrors* grBSLRMS = new TGraphErrors(ngpoints[ci],EVTIME[ci],BSL[ci],0,RMS[ci]);
    grBSLRMS->SetTitle("BSLRMS " + cname);
    grBSLRMS->SetLineColor(4);
    grBSLRMS->Draw("AP");
    TH1F *h3 = grBSLRMS->GetHistogram();
    h3->SetBins(100000,EVTIME[ci][0]-0.05*(EVTIME[ci][ngpoints[ci]-1]-EVTIME[ci][0]),EVTIME[ci][ngpoints[ci]-1]+0.05*(EVTIME[ci][ngpoints[ci]-1]-EVTIME[ci][0]));
  //   TH1D *h1 = grBSL->GetHistogram();
    grBSLRMS->GetXaxis()->SetTimeDisplay(1);
    grBSLRMS->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%d/%m/%y}");
    grBSLRMS->GetXaxis()->SetLabelOffset(0.04);
    grBSLRMS->GetYaxis()->SetTitle("[mV]");
  //   grBSLRMS->SetMinimum(0);
    grBSLRMS->Write("BSLandRMSevolution"+channel);

  //   gStyle->SetOptFit(110);
    TCanvas *ctotbsl = new TCanvas("Pedestals"+cname,"Pedestals "+cname,1600*4/4,1200*3/4);
    ctotbsl->Divide(2,2);
  //   ctotbsl->UseCurrentStyle();
    ctotbsl->cd(1);
    hbslall[ci]->Draw();
    ctotbsl->cd(2);
    grBSL->Draw("AP");
    ctotbsl->cd(3);
    grRMS->Draw("AP");
    ctotbsl->cd(4);
    grBSLRMS->Draw("AP");
    ctotbsl->Modified();
    ctotbsl->Update();
    ctotbsl->Write();

    char tmpdir[500];
    sprintf(tmpdir,"%s",gSystem->pwd());
    cout<<"Actual directory: "<<tmpdir<<endl; 
    char plotdirname[500];
    sprintf(plotdirname,"%s/Run%03d-%s/",PLOTDIR,runNo,poolid);//,abs(threshold));
    
    gSystem->mkdir(plotdirname,kTRUE);
    gSystem->ChangeDirectory(plotdirname); 
    ctotbsl->SaveAs(".pdf");
    ctotbsl->SaveAs(".png");
//     gSystem->ChangeDirectory("./../..");
    gSystem->ChangeDirectory(tmpdir);
  }
  
  cout<<"End of file processing.\n"<<outtree->GetEntries()<<" events were found"<<endl;
  ofile->Write("",TObject::kOverwrite);
  

  return 0;
}

//#endif
