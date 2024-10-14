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
// #include "MyFunctions.h"
#include "MyFunctions.C"


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

int FindZeros(int n, double* data)
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

int month2int_littleendian(const char *month)
{
  const char hash[] = { 3, 12, 8, 2, 1, 11, 7, 5, 0, 10, 4, 9, 6 };
  return hash[(* (const int32_t *) month & ~0x20202020 + 146732) % 13];
}


int main(int runNo=1, int poolNo=1, int trigger=0, string filetype = "")  ///put trigger = ChannelNo in case of external trigger
{
  gROOT->LoadMacro("MyFunctions.C");
  
  if (trigger>4){
    cout<<"Wrong trigger channel : "<<trigger<<endl;
    return -15;
  }
  
  const double mV =1000.;
  
  int detNo=0;
  
  
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
//    cout<<command<<endl<<"returned: "<<tst<<endl;
    if (tst !=0)
    {
      cout<<command<<endl<<"returned: "<<tst<<endl;
      cout<<"Probably the tree of the run "<<runNo<<" in Pool "<<poolNo<<" was not found in directory "<<datadirname<<endl<<"Exiting..."<<endl;
      return tst;
    }
    ftmp=fopen(afile,"r");
    if (fgets(fnametmp,200,ftmp) == NULL)
    {
      cout<<"Failed to read the filename for Pool "<<poolNo<<" run "<<runNo<<" at "<<datadirname<<endl<<"Exiting..."<<endl;
      return -2;
    }
    int rtmp = runNo, dtmp = 1;
    float Rd;
    char ftypetmp[500];
    strcpy(ftypetmp,"");
//     int stst = sscanf(fnametmp,"%3d-%02d-%d-%d_%30[^ /,\n\t]",&rtmp,&dtmp,&vm,&vd,ftypetmp);
    int stst=0;
    cout <<"___________________________________________________\n\narguments read = "<<stst<<endl;
    filetype.assign(ftypetmp);
    Rd=0;
    cout<<"run = "<<rtmp<<" ("<<runNo<<")"<<endl;
    cout<<"amplifier = "<<dtmp<<endl;
    detNo=dtmp;
    cout<<"Vm = "<<vm<<endl;
    cout<<"Vd = "<<vd<<endl;
    cout<<"distance = "<<Rd<<endl;
    dgap=Rd;
//       cout<<"Drift gap = "<<dgap<<endl;
    if (stst>5) 
        cout<<"filetype = "<<filetype<<endl;
    fclose(ftmp);
    sprintf(command,"rm %s\n",afile);
    tst=system(command);
  }
  else
  {
    cout<<"Available runs: "<<MINRUN<<" - "<<MAXRUN<<endl;
    return (-2);
  }

  char amplifier[20];
  if (detNo==1)
    sprintf(amplifier,"cividec");
  else if (detNo==2)
    sprintf(amplifier,"ortec");
  else if (detNo==3)
    sprintf(amplifier,"ATHL");
  else if (detNo==4)
    sprintf(amplifier,"CERN");
  else 
    sprintf(amplifier,"other");
  
  
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
  cout<<"output filename = "<<ofname<<endl;

  
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
    cout<<"Failed to open \""<<ifilename<<"\"\ncorresponding to the run "<<detNo<<"\nExiting..."<<endl;
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
    bname = "amplC"+channel;
    cout << bname << endl;
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

  
  const int nadiv = 2;
  
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
      cout<<"\n\nFramesize = "<<framesize<<" points"<<endl;
      cout<<"Gain = "<<gain[i]<<" V"<<endl;
      cout<<"Offset = "<<offset[i]<<" V"<<endl;
      cout<<"Range Max = "<< rmax[i]*mV<<" mV";
      cout<<"\t\tRange Min = "<< rmin[i]*mV<<"\t";
      range[i]=rmax[i]-rmin[i];
      cout<<"Full range = "<<range[i]*mV<<" mV  -->  " << range[i]/8 *mV <<" mV / div  --> "<<range[i]/256*mV<<" mV / bit"<<endl<<endl;
      bslmax[i] = rmax[i];
      bslmin[i] = rmin[i] + range[i]*((8.-nadiv)/8.);
//       cout<<"Corresponding to "<<endl<<endl<<endl;
//       break;  /// only the first active channel is used!!!!!!
       meanoffset[i]= (rmax[i] + rmin[i])/2. ;
      cout<<"Mean offset C"<<i+1<<" = "<< meanoffset[i]*mV <<" mV"<<endl<<endl;
      vpb[i] = range[i]/256.;
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

  TString rtypes(runid);
  rtypes+="-";
  rtypes+=amplifier;
//   rtypes+="-";
//   rtypes+=mesh;
//   rtypes+="-";
//   rtypes+=drift;
//   rtypes+="-";
//   rtypes+=ctrigger;
  rtypes+="-";
  rtypes+=rtype;
  
  char catt[30];
  
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
  sprintf(ofilename,"%s/%s",outdirname,ofname);
  TFile *ofile = new TFile(ofilename,"RECREATE");

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
    if (trigger==ci+1) continue;
    sprintf(fname2,"Signal_bsl_%d",ci+1);
    hbsl[ci] = new TH1F(fname2,fname2,nbins,bslmin[ci],bslmax[ci]);
    hbsl[ci]->GetXaxis()->SetTitle("[V]");
    cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;
    sprintf(fname2,"Signal_bsl_tot_%d",ci+1);
    sprintf(htitle,"Signal baseline from all events");
    hbslall[ci] = new TH1D(fname2,htitle,nbins,bslmin[ci],bslmax[ci]);
    hbslall[ci]->GetXaxis()->SetTitle("[V]");
    cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;

    sprintf(fname2,"Subtracted_bsl_%d",ci+1);
    hbsl0[ci] = new TH1F(fname2,fname2,(int)nbins*binscale,-range[ci]/8.,range[ci]/8.);
    hbsl0[ci]->GetXaxis()->SetTitle("[V]");
    cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;
    sprintf(fname2,"Subtracted_bsl_tot_%d",ci+1);
    hbsl0all[ci] = new TH1D(fname2,fname2,(int)nbins*binscale,-range[ci]/8.,range[ci]/8.);
    hbsl0all[ci]->GetXaxis()->SetTitle("[V]");
    cout<<"defined histogram["<<ci<<"] for channel C"<<ci+1<<" : "<< fname2 <<endl;
  }

  double *amplCo[4];
  for (int i=0; i<4; i++)
  {
    if (trigger!=i+1 && active[i])
      amplCo[i] = new double[FRAMESIZE]; 
  }

//   double *amplSum; //in this case (lina4) is just ampl - baseline 
//     amplSum = new double[FRAMESIZE];
  
  double *tbase = new double[FRAMESIZE];
  for (int i=0;i<FRAMESIZE;i++)
    tbase[i]=i*dt;
    
  double *EVTIME[4], *BSL[4], *RMS[4], *sBSL[4], *sRMS[4];
  
  
  for (int i=0; i<4; i++)
  {
      if (trigger!=i+1 && active[i])
      {
	EVTIME[i] = new double[nevents];  
	BSL[i] = new double[nevents];  
	RMS[i] = new double[nevents];  
	sBSL[i] = new double[nevents];  
	sRMS[i] = new double[nevents];  
      }
  }
  
  
    
  int maxpoints=0;
  int srsNo=0;

  TTree *outtree;
  
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

  outtree->Branch("sumpoints", &maxpoints, "sumpoints/I");
  for (int i=0; i<4; i++)
  {
    if (trigger!=i+1 && active[i])
    {
       TString channel = TString::Itoa(i+1,10);
       TString bname = "amplC"+channel;
       TString btype = bname+"[sumpoints]/D";
       outtree->Branch(bname, &amplCo[i][0], btype);
    }
  }
//   outtree->Branch("amplSum2", amplSum[1], "amplSum2[sumpoints]/D");
  double bslsum[]={0,0,0,0};
  double rmssum[]={0.1,0.1,0.1,0.1};
  outtree->Branch("bslSum", bslsum, "bslSum[4]/D");
  outtree->Branch("rmsSum", rmssum, "rmsSum[4]/D");
  
  cout<<"\n\n Tree is ready!"<<endl;

#ifdef DEBUG 
//    TCanvas *cbsl = new TCanvas("baseline"+rtypes,"baseline "+rtypes,1600*4/4,1200*3/4);
    if (trigger>0) actch-=1;
    TCanvas *cbsl = new TCanvas("baseline"+rtypes,"baseline "+rtypes,1600*4/4,1200*actch*2/4);
    cout<<" Atchu = "<<actch<<endl;
    cbsl->Divide(2,actch);
    int linep = 0;
#endif  

  double tshift[2]={0.,0.};
//   int dd,mm,yy,hour,min,sec;
//   char month[20];
  int nseg;

  double bsl[2]={0.,0.};
  double rms[2]={0.,0.};
  
  
  sprintf(fname2,"gausC");
  TF1* fbsl = new TF1(fname2,"gaus",-1.,1.);
  fbsl->SetNpx(2000);
  sprintf(fname2,"gaus0");
  TF1* fbsl0 = new TF1(fname2,"gaus",-1.,1.);
  fbsl0->SetNpx(2000);

  sprintf(fname2,"gausCtot");
  TF1* fbslall = new TF1(fname2,"gaus",-1.,1.);
  fbslall->SetNpx(2000);
  sprintf(fname2,"gaus0tot");
  TF1* fbsl0all = new TF1(fname2,"gaus",-1.,1.);
  fbsl0all->SetNpx(2000);

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
//   cout<<"Loading tree baskets in virtual memory"<<endl;
//   int nbaskets=intree->LoadBaskets(4000000000);
//   cout<<"Loaded "<<nbaskets<<" baskets "<<endl;
// #endif
  
  int oveeflow=0;
  int lastevent=-1;
  double tdt = -1.;  
  for (int nEv=eventS; nEv<nevents && nEv<eventF+1; nEv+=1)
  { 
#ifdef DEBUG 
//     eventS = 16438;
    cout<<"Give event No to debug: ";
    cin>>nEv;
    if (nEv<0) return (-1);
    linep = 0;
//     eventF = eventS;
#endif

    intree->GetEntry(nEv);

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
  
    for (int ci=0; ci<4;ci++)
    {
//       cout<<"Searching for trigger channel, c"<<ci+1<<endl;
      if (!active[ci]) continue;
      
      itrigger = -100;

      if (trigger == ci+1)   /// calculate the SRS event number here
      {
        // 	srsNo = CalculateSRSEventNo(amplC[ci],inpoints,overflow,lastevent);
//         cout << srsNo<<endl;
        lastevent=srsNo;

        for (int k=0; k<spoints[ci]-2; k++ )  ///processing trigger (ci=0)
        {
            if(amplC[ci][k]<-0.3) 
            {
                itrigger =k;
                break;
            }
        }
        ttrig = itrigger * dt;
//         cout<<"itrigger = "<<itrigger<<" , ttrig = "<<ttrig<<" ns"<<endl;
        continue;
      }
  
      
      frangemax[ci] = rmax[ci];
      frangemin[ci] = rmin[ci];


      hbsl[ci]->Reset("");

//       return 3;

      long double sumsig=0.;
      long double diff=0.;
      for (int k=0; k<spoints[ci]; k++ ) /// read the segment data points (ci=1)
      {
        amplCo[ci][k]=amplC[ci][k];       /// storing original waveform
        hbsl[ci]->Fill(amplCo[ci][k]);    /// event baselin calculation
        sumsig+=(amplCo[ci][k]-meanoffset[ci]);
#ifdef DEBUGMSG 
    diff = (amplCo[ci][k]-meanoffset[ci]) + diff;
	cout<<k<<" = "<<setw(15)<<amplCo[ci][k]<<" --> "<<setprecision(10)<< setw(15)<<(amplCo[ci][k]-meanoffset[ci])  <<" mV  --> sum = "<<sumsig<<"  diff = "<<diff<<endl;
	diff = (amplCo[ci][k]-meanoffset[ci]);
#endif 
          
       }
      
#ifdef DEBUG 
#ifdef DEBUGMSG 
      cout<<"Sumsig = "<<sumsig  <<" --> "<<sumsig/502<<"   gain = "<<gain[ci]<<" --> gain/2 = "<<gain[ci]/2.<<endl ;
      cout<< " V/bit = "<< vpb[ci] << " --> V/bit/2 = " <<vpb[ci]/2. <<endl ;
#endif
      cbsl->cd(1+linep*2);
      cout<<"changing to pad #"<<1+linep*2<<endl;
      TGraph* gr = new TGraph(spoints[ci],tbase,amplCo[ci]);
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
	
	fbsl->SetParLimits(0,0.5*maxAmp,2.*maxAmp);
	fbsl->SetParLimits(1,meant-var,meant+2*var);
	fbsl->SetParLimits(2,vpb[ci]/2.,10.*vpb[ci]);
	fbsl->SetParameter(0,maxAmp);
	fbsl->SetParameter(2,rmst);
	fbsl->SetParameter(1,meant+var/2.);
// // 	cbsl->cd(1);

	sprintf(fname2,"channel %d event %d baseline",ci+1,evNo );
	hbsl[ci]->SetTitle(fname2);
	fitstatus1[ci] = hbsl[ci]->Fit(fbsl,"BQN","",meant-3*rmst,meant+3*rmst);
      
	bsl[1]=fbsl->GetParameter(1);
	rms[1]=fbsl->GetParameter(2);
	
	EVTIME[ci][ngpoints[ci]]=epochF;
	BSL[ci][ngpoints[ci]]=fbsl->GetParameter(1)*1000.;
	RMS[ci][ngpoints[ci]]=fbsl->GetParameter(2)*1000.;
	sBSL[ci][ngpoints[ci]]=fbsl->GetParError(1)*1000.;
	sRMS[ci][ngpoints[ci]]=fbsl->GetParError(2)*1000.;
	ngpoints[ci]++;

#ifdef DEBUG
	fbsl->Draw("same");
    cbsl->Modified();
	cbsl->Update();
#endif
	hbslall[ci]->Add(hbsl[ci]);

	bsl[0]=0.0;
	rms[0]=0.0;

      }

      maxpoints = inpoints;
      
      if (maxpoints <= 0) continue;


      int tst = SubtractBaseline(maxpoints,amplC[ci],amplCo[ci],bsl[1]);//amplC corresponds to initial signal and is amplCo is saved in the output tree

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
	cout<<" bsl = " <<bsl[1]<<" rms = " << rms[1];
	cout<<" bsl0 = " <<bslsum<<" rms0 = " << rmssum<<endl;
	cout<<" SumSig = " <<sumsig<<" histogram Integral = " << hentries<<endl;
// //         cbsl->Modified();
// //         cbsl->Update();
// 	return (fitstatus1);
      }
	
//       hbsl0->GetXaxis()->SetRangeUser(floor(bslsum*100.+0.5)/100-0.02,floor(bslsum*100.+0.5)/100+0.02);
#ifdef DEBUG 
      cbsl->Modified();
      cbsl->Update();
      continue;
      return -4;  
#endif
/// fit and draw the acumulated baseline of the signals
      if ((nEv==eventF-1) || (framesize<10000 && nEv%50==0))
      {
// //         cbsl->cd(3);
	double maxbt0 = hbslall[ci]->GetMaximumBin();
 	double maxv = hbslall[ci]->GetBinContent(maxbt0);
        double meantt0 = hbslall[ci]->GetBinLowEdge(maxbt0);
        double rmstt0 = hbslall[ci]->GetRMS();
	
	fbslall->SetParameter(0,maxv);
	fbslall->SetParameter(1,meantt0);
	fbslall->SetParameter(2,rmstt0);
	fbslall->SetParLimits(0,0.5*maxv,1.5*maxv);
	fbslall->SetParLimits(1,bsl[1]-0.01,bsl[1]+0.01);
	fbslall->SetParLimits(2,0.0005,0.03);

	hbslall[ci]->Fit(fbslall,"BQN","",bsl[1]-0.01,bsl[1]+0.01);
// // 	hbslall->Draw("E0");
	bsl[1]=fbslall->GetParameter(1);
// // 	cbsl->Modified();
// // 	cbsl->Update();
      }

    }
/// end of the 4 channel loop      
      
      
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
  
  cout<<"End of file processing.\n"<<evNo<<" events were found"<<endl;
  ofile->Write("",TObject::kOverwrite);
  

  return 0;
}

//#endif
