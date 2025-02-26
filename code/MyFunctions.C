#ifndef MYFUNCTIONS_C
#define MYFUNCTIONS_C 1

#include "MyFunctions.h"

ClassImp(PEAKPARAM);

void MyFunctions(){
  cout<<"Loading functions library..."<<endl;}


int GetOsciloscopeSetup(int runNo, int poolNo, const char* dirname, OSCSETUP* oscsetup)
{

  char filename[1000];
  sprintf(filename,"%s/OsciloscopeSetup.txt",dirname);
  ifstream inf;
  inf.open(filename); 
  int entries=0;
  if (!inf.is_open()) /// running for the first time!
  {
     cout<<RED<<"\n\n\nParameter file "<<filename<<" doesn't exist. Exiting ...\n\n\n"<<endl;
     return -1;
  }
  else /// read existing table entries 
  {
    inf.ignore(1000,'\n');  /// skip the headerline
    while (!inf.eof() )
    {
       int runtmp=-1;
       int pooltmp=-1;
       inf>>runtmp;
       if (inf.eof()) break;
       inf>>pooltmp;
       oscsetup->poolNo=pooltmp;
       oscsetup->runNo=runtmp;
       inf>>oscsetup->srsCh;
       for (int i=0;i<4;i++)
       {
           inf>>oscsetup->DetName[i];
           inf>>oscsetup->V1[i];
           inf>>oscsetup->V2[i];
           inf>>oscsetup->Z[i];
           inf>>oscsetup->Photocathode[i];
           inf>>oscsetup->Amplifier[i];
//            cout <<BLUE<<"Read ["<<i<<"]:"<<endl;
//            cout<<oscsetup->DetName[i]<<"\t";
//            cout<<oscsetup->V1[i]<<"\t";
//            cout<<oscsetup->V2[i]<<"\t";
//            cout<<oscsetup->Photocathode[i]<<"\t";
            cout<<MAGENTA<<"Amplifier ->"<<oscsetup->Amplifier[i]<<"<-"<<endlr;
           oscsetup->AmplifierNo[i]=WhichAmplifier(oscsetup->Amplifier[i]);
       }
       entries++;
       if (runtmp==runNo && pooltmp==poolNo)
       {
           inf.close();
           return entries; 
       }
       
       if (entries>100000) return -11;
    }
    
    cout<<RED<<endl<<"Run "<<runNo<<" Pool"<<poolNo<<" NOT found in the setup file "<<filename<<" which contains "<<entries <<" entries!"<<endl<<endlr;
    inf.close();
    return 0;    
  }
  
}
    
void AssignRun(RUNPAR* allruns,RUNPAR *runpar,int irun)
{
  allruns[irun].detNo=runpar->detNo;
  allruns[irun].runNo=runpar->runNo;
  allruns[irun].ampl=runpar->ampl;
  allruns[irun].sampl=runpar->sampl;
  allruns[irun].charge=runpar->charge;
  allruns[irun].scharge=runpar->scharge;
  allruns[irun].rate=runpar->rate;
  allruns[irun].srate=runpar->srate;
  allruns[irun].grate=runpar->grate;
  allruns[irun].sgrate=runpar->sgrate;
  allruns[irun].tot=runpar->tot;
  allruns[irun].stot=runpar->stot;
  allruns[irun].risetime=runpar->risetime;
  allruns[irun].srisetime=runpar->srisetime;
  allruns[irun].width=runpar->width;
  allruns[irun].swidth=runpar->swidth;
  allruns[irun].chovampl=runpar->chovampl;
  allruns[irun].schovampl=runpar->schovampl;
   
}
void RegisterRunParameters(RUNPAR* runpar, const char* dirname)
{
  const char title[18][12]={"detNo","runNo","ampl","sampl","charge","scharge","rate","srate","grate","sgrate","tot","stot","risetime","srisetime","width","swidth","chovampl","schovampl"};
  RUNPAR *allruns = new RUNPAR[TOTRUNS*100];
  for(int i=0;i<TOTRUNS*100;i++) {
    allruns[i].runNo=-1;
    allruns[i].detNo=-1;
  }
  char filename[1000];
  sprintf(filename,"%s/ParameterTable.txt",dirname);
  ifstream inf;
  inf.open(filename); 
  int entries=0;
  if (!inf.is_open()) /// running for the first time!
  {
     cout<<"\n\n\nParameter file "<<filename<<" doesn't exist. Creating new ...\n\n\n"<<endl;
     FILE* outf = fopen(filename,"w");
     for (int i=0;i<18;i++) fprintf(outf,"%12s ",title[i]);
     fprintf(outf,"\n");
     fclose(outf);
  }
  else /// read existing table entries 
  {
    inf.ignore(1000,'\n');  /// skip the headerline
    while (!inf.eof() )
    {
       int dettmp=-1;
       int runtmp=-1;
       inf>>dettmp;
       if (inf.eof()) break;
       inf>>runtmp;
//        int irun = runtmp-MINRUN;
       int irun = 100*dettmp + runtmp;
       
       if (irun<0) break;
       allruns[irun].detNo=dettmp;
       allruns[irun].runNo=runtmp;
       inf>>allruns[irun].ampl;
       inf>>allruns[irun].sampl;
       inf>>allruns[irun].charge;
       inf>>allruns[irun].scharge;
       inf>>allruns[irun].rate;
       inf>>allruns[irun].srate;
       inf>>allruns[irun].grate;
       inf>>allruns[irun].sgrate;
       inf>>allruns[irun].tot;
       inf>>allruns[irun].stot;
       inf>>allruns[irun].risetime;
       inf>>allruns[irun].srisetime;
       inf>>allruns[irun].width;
       inf>>allruns[irun].swidth;
       inf>>allruns[irun].chovampl;
       inf>>allruns[irun].schovampl;
       entries++;
    }
    cout<<endl<<"Found "<<entries<<" existing entries in the table file "<<filename<<endl<<endl;
    inf.close();
  }
  /// adding new parameters now!
//  int irun = runpar->runNo - MINRUN;
   int irun = 100*runpar->detNo + runpar->runNo ;
   AssignRun(allruns,runpar,irun);
   FILE* outf = fopen(filename,"w");
   if (outf==NULL){ cout<<"FAILED TO OPEN OUTFILE"<<endl; }
   fprintf(outf,"%6s ",title[0]);
   fprintf(outf,"%6s ",title[1]);
   for (int i=2;i<18;i++) fprintf(outf,"%14s ",title[i]);
   fprintf(outf,"\n");
   for (int i=0;i<TOTRUNS*100;i++)
   {
//      irun = allruns[i].runNo - MINRUN;
     irun = i;
     if (allruns[irun].detNo <= 0 || allruns[irun].runNo<0) continue;
//      cout<<"Writing entry "<<irun<<" corresponding to "<< allruns[irun].runNo<<endl;
     fprintf(outf,"%6d %6d %14e %14e %14e %14e %14e %14e %14e %14e ",allruns[irun].detNo, allruns[irun].runNo, allruns[irun].ampl, allruns[irun].sampl, allruns[irun].charge, allruns[irun].scharge, allruns[irun].rate, allruns[irun].srate, allruns[irun].grate, allruns[irun].sgrate);
     fprintf(outf,"%14e %14e %14e %14e %14e %14e %14e %14e\n",allruns[irun].tot, allruns[irun].stot, allruns[irun].risetime, allruns[irun].srisetime, allruns[irun].width, allruns[irun].swidth, allruns[irun].chovampl, allruns[irun].schovampl);
   }
   fclose(outf);
  
}

void DrawPlot(int run, TH1F *histo, int print)
{
///     const char *txt = histo->GetName();
     const char *txt = histo->GetTitle();
     char txt2[1000];
//      TPostScript *ps;

     sprintf(txt2,"plots%d/%s.eps",run,txt);

     if (print)
     {
//         ps = new TPostScript(txt2,113);
//         ps->NewPage();
     }

     TCanvas* canv = new TCanvas(txt,txt,1);
     canv->SetFillColor(10);
     canv->SetFillColor(10);
     canv->SetBorderSize(2);
     canv->SetGridx();
     canv->SetGridy();
//     canv->SetLogy();
     canv->SetFrameBorderSize(12);
    gStyle->SetOptStat(1000010);
    gStyle->SetPalette(1,0);
//    histo2->SetMinimum(0.1);
    histo->GetXaxis()->SetTitleFont(22);
    histo->GetYaxis()->SetTitleFont(22);
    histo->GetXaxis()->SetLabelFont(22);
    histo->GetYaxis()->SetLabelFont(22);
    histo->GetXaxis()->SetTitleColor(1);
    histo->GetXaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetTitleOffset(1.0);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->SetMarkerStyle(21);
    histo->SetMarkerSize(0.34);
    histo->SetMarkerColor(4);
    histo->SetFillColor(0);
    histo->SetLineColor(4);
    histo->SetLineWidth(1);
    histo->Draw();
    canv->Update();

//     if (print)
//       ps->Close();
}

void DrawPlot2D( TH2D *histo, int print)
{
//     const char *txt = histo->GetName();
     const char *txt = histo->GetTitle();
     char txt2[1000];
//      TPostScript *ps;

        sprintf(txt2,"plots/%s.eps",txt);
     
     if (print)
     {
//         ps = new TPostScript(txt2,113);
//         ps->NewPage();
     }
     
     sprintf(txt2,"%s",txt);

     TCanvas *canv = new TCanvas(txt2,txt,845,800);
     canv->SetFillColor(10);
     canv->SetFillColor(10);
     canv->SetBorderSize(2);
     canv->SetGridx();
     canv->SetGridy();
//     canv->SetLogy();
     canv->SetFrameBorderSize(12);
     gStyle->SetOptStat(1000010);
     gStyle->SetPalette(1,0);
    
//    histo2->SetMinimum(0.1);
    histo->GetXaxis()->SetTitleFont(22);
    histo->GetXaxis()->SetTitleColor(1);
    histo->GetYaxis()->SetTitleFont(22);
    histo->GetXaxis()->SetLabelFont(22);
    histo->GetYaxis()->SetLabelFont(22);
    histo->GetXaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetTitleOffset(1.0);
//    histo->GetXaxis()->SetLabelOffset(0.02);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->SetMarkerStyle(21);
    histo->SetMarkerSize(1.4);
    histo->SetFillColor(0);
    histo->SetLineColor(4);
    histo->SetLineWidth(1);
    histo->Draw("colz");
    
    canv->Update();
  
//     if (print)
//       ps->Close();
}

double distance(double x1,double y1,double x2,double y2)
{
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

float calculateTime(double Rstart, double Rend, int SCpoints, double *lt, int *Condition)
{
	float time=0.;

	for (int i=1; i<SCpoints; i++)
	{
		if (Rstart > lt[i]) continue;
		if (Rend < lt[i]) return (time/60.);
		if (Condition[i] && Condition[i-1])
		{
			if (lt[i]-lt[i-1]<70.)
				time+=(lt[i]-lt[i-1]);
		};
	};
	return (time/60.);
}
float calculateBGTime(double Rstart, double Rend, int SCpoints, double *lt, int *Condition, int *Tracking)
{
	float time=0.;

	for (int i=1; i<SCpoints; i++)
	{
		if (Rstart > lt[i]) continue;
		if (Rend < lt[i]) return (time/60.);
		if (Condition[i] && Condition[i-1] && !Tracking[i])
		{
			if (lt[i]-lt[i-1]<70.)
				time+=(lt[i]-lt[i-1]);
		};
	};
	return (time/60.);
}


void DrawaGraph(TGraph* gr, char* txt, int same)
{

  gr->SetTitle(txt);
  if (same)
  {
    gr->SetMarkerStyle(25);
    gr->SetMarkerSize(.54);
    gr->SetMarkerColor(same);
    gr->SetFillColor(0);
    gr->SetLineColor(same);
    gr->SetLineWidth(1);
    gr->Draw("PL");
  }
  else
  {
     char txt2[2000];
     strcpy(txt2,txt);
     char *ch;
     while ((ch=strrchr(txt2,' '))!=NULL)
       *ch='_';
     
     TCanvas* canv = new TCanvas(txt2,txt2,1);
     canv->SetFillColor(10);
     canv->SetFillColor(10);
     canv->SetBorderSize(2);
     canv->SetGridx();
     canv->SetGridy();
//     canv->SetLogy();
     canv->SetFrameBorderSize(12);
    gStyle->SetOptStat(1000010);
    gStyle->SetPalette(1,0);
    gr->Draw("APL");
    TH1F *histo=gr->GetHistogram();
//    histo->SetMinimum(0.1e-5);
   // histo->SetMaximum(1e2);
    histo->GetXaxis()->SetTitleFont(22);
    histo->GetYaxis()->SetTitleFont(22);
    histo->GetXaxis()->SetLabelFont(22);
    histo->GetYaxis()->SetLabelFont(22);
    histo->GetXaxis()->SetTitleColor(1);
    histo->GetXaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetTitleOffset(1.0);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitle(txt);
     
    gr->SetMarkerStyle(24);
    gr->SetMarkerSize(0.55);
    gr->SetMarkerColor(4);
    gr->SetFillColor(0);
    gr->SetLineColor(4);
    gr->SetLineWidth(1);
    canv->Update();

  }
}


void DrawaGraphErr(TGraphErrors* gr, char* txt, int same, int time=0)
{

  gr->SetTitle(txt);
  if (same)
  {
    gr->SetMarkerStyle(25);
    gr->SetMarkerSize(.54);
    gr->SetMarkerColor(same);
    gr->SetFillColor(0);
    gr->SetLineColor(same);
    gr->SetLineWidth(1);
    gr->Draw("PL");
  }
  else
  {
     char txt2[2000];
     strcpy(txt2,txt);
     char *ch;
     while ((ch=strrchr(txt2,' '))!=NULL)
       *ch='_';
     
     TCanvas* canv = new TCanvas(txt2,txt2,1);
     canv->SetFillColor(10);
     canv->SetFillColor(10);
     canv->SetBorderSize(2);
     canv->SetGridx();
     canv->SetGridy();
//     canv->SetLogy();
     canv->SetFrameBorderSize(12);
    gStyle->SetOptStat(1000010);
    gStyle->SetPalette(1,0);
    gr->Draw("APL");
    TH1F *histo=gr->GetHistogram();
//    histo->SetMinimum(0.1e-5);
   // histo->SetMaximum(1e2);
    histo->GetXaxis()->SetTitleFont(22);
    histo->GetYaxis()->SetTitleFont(22);
    histo->GetXaxis()->SetLabelFont(22);
    histo->GetYaxis()->SetLabelFont(22);
    histo->GetXaxis()->SetTitleColor(1);
    histo->GetXaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetTitleOffset(1.0);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitle(txt);
    if (time>0)
    {
       histo->GetXaxis()->SetTimeDisplay(time);
       histo->GetXaxis()->SetTimeDisplay(time);     
       histo->GetXaxis()->SetTimeFormat("#splitline{%02d/%02m}{%H:%M}");
       histo->GetXaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%d-%m-%y}");
       histo->GetXaxis()->SetLabelOffset(0.021);
       histo->GetXaxis()->SetLabelSize(0.03);
    }
     
    gr->SetMarkerStyle(24);
    gr->SetMarkerSize(0.55);
    gr->SetMarkerColor(4);
    gr->SetFillColor(0);
    gr->SetLineColor(4);
    gr->SetLineWidth(1);
    if (time<-5) 
    {
       canv->BuildLegend();
       cout <<" Build Legend!!!!"<<endl;
    }
    canv->Update();

  }
}


int FitaGraph(double* x, double*y, int nbins, double* par, double endval)
{
   //int npar=0;
   TGraph* gr=new TGraph(nbins,x,y);
   TF1 *f2= new TF1("FitChi2","pol9",0.001,endval);
   gr->Fit(f2,"NB0Q","",0.001,endval);
   for (int i=0;i<f2->GetNpar();i++)
     par[i] = f2->GetParameter(i);
   return (f2->GetNpar());



}

int GetProjectionFitParameters(TH2D* histo2d, double* mean, double* sigma, double* gmean, double* gsigma, double* runs, int maxr)
{    
    TH1D *proj;
//     TH1D* proj2;
    int nbins=histo2d->GetNbinsX();
//     int nybins=histo2d->GetNbinsY();
    //double Energy[dbins], rms[dbins], sigma[dbins];
   
   double nothing=gmean[0];
   nothing=gsigma[0];
    
    int gpoints =0;
// //      TCanvas* c1= new TCanvas("DL_fit","DL_fit");
// //      c1->SetFillColor(10);
// //      c1->SetBorderSize(1);
// //      c1->SetFrameBorderSize(0);
// //      gPad->SetFillColor(10);
// //      gPad->SetBorderSize(1);
// //      gPad->SetGridx();
// //      gPad->SetGridy();
// //      gPad->SetFrameBorderSize(0);
// //      gStyle->SetOptStat(1111);
// //      gStyle->SetOptFit(111);
// //      TH1::StatOverflows(kTRUE);
    proj->StatOverflows(kTRUE);
    for (int n=1; n<nbins; n++)
    {
       
       proj = histo2d->ProjectionY("proj",n,n+1,"");
       
       double runNo=histo2d->GetXaxis()->GetBinLowEdge(n);
       if (runNo>maxr) break;
       
       if (proj->Integral()<=0) continue;
       
       runs[gpoints] = runNo;
       

       sigma[gpoints] = proj->GetRMS();
       if ( sigma[gpoints]<=0. ) continue;       
       
       mean[gpoints] = proj->GetMean();
       if ( mean[gpoints]<=0. ) continue;       
       
       cout <<"fitting slice " <<n<<endl;///proj->GetMean()<< endl;;
     
//        double miny=proj->GetBinLowEdge(1);
//        double maxy=proj->GetBinLowEdge(nbins-1);

////////////gaus is here!!!!!!       
// //        TF1* f1 = new TF1("gauss","gaus",miny,maxy);
// //        proj->Fit(f1,"NQR");
// //        
// // 
// //        gsigma[gpoints] = f1->GetParameter(2);
// //        
// //        gmean[gpoints] = f1->GetParameter(1);
// //        if ( gmean[gpoints]<=0 ) continue;       

///       gStyle->SetOptStat(1111111);
//         c1->Update();
///     cin >> n;
       gpoints++;
     }
    
    return (gpoints);   
}

void DrawGraphsFrom(TH2D* pthc, TH2D* pth, int maxr, char* txt1, char* txt2)
{
  
  double mean[10000], sigma[10000], gmean[10000], gsigma[10000], fitruns[10000];
  double bmean[10000], bsigma[10000], bgmean[10000], bgsigma[10000], bfitruns[10000];
  double err0[10000];
  char txt[1000];
  
  int ntst=GetProjectionFitParameters(pthc, mean, sigma, gmean, gsigma, fitruns, maxr);
  int bntst=GetProjectionFitParameters(pth, bmean, bsigma, bgmean, bgsigma, bfitruns, maxr);
  
  
  TGraphErrors* gr1= new TGraphErrors(ntst,fitruns,mean,err0,err0);
//   TGraphErrors* gr2= new TGraphErrors(ntst,fitruns,gmean,err0,err0);
  TGraphErrors* bgr1= new TGraphErrors(bntst,bfitruns,bmean,err0,err0);
//   TGraphErrors* bgr2= new TGraphErrors(bntst,bfitruns,bgmean,err0,err0);
  
  double max1=TMath::MaxElement(ntst,mean);
//   double max2=TMath::MaxElement(ntst,gmean);
  double max1b=TMath::MaxElement(bntst,bmean);
//   double max2b=TMath::MaxElement(bntst,bgmean);
  double max = TMath::Max(max1,max1b);
//   max1 = TMath::Max(max2,max2b);
//   max = TMath::Max(max1,max);
  
  double min1=TMath::MinElement(ntst,mean);
//   double min2=TMath::MinElement(ntst,gmean);
  double min1b=TMath::MinElement(bntst,bmean);
//   double min2b=TMath::MinElement(bntst,bgmean);
  double min = TMath::Min(min1,min1b);
//   min1 = TMath::Min(min2,min2b);
//   min = TMath::Min(min1,min);
  
  TH1F *htmp=(TH1F*) gr1->GetHistogram();
  htmp->SetMaximum(1.1*max);
  htmp->SetMinimum(0.9*min);
  
  sprintf(txt,"%s",txt1);
  DrawaGraphErr(gr1, txt, 0, 0);
  sprintf(txt,"Gaus %s",txt1);
//   DrawaGraphErr(gr2, txt, 3, 0);
 
  sprintf(txt,"%s bkg",txt1);
  DrawaGraphErr(bgr1, txt, 2, 0);
  sprintf(txt,"Gaus %s bkg",txt1);
//   DrawaGraphErr(bgr2, txt, 6, -13);
  
      TGraph *graph = new TGraph(4);
      graph->SetName("LowBkgRuns");
      graph->SetTitle("Low bkg Runs");
      graph->SetFillColor(5);
      graph->SetFillStyle(3005);
      graph->SetPoint(0,15075.,0.9*min);
      graph->SetPoint(1,15075.,1.1*max);
      graph->SetPoint(2,15110.,1.1*max);
      graph->SetPoint(3,15110.,0.9*min);
     /// graph->Draw("f");
    
    graph->Draw("F");
  
  
  gr1= new TGraphErrors(ntst,fitruns,sigma,err0,err0);
//   gr2= new TGraphErrors(ntst,fitruns,gsigma,err0,err0);
  bgr1= new TGraphErrors(bntst,bfitruns,bsigma,err0,err0);
//   bgr2= new TGraphErrors(bntst,bfitruns,bgsigma,err0,err0);
  
  max1=TMath::MaxElement(ntst,sigma);
//   max2=TMath::MaxElement(ntst,gsigma);
  max1b=TMath::MaxElement(bntst,bsigma);
//   max2b=TMath::MaxElement(bntst,bgsigma);
  max = TMath::Max(max1,max1b);
//   max1 = TMath::Max(max2,max2b);
//   max = TMath::Max(max1,max);
  
  min1=TMath::MinElement(ntst,sigma);
//   min2=TMath::MinElement(ntst,gsigma);
  min1b=TMath::MinElement(bntst,bsigma);
//   min2b=TMath::MinElement(bntst,bgsigma);
  min = TMath::Min(min1,min1b);
//   min1 = TMath::Min(min2,min2b);
//   min = TMath::Min(min1,min);
  
  htmp=(TH1F*) gr1->GetHistogram();
  htmp->SetMaximum(1.1*max);
  htmp->SetMinimum(0.9*min);
  
  sprintf(txt,"%s",txt2);
  DrawaGraphErr(gr1, txt, 0, 0);
  sprintf(txt,"Gaus %s",txt2);
//   DrawaGraphErr(gr2, txt, 3, 0);
 
  sprintf(txt,"%s bkg",txt2);
  DrawaGraphErr(bgr1, txt, 2, 0);
  sprintf(txt,"Gaus %s bkg",txt2);

      TGraph *graph2 = new TGraph(4);
      graph2->SetName("LowBkgRunsW");
      graph2->SetTitle("Low bkg Runs");
      graph2->SetFillColor(5);
      graph2->SetFillStyle(3005);
      graph2->SetPoint(0,15075.,0.9*min);
      graph2->SetPoint(1,15075.,1.1*max);
      graph2->SetPoint(2,15110.,1.1*max);
      graph2->SetPoint(3,15110.,0.9*min);
      graph2->SetFillStyle(3004);
     /// graph->Draw("f");
      graph2->Draw("F");
//   DrawaGraphErr(bgr2, txt, 6, -13);
}  

int GetProjectionTParameters(TH2D* histo2d, double* mean, double* sigma, double* gmean, double* gsigma, double* runs)
{    
    TH1D *proj;
///    TH1D* proj2;
    int nbins=histo2d->GetNbinsX();
    int nybins=histo2d->GetNbinsY();
    nybins+=0*gsigma[0];
    //double Energy[dbins], rms[dbins], sigma[dbins];
   
  TCanvas *c1= new TCanvas("demo","demo");
  
    
    int gpoints =0;
    
  //  proj->StatOverflows(kFALSE);
    for (int n=1; n<=nbins; n++)
    {
       
       proj = histo2d->ProjectionY("proj",n,n+1,"");
       
       double dt=histo2d->GetXaxis()->GetBinWidth(1);
       
       double runNo=histo2d->GetXaxis()->GetBinLowEdge(n)+dt/2.;
       
       proj->GetXaxis()->SetRange(2,proj->GetNbinsX()-2); 
       
       if (proj->Integral()<=0) continue;
       
//        proj->Sumw2();
       runs[gpoints] = runNo;
       gmean[gpoints] = dt/2.;
       

       sigma[gpoints] = proj->GetRMS();
       if ( sigma[gpoints]<=0. ) continue;       
       
       mean[gpoints] = proj->GetMean();
       if ( mean[gpoints]<=0. ) continue;       
    
       double miny= mean[gpoints]-1.5*sigma[gpoints];
       double maxy= mean[gpoints]+1.2*sigma[gpoints];

       double sumy=0,sumxy=0;
       for (int i=1;i< proj->GetNbinsX()-1; i++)
       {
           sumy+=proj->GetBinContent(i);
	   sumxy+=proj->GetBinContent(i)+proj->GetBinCenter(i); 
       }
       mean[gpoints] = sumxy/sumy;
       
          
       cout <<"fitting slice " <<n<<endl;///proj->GetMean()<< endl;;
     
//        double miny=proj->GetBinLowEdge(1);
//        double maxy=proj->GetBinLowEdge(nbins-1);
/// /////////gaus is here!!!!!!       
       TF1* f1 = new TF1("gauss","gaus",miny,maxy);
       proj->Fit(f1,"QRW0");
       
       sigma[gpoints] = (2.3*f1->GetParameter(2)) / f1->GetParameter(1);
//        sigma[gpoints] = (f1->GetParameter(2));
       
       mean[gpoints] = f1->GetParameter(1);///*174.94;
       if ( gmean[gpoints]<=0 ) continue;       

       proj->SetMarkerStyle(1);
       proj->SetMinimum(0);
       gStyle->SetOptStat(1111111);
          c1->Update();
// int tst;
// cin>>tst;
//      cin >> n;
       gpoints++;
     }
//     delete proj;
    return (gpoints);   
}

void DrawXsliceFitsFrom(TH2D* pthc, char* txt1, char* txt2, int timef, double*mean, double* sigma, double* times)
{
  double  xerr[10000], gsigma[10000], fitruns[10000];
  double err0[10000];
  char txt[1000];
  
  int ntst=GetProjectionTParameters(pthc, mean, sigma, xerr, gsigma, fitruns);
  
  TGraphErrors* gr1= new TGraphErrors(ntst,fitruns,mean,xerr,err0);
  TGraphErrors* gr2= new TGraphErrors(ntst,times,mean,xerr,sigma);
  
  sprintf(txt,"%s",txt1);
  DrawaGraphErr(gr1, txt, 0, timef);
  sprintf(txt,"%s",txt1);
  DrawaGraphErr(gr2, txt, 0, 1);
 
  
  
  gr1= new TGraphErrors(ntst,fitruns,sigma,xerr,err0);
  gr2= new TGraphErrors(ntst,times,sigma,xerr,err0);
  
  
  sprintf(txt,"%s_graph",txt2);
  DrawaGraphErr(gr1, txt, 0, timef);
  sprintf(txt,"%s",txt2);
  DrawaGraphErr(gr2, txt, 0, 1);
 


}

double SumArray(double* ampl,int ni,int nf)
{
   double s=0.;
   for (int i=ni;i<=nf;i++)
     s+=ampl[i];
   return s;
}

void SmoothArray(double *arr, double *sarr, int apoints, int np, int inverse)
{
   if (np<=0) np=1;
   np = (np / 2) * 2 + 1; // make it odd >= np
   for (int i=0;i<np/2-1;i++)
     sarr[i]=arr[i]*inverse;
   double sum = SumArray(arr,0,np-1);
   sarr[np/2]=sum/np*inverse;
   
   for (int i=np/2+1;i<apoints-np/2;i++)
   {
     sum -= arr[i-(np/2+1)];
     sum += arr[i+np/2];
     sarr[i] = sum/np*inverse;

//      cout<<i<<endl;
  }
   for (int i=np/2+1-1;i>=0;i--)
   {
     sarr[i] = sarr[i+1];
   }
   for (int i=apoints-np/2;i<apoints;i++)
   {
     sarr[i] = sarr[i-1];
   }

}

void DerivateArray(double *arr, double *sarr, int apoints, double dt, int np, int neg=1)
{
   if (neg<0)
     neg=-1;
   else
     neg=1;
   if (np>20) np=20;
   if (np<1) np=1;  /// 1 <= np <= 10
//    dt*=1e9;
//    double sum = 0.;
   for (int i=0;i<np;i++)
     sarr[i]=0.;
   for (int i=apoints-np;i<apoints;i++)
     sarr[i]=0.;
   
   for (int i=np;i<apoints-np;i++)
//    for (int i=0;i<apoints-np;i++)
   {
     sarr[i] = neg*(arr[i+np] - arr[i-np])/(2*np*dt);
//      sarr[i] = (arr[i+np] - arr[i-np])/(2.*np);
   }
}

void DerivateArrayT(double *arr, double *sarr, int apoints, double dt, double DT, int neg=1)
{ 
    /// with the time window
   if (neg<0)
     neg=-1;
   else
     neg=1;
   int np = TMath::FloorNint(DT/dt);
   if (np>12) np=12;
   if (np<1) np=1;  /// 1 <= np <= 12
   
//    double sum = 0.;
   for (int i=0;i<np;i++)
     sarr[i]=0.;
   for (int i=apoints-np;i<apoints;i++)
     sarr[i]=0.;
   
   for (int i=np;i<apoints-np;i++)
//    for (int i=0;i<apoints-np;i++)
   {
     sarr[i] = neg*(arr[i+np] - arr[i-np])/(2*np*dt);
//      sarr[i] = (arr[i+np] - arr[i-np])/(2.*np);
   }

}

double FindFraction(TSpline *spline, double zero, double min, double ampl, double fraction)
{ 
   const double precision = 1e-6;  /// 1 ps
   double mid = (zero+min)/2.;
//    if (ampl>1.) ampl = 1.;
   if (fabs(min-zero)<precision)    // Check if the difference between min and zero is below precision
      return (mid);
   if (fabs(spline->Eval(min) - ampl*fraction)<precision)    // Check if the spline value at min is close to the desired                   fraction of the amplitude
      return (min);
   else if (fabs(spline->Eval(zero) - ampl*fraction)<precision)    // Check if the spline value at zero is close to the desired fraction of the amplitude
      return (zero);
   
   if (spline->Eval(mid) < ampl*fraction)    // If the spline value at mid is less than the desired fraction, search in the upper half
   {
       min = mid;
       double frac = FindFraction(spline, zero, min, ampl, fraction);
       return (frac);
   }
   else    // If the spline value at mid is greater than or equal to the desired fraction, search in the lower half
   {
       zero = mid;
       double frac = FindFraction(spline, zero, min, ampl, fraction);
       return (frac);
   }
} 

void FindTimesS(TGraph* graph, DPARAM *par)
{
   double *data = graph->GetY();
   double *T = graph->GetX();
   int points = graph->GetN();
   
   /// find maximum of pulse derivative...
   int max=TMath::LocMax(points,data);
   par->maxtime=T[max];
   par->ampl=data[max];
   double ampl=data[max];

   /// find minimum of pulse derivative...
   int min=TMath::LocMin(points,data);
   double amplmin=data[min];
   par->t3 = T[min];
 
//    cout <<"____________________________\n min point at: "<<min<<endl;
   ///Find start of pulse
   int zero=max;
   while (data[zero]>0. && zero>10) zero--;
   par->t1 = T[zero];
//    cout <<" zero crossing at: "<<zero<<endl;
   /// Find end of Pulse
   int zero2=max;
   while (data[zero2]>0. && zero2<points) zero2++;
   
   int zero3=min;
   while (data[zero3]<0. && zero3>1)
   {
     zero3--;
     if (fabs(data[zero3])<fabs(data[min])*0.2) 
       break;
   }
   par->t2=T[zero3];
   
///   TGraph *graph = new TGraph(points,t,data);
   TGraph *graphS = new TGraph(zero2-zero+1,&T[zero],&data[zero]);
   TSpline *spline = new TSpline3("spline",graphS);
   double smax = 0.;
   for (double x = T[zero]; x< T[zero2]; x+=0.0005)  //times in ns!
   {
     if (spline->Eval(x)>smax)
     {
       par->maxtimeS = x;
       smax = spline->Eval(x);
       par->sampl=smax;
     }
   }
//    char *ch = "test";
//    DrawaGraph(graphS, ch, 0);
//    spline->Draw("same");

   for (double x = T[zero2]; x > par->maxtimeS; x-=0.0005)
   {
     if (spline->Eval(x)*spline->Eval(x-0.0005)<=0.)
     {
       par->ftime = x-0.00025;
       break;
     }
   }
   double sum=0.;
   for (double x = T[zero]; x <= par->ftime; x+=0.0005)
     sum += spline->Eval(x);
   
   sum*=0.0005 ;
   
   par->intg = sum;
  
  
///   double miny = f1->Eval(mean);
   
   par->time50 = FindFraction(spline,T[zero],T[max],ampl, 0.5);
///    cout<<"T(50%) = "<<par[2] <<" ampl = "<<ampl<<" tmax = "<< par[0]<<"   spline tmax = "<< par[3]<<"   spline t20 = "<< par[4]<<endl;
   delete spline;
   delete graphS;
///   delete graph;

}

void FindTimesSD(TGraph* graph, DPARAM *par, double tstart, double tend, double *tf)
{
   double *data = graph->GetY();
   double *T = graph->GetX();
   double dt = T[1]-T[0];
   int points = graph->GetN();
   
   /// find maximum of pulse
   int max=TMath::LocMax(points,data);
   par->maxtime=T[max];
   par->ampl=data[max];
   double ampl=data[max];

   ///Find start of pulse
   int istart=max;
   while (data[istart] > par->rms && istart>10) istart--;
   par->t1 = T[istart];

   /// find end of pulse 
   int iend = max;
   while (data[iend] > par->rms && iend<points) iend++;
   par->t2 = T[iend];
   
   int spoints = 0;
   int ppoints = 0;
   double intg = 0.;
   double charge = 0.;
   int its=0,ite=0;
   int Start=1, End=1;   
   for(int i=0; i<points; i++)
   {
     if (Start && T[i]>=tstart)
     {
       its = i;
       Start = 0;
     }
     if (End && T[i]>tend)
     {
       ite = i;
       End = 0;
     }
     if ((!Start) && End)
     {
       intg += data[i];
       spoints++;
     }
     if (i>=istart && i<=iend)
     {
       charge += data[i];
       ppoints++;
     }
   }

///   recalculate the maximum from the spline (maxtimeS and smax)
   TGraph *graphS = new TGraph(spoints,&T[its],&data[its]);
   TSpline *spline = new TSpline3("spline",graphS);
   double smax = 0.;
   for (double x = T[istart]; x< T[iend]; x+=0.0005)  //times in ns!
   {
     if (spline->Eval(x)>smax)
     {
       par->maxtimeS = x;
       smax = spline->Eval(x);
       par->sampl=smax;
     }
   }

   par->intg = intg*dt/(par->ampl) *0.1;//*dt;
   par->charge = charge/80.; /// (normalize to ~80 points)

//    char *ch = "test";
//    DrawaGraph(graphS, ch, 0);
//    spline->Draw("same");

//    for (double x = T[istart2]; x > par->maxtimeS; x-=0.0005)
//    {
//      if (spline->Eval(x)*spline->Eval(x-0.0005)<=0.)
//      {
//        par->ftime = x-0.00025;
//        break;
//      }
//    }

   if (fabs(smax-par->ampl)/(smax+par->ampl)>0.03) smax=par->ampl;
   
   double t0 = par->maxtimeS;
   double t1 = par->t1;
   double dx = 0.00025;
    
   for (int i=6; i>=0; i--)
   {
     double ampl0 = smax * (i+1.) / 10.;
     for (double x = t0; x>t1; x-=dx)
     {
       double y1 = spline->Eval(x-dx) - ampl0;
       double y2 = spline->Eval(x) - ampl0;
       if ( y1*y2 <=0. )
       {
        if (y1==0)
            tf[i] = x-dx;
        else
            tf[i] = x-dx + dx * y1/(fabs(y1)+fabs(y2));
        
        t0 = tf[i];
    // 	 cout<<tf[i]<<" ";
        break;
       }
     }
   }
/// find the threshold crossing point
   par->ttrig=100.;
   double ampl0 = Threshold*1.75;
   for (double x = par->maxtimeS; x>t1; x-=dx)
     {
       double y1 = spline->Eval(x-dx) - ampl0;
       double y2 = spline->Eval(x) - ampl0;
       if ( y1*y2 <=0. )
       {
	 if (y1==0)
	   par->ttrig = x-dx;
	 else
	 par->ttrig = x-dx + dx * y1/(fabs(y1)+fabs(y2));
	 
// 	 cout<<tf[i]<<" ";
	 break;
       }
     }


  
///   double miny = f1->Eval(mean);
   
//    par->time50 = FindFraction(spline,T[istart],T[max],ampl, 0.5);
///    cout<<"T(50%) = "<<par[2] <<" ampl = "<<ampl<<" tmax = "<< par[0]<<"   spline tmax = "<< par[3]<<"   spline t20 = "<< par[4]<<endl;
   delete spline;
   delete graphS;
///   delete graph;

}

void FindTimesF(TF1* f1,double t1,double t2,double* tf)
{
   double ampl = f1->GetMaximum(t1,t2);
   double xmax = f1->GetMaximumX(t1,t2);
   double dx = 0.00025;
   double t0 = xmax;
   double y1,y2;
   for (int i=6; i>=0; i--)
   {
     double ampl0 = ampl * (i+1.) / 10.;
     for (double x = t0; x>t1; x-=dx)
     {
       y1 = f1->Eval(x-dx) -ampl0;
       y2 = f1->Eval(x) -ampl0;
       if ( y1*y2 <=0. )
       {
	 if (y1==0)
	   tf[i] = x-dx;
	 else
	 tf[i] = x-dx + dx * y1/(fabs(y1)+fabs(y2));
	 
	 t0 = tf[i];
	 break;
       }
     }
   }
  
}

void FindTimesG(TGraph* gr1,double t1,double t2,double* tf)
{
   double dx = 0.00025;
   double y1,y2;

   double *data = gr1->GetY();
   double *T = gr1->GetX();
   int points = gr1->GetN();
   int npoints = 0;
   int ts=0,te=0;
   for (int i=0; i<points; i++)
   {
     if (T[i]>=t1-dx){
       ts=i;
       break;
     }
   }
   for (int i=ts; T[i]<t2; i++)
   {
     te=i;
     npoints++;
   }
       
   TGraph *graphS = new TGraph(npoints,&T[ts],&data[ts]);
   TSpline *f1 = new TSpline3("spline",graphS);

   double ampl = f1->Eval(t1);
   double xmax = t1;

   for (double x = t1; x<t2; x+=dx)
   {
     y1=f1->Eval(x);
     if (y1>ampl)
     {
       ampl=y1;
       xmax=x;
     }
   }
//       char ch[20];
//       sprintf(ch,"testTrig");
//       DrawaGraph(graphS, ch, 0);
//       f1->Draw("same");
//   return ;   
   double t0 = xmax;
   
   for (int i=6; i>=0; i--)
   {
     double ampl0 = ampl * (i+1.) / 10.;
     for (double x = t0; x>t1; x-=dx)
     {
       y1 = f1->Eval(x-dx) -ampl0;
       y2 = f1->Eval(x) -ampl0;
       if ( y1*y2 <=0. )
       {
	 if (y1==0)
	   tf[i] = x-dx;
	 else
	 tf[i] = x-dx + dx * y1/(fabs(y1)+fabs(y2));
	 
	 t0 = tf[i];
// 	 cout<<tf[i]<<" ";
	 break;
       }
     }
   }
//    cout<<endl;
   delete f1; 
   delete graphS;
}


void IntegrateG(TGraph* graph, DPARAM *par)
{
 par->charge = 0.;
 if (graph) {
   double *data = graph->GetY();
   double *T = graph->GetX();
   int points = graph->GetN();
//    cout<<" Graph has "<<points <<" points"<<endl;
   double tstart = par->t1;  
   int t0=0;
   while (T[t0]<tstart && t0<points) {t0++;}
   
   double sum = 0.;
   for (int i=t0 ; T[i]<par->t3 || data[i]>0.; i++)
   {
     if (i> points)
     {
//        cout<<" end of points reached"<<endl;
//        cout <<i <<endl;
//        cout  <<endl;
       break;
     }
     double dt = T[i+1]-T[i];
     sum += (data[i]+data[i+1]) / (2.*dt);
   }
   par->charge = sum;
 }
}

double IntegrateA(int npoints, double* data, double* integral,double dt)
{
  double sum=0.;
//   dt*=1e9;
  if (npoints<4) return sum;
  
  for (int i=0;i<4;i++) integral[i]=0;
  
  for (int i=4;i<npoints;i++)
  {
    integral[i] = integral[i-4] + 2*(data[i-3]+data[i-1])*dt;
    sum+=integral[i];
  }
  return (sum);
}
    
double IntegratePulse(int npoints, const double* data, double* integral,double dt, double tint)
{
//   dt*=1e9;
  double sum=data[0]*dt;
  if (npoints<2) return sum;
  
  int nint = TMath::FloorNint(tint/dt); //round to the nearest integer of the tint/dt
  ///tint is a set time window for integration 
  integral[0]=data[0]*dt;
  if (nint<2) nint = 2;
 
  for (int i=1; i<nint; i++){
    integral[i]=data[i]*dt+integral[i-1];
    //this could be also integral[i] = integral[i]+integral[i-1]; ????
    sum+=data[i]*dt;
  }
  
  for (int i=nint; i<npoints; i++){
    integral[i]=data[i]*dt+integral[i-1]-data[i-nint]*dt;
    sum+=data[i]*dt;
  }
  return (sum);
}
    
  
int AnalyseIntegratedPulse(int points, double* data, IPARAM *par, double threshold, double dt)
{
  /// use the integrated+filtered pulse to define a region where a trigger occured. (integral above threshold) 
  ///pulses are considered negative!!!
//   dt*=1e9;
  if (points < 20) return -1;
  int ntrig=0;
  int tpoint=0;
  par->tot=0;
  for (int i=0; i<points; i++)   {
    if (data[i]<=threshold) {
      tpoint = i;
      ntrig=1;
      par->tot=1;
      break;
    }
  }
//   cout<<"tpoint = "<<tpoint*dt<<endl;
  if (tpoint>=points-10 || tpoint == 0) return 0;

  double miny = data[tpoint];
  par->charge=0.;

  for (int i=tpoint; i<points; i++)
  {
    if (data[i]<miny)
    {
      par->ampl=data[i];
      par->maxtime=i;
      miny=data[i];
    }
    par->charge+=data[i];
    if (data[i]<=threshold)
      par->tot++;
    if (data[i]>threshold/10.)
    {
      par->ftime=i;
      break;
    }
  }
  /// fast scan for risetime, risecharge and t_start
  par->t90=tpoint;
  par->t10=tpoint;
  par->stime=tpoint;
  for (int i=par->maxtime; i>0; i--)
  {
     if (data[i]>=par->ampl*0.9)
     {
       par->t90=i;
       break;
     }
  }
  par->risecharge=0.;
  for (int i=par->t90; i>0; i--)
  {
    par->risecharge+=data[i]; 
    if (data[i]>=par->ampl*0.1)
    {
      par->t10=i;
      break;
    }
  }
  
  for (int i=par->t10; i>0; i--)
  {
    if (i<tpoint)
      par->charge+=data[i];
      par->stime=i;
    if (data[i]>threshold/100. || (data[i+1]-data[i-1])/dt >=-0.0001 )
    {
      break;
    }
  }
  par->width = par->ftime-par->stime;
  par->ampl*=-1.;
  par->charge*=-1.*dt;
//   cout<<"tstart = "<<par->stime*dt<<endl;
//   cout<<"tend = "<<par->ftime*dt<<endl;
//   cout<<"tmax = "<<par->maxtime*dt<<endl;
//   cout<<"ampl = "<<par->ampl<<endl;
//   cout<<"t10 = "<<par->t10*dt<<endl;
//   cout<<"t90 = "<<par->t90*dt<<endl;
//   cout<<"rt = "<<(par->t90-par->t10)*dt<<endl;
//   cout<<"charge = "<<par->charge*dt/N_INTEGRATION_POINTS<<endl;
//   cout<<"risecharge = "<<par->risecharge*dt/N_INTEGRATION_POINTS<<endl;
  return (ntrig);
}
      
int SubtractRefferenceChannel(int* spoints,double** data, double* sampl,int cref, double* bsl)
{
  int points = spoints[cref];
  if (points<0)
  {
    cout<<"Pulse to subtract is empty!!!"<<endl;
    return -1;
  }
  double f[] = {0.,0.,0.,0.};
  for (int i=0;i<4;i++)
  {
    if (i!=cref && spoints[i]>0)
    {
      f[i]++;
      f[cref]--;
    }
  }
    
  for (int i=0;i<points;i++)
  {
      sampl[i]=f[0]*(data[0][i]-bsl[0]) + f[1]*(data[1][i]-bsl[1]) + f[2]*(data[2][i]-bsl[2]) + f[3]*(data[3][i]-bsl[3]); 
  }
  
  return (-f[cref]);
}
 
int SubtractMaxChannel(int points,double** data, double* sampl, double* bsl, int* mask)
{
  if (points<=0)
  {
    cout<<"Nothing to subtract!!!"<<endl;
    return -1;
  }
//   cout<<"points  = "<<points<<endl;
  double max;  
  for (int i=0;i<points;i++)
  {
    max = -1000.;
    for (int ci=0;ci<4;ci++)
    {
      if (mask[ci]) continue;
      data[ci][i] -= bsl[ci];
      if (data[ci][i]>max) 
	max = data[ci][i];
    }
    sampl[i]=0;
//     if (i<100) 
//      cout<<" max = "<<max<<endl;
    for (int ci=0;ci<4;ci++)
    {
      if (mask[ci]) continue;
      data[ci][i]-=max;
      sampl[i]+=data[ci][i];
    }
//     if (i<30) cout<<i<<endl;
  }
  
  return (max);
}
 
int SubtractTwoChannels(int points,double** data, double* sampl, double* bsl, int* mask)
{
  if (points<=0)
  {
    cout<<"Nothing to subtract!!!"<<endl;
    return -1;
  }
//   cout<<"points  = "<<points<<endl;
//   double max;  
  for (int i=0;i<points;i++)
  {
    int i1 = mask[0];
    int i2 = mask[1];
    double y1=data[i1][i]-bsl[i1];
    double y2=data[i2][i]-bsl[i2];
    if (y2-y1<0)
    {
      sampl[i]=y2-y1;
    }
    else
      sampl[i]=y1-y2;
  }
//     max = -1000.;
//     for (int ci=0;ci<2;ci++)
//     {
//       data[mask[ci]][i] -= bsl[ci];
//       if (data[mask[ci]][i]>max) 
// 	max = data[mask[ci]][i];
//     }
//     sampl[i]=0;
// //     if (i<100) 
// //      cout<<" max = "<<max<<endl;
//     for (int ci=0;ci<2;ci++)
//     {
//       data[mask[ci]][i]-=max;
//       sampl[i]+=data[mask[ci]][i];
//     }
// //     if (i<30) cout<<i<<endl;
//   }
  
  return (points);
}
 
int SubtractBaselineOld(int points,double** data, double* sampl, double* bsl, int* mask)
{
  if (points<=0)
    {
      cout<<"Nothing to subtract!!!"<<endl;
      return -1;
    }
//   cout<<"points  = "<<points<<endl;
//   double max;  
  for (int i=0;i<points;i++)
    {
      int i1 = mask[0];
      int i2 = mask[1];
      double y1=data[i2][i]-bsl[i2];
      sampl[i]=y1;
    }
  
  return (points);
} 

int SubtractBaseline(int points,double* data, double* sampl, double bsl)
{
  if (points<=0)
    {
      cout<<"Nothing to subtract!!!"<<endl;
      return -1;
    }
//   cout<<"points  = "<<points<<endl;
//   double max;  
  for (int i=0;i<points;i++)
    {
      double y1=data[i]-bsl;
      sampl[i]=y1;
    }
  
  return (points);
} 

void Pallette1()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

   Double_t Red[3]    = { 1.00, 0.00, 0.00};
   Double_t Green[3]  = { 0.00, 1.00, 0.00};
   Double_t Blue[3]   = { 1.00, 0.00, 1.00};
   Double_t Length[3] = { 0.00, 0.50, 1.00 };

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}

void Pallette2()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

   Double_t Red[3]    = { 1.00, 0.50, 0.00};
   Double_t Green[3]  = { 0.50, 0.00, 1.00};
   Double_t Blue[3]   = { 1.00, 0.00, 0.50};
   Double_t Length[3] = { 0.00, 0.50, 1.00 };

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}


int AnalyseLongPulse(int points, double* data, double* drv, PEAKPARAM *par, double threshold, double dt, int tshift)
{
  /// use the integrated+filtered pulse to define a region where a trigger occured. (integral above threshold) 
  ///pulses are considered negative!!!
  if (points - tshift < 20) return -1;
  int ntrig=0;
  int tpoint=0;
  double drvtrig = 0.00002;
  if (threshold <0.0025)
    drvtrig = 0.000005;
  par->tot[0]=0;
  
  for (int i=tshift; i<points-20; i++)
  {
    if (data[i]<=threshold)
    {
      tpoint = i;
      ntrig=1;
//       par->tot[0]=1;
      break;
    }
    else 
      tpoint=i;
  }
  if (ntrig<=0) return (-1); // cout<<"No trigger in event!"<<endl;
  
//    cout<<"tpoint = "<<tpoint*dt<<endl;
  if (tpoint>=points-10) return (-1);

  double miny = data[tpoint];
  par->charge=0.;
  par->maxtime_pos=tpoint;
  par->ampl=data[tpoint];

  for (int i=tpoint; i<points; i++)
  {
    if (data[i]<miny)
    {
      par->ampl=data[i];
      par->maxtime_pos=i;
      miny=data[i];
    }
    if (data[i]<=threshold)
    {
      par->tot[0]++;
      par->ftime_pos=i; /// this is added to avoid a pulse at the end of the data that does not return to 0!!!
    }
    else //if (data[i]>threshold)
    {
      par->ftime_pos=i;
      par->tot[0]--;
      break;
    }
    /// note down the point the signal has gone above the threshold
  }
  /// fast scan for risetime, risecharge and t_start
  par->t90=tpoint;
  par->t10=tpoint;
  par->stime_pos=tpoint;
  par->ttrig=tpoint; 
  for (int i=par->maxtime_pos; i>0; i--)
  {
     if (data[i]>=par->ampl*0.9)
     {
       par->t90=i;
       break;
     }
  }
  par->risecharge=0.;
  for (int i=par->t90; i>0; i--)
  {
    par->risecharge+=data[i]; 
    if (data[i]>=par->ampl*0.1 || (data[i]>0.5*threshold && fabs(drv[i])<=drvtrig ))
    {
      par->t10=i;
      break;
    }
  }

  for (int i=par->maxtime_pos; i<points; i++)
  {
    if (data[i]>=par->ampl*0.1 || (data[i]>0.5*threshold && fabs(drv[i])<=drvtrig ))
    {
      par->tb10=i;
      break;
    }
  }
  
  for (int i=par->t10; i>0; i--)
  {
    cout<<data[i]<<"  "<<fabs(drv[i])<<"   "<<threshold<<endl;
    par->stime_pos=i;
//     if (i<tpoint)
//     {
//       par->charge+=data[i];
//     }
    if (data[i]>threshold/5. || ((fabs(drv[i])<=drvtrig && data[i]>threshold*0.8 ) ) )
    {
      par->stime_pos=i;
      break;
    }
  }
  for (int i=par->ftime_pos; i<points; i++)
  {
    par->ftime_pos=i;
    if (data[i]>threshold/5. || (fabs(drv[i])<=drvtrig && data[i]>threshold*0.8 ) )
    {
      break;
    }
  }
  STARTPOINT start_point;
  start_point.pos=par->stime_pos;
  
  par->charge=0.;
  for (int i=par->stime_pos;i<=par->ftime_pos;i++)
    par->charge+=data[i];

  par->sampl = data[par->stime_pos];
  par->fampl = data[par->ftime_pos];
  par->bslch = -0.5 * (data[par->stime_pos] + data[par->ftime_pos])*(par->ftime_pos - par->stime_pos +1.);
  par->width = par->ftime_pos-par->stime_pos;
  par->ampl*=-1.;
  par->charge*=-1.*dt;   ///charge is calculated in V * ns. 
//   cout<<"tstart = "<<par->stime_pos*dt<<endl;
//   cout<<"t10 = "<<par->t10*dt<<endl;
//   cout<<"t90 = "<<par->t90*dt<<endl;
//   cout<<"tmax = "<<par->maxtime_pos*dt<<endl;
//   cout<<"tb10 = "<<par->tb10*dt<<endl;
//   cout<<"tend = "<<par->ftime_pos*dt<<endl;
//   
//   cout<<"rt = "<<(par->t90-par->t10)*dt<<endl;
//   cout<<"tot = "<<(par->tot)*dt<<endl;
//   cout<<"DT = "<<(par->tb10-par->t10)*dt<<endl;
//   cout<<"DTall = "<<(par->ftime_pos-par->stime_pos)*dt<<endl;
//   
//   cout<<"ampl = "<<par->ampl<<endl;
//   cout<<"charge = "<<par->charge*dt/N_INTEGRATION_POINTS<<endl;
//   cout<<"risecharge = "<<par->risecharge*dt/N_INTEGRATION_POINTS<<endl;
  return (par->ftime_pos);
}

GLOBALMAXIMUM FindGlobalMaximum(double *arr, int arrpoints, double *arrt)
{        
    GLOBALMAXIMUM thisMaximum;
    thisMaximum.ampl = 111111;
    thisMaximum.pos = 5; 
    thisMaximum.x = -3;
    for (int i = 0; i < arrpoints ; i++) 
    { 
      if (arr[i] < thisMaximum.ampl)
      {
        thisMaximum.ampl = arr[i];
        //cout<< RED << GLOBALMAXIMUM.ampl<<endlr;
        thisMaximum.pos = i;
        thisMaximum.x = arrt[i];
        //cout<<BLUE<< GLOBALMAXIMUM.pos<<endlr; 
      }

      //cout<<BLUE<< i << " "<< arr[i]<<"\n";
    }

    //cout<<"-----------------------------"<<endl; 
    cout<<BLUE<<"Global maximum position :" <<thisMaximum.pos<<endl;
    return thisMaximum;
}

double FindBaselineLevel(double *arr, int baseline_region_end)
{
    double baseline_sum = 0;
    double baseline_sqr_sum = 0;
    double baseline_rms, baseline_level;
    for(int i = 0; i < baseline_region_end; i++)
    {
        double y = arr[i];
        baseline_sum +=y;
        baseline_sqr_sum+=y*y;
    }
          
    baseline_level = baseline_sum/baseline_region_end;
    double baseline_rms_sqr = baseline_sqr_sum/baseline_region_end - baseline_level*baseline_level;
    baseline_rms = baseline_rms_sqr > 0 ? TMath::Sqrt(baseline_rms_sqr):0;
    //cout<<BLUE<<"baseline rms = "<< baseline_rms<<endlr;

    return baseline_rms;
}

STARTPOINT FindStartPoint(double *arr, double *arrt, double baseline_rms, int start, int pos, float dt)
{
    STARTPOINT start_point;
    start_point.ampl = 0;
    start_point.pos = 1; 
    start_point.x = -1; 
    for (int i = pos; i >=start ; --i)
    {
        //if(arr[i] - baseline_rms > 0)
        if(arr[i] > baseline_rms)
        {
            start_point.ampl = arr[i];
            cout<< RED << start_point.ampl<<endlr;
            start_point.pos = i;
            //cout<<BLUE<< start_point.pos<<endlr; 
            start_point.x = arrt[i];
            //cout<<GREEN<< start_point.x<<endlr;
            break;
        }
    }


    //cout<<"Starting point position: \t" <<start_point.x<<"\t in the point of: \t"<< start_point.pos<<endl;
    return start_point;
}

ENDPOINT FindEndPoint(double *arr, double *arrt, double baseline_rms, int start, float dt, int arrpoints, GLOBALMAXIMUM global_maximum)
{
    ENDPOINT end_point;
    //GLOBALMAXIMUM global_maximum;
    end_point.ampl = 0;
    end_point.pos = arrpoints-1; 
    
    //int j = type > 0 ? start : global_maximum.pos;  //returns a value based on the condition

    for (int i = global_maximum.pos; i < arrpoints; ++i)
    {
        if(arr[i] > baseline_rms)
        {
            end_point.ampl = arr[i];
            cout<< RED << end_point.ampl<<endlr;
            end_point.pos = i;
            //cout<<BLUE<< end_point.pos<<endlr; 
            end_point.x = arrt[i];
            cout<<GREEN<< end_point.x<<endlr;
            break;
        }
    }


    //cout<<"Ending point position: \t" <<end_point.x<<"\tin the point of: \t"<< end_point.pos<<endl;
    return end_point;
}

SEARCHEPEAKENDPOINT SearchEpeakEndPoint(double *arr, double *arrt, GLOBALMAXIMUM global_maximum, ENDPOINT end_point, int arrpoints)
{
        SEARCHEPEAKENDPOINT e_peak_end;
        e_peak_end.ampl = 11111111.;
        e_peak_end.pos = global_maximum.pos; 

        int j = global_maximum.pos +10; 
        if (j>arrpoints-1)
            j = arrpoints-2;
        for (int i = global_maximum.pos + 1.5; i < j && i+1<arrpoints; ++i)
        {
            if (arr[i] < e_peak_end.ampl)
            {
                e_peak_end.ampl = arr[i];
                e_peak_end.pos = i;
                if(arr[i+1]>0)
                    break;
            }
        }
        e_peak_end.x = arrt[e_peak_end.pos];
        //cout<<"electron peak end point @ "<< e_peak_end.x <<endl;

    return e_peak_end;
}

double fermi_dirac_general(double *x, double *par)
{
    double fdreturn = par[0]/TMath::Power((1.+TMath::Exp(-(x[0]-par[1])*par[2])),par[3]);
    //double fdreturn = par[0]/TMath::Power((1+TMath::Exp(-(x[0]-par[1])*par[2])),par[3]); 
    //double fdreturn = par[0]/(1+par[3]*TMath::Exp(-(x[0]-par[1])*par[2]));
    return fdreturn;
}

double fermi_dirac(double *x, double *par)
{
    //double fdreturn = par[0]/(par[1]+TMath::Exp(-par[2]*(x[0]-par[3])))+par[4];
    double fdreturn = par[0]/(1.+TMath::Exp(-par[2]*(x[0]-par[1])))+par[3];
    return fdreturn;
}

double fermi_dirac_generalsub(double *x, double *par)
{
    //double fdreturn = (par[0]/TMath::Power((1+TMath::Exp(-(x[0]-par[1])*par[2])),par[3])) - (par[0]/TMath::Power((1+TMath::Exp(-(x[0]-par[4])*par[5])),par[6]));
    double fdreturn = (par[0]/(1.+TMath::Exp(-par[2]*(x[0]-par[1])))+par[3]) + (par[0]/TMath::Power((1.+TMath::Exp(-(x[0]-par[4])*par[5])),par[6]));
    //double fdreturn = par[0]/(1+TMath::Exp(-par[2]*(x[0]-par[1])))+par[3] - par[0]/(1+TMath::Exp(-par[4]*(x[0]-par[5])))+par[6];
    //double fdreturn = (par[0]/(1+par[3]*TMath::Exp(-(x[0]-par[1])*par[2]))) - (par[0]/(1+par[6]*TMath::Exp(-(x[0]-par[4])*par[5])));
    return fdreturn;
}

double fermi_dirac_sym_1(double *x, double *par)
{
    double fdreturn = (1./(1.+TMath::Exp(-par[2]*(x[0]-par[1]))))+0.*par[0];
    //double fdreturn = par[0]/(1+TMath::Exp((x[0]-par[1])*par[2]))+par[3];
    return fdreturn;
}
double fermi_dirac_sym_double(double *x, double *par)
{
    // double fdreturn = (par[0]/(1.+TMath::Exp(-par[2]*(x[0]-par[1])))+par[3])*(1./(1.+TMath::Exp(-par[5]*(x[0]-par[4]))));
    double fdreturn = (par[0]/(1.+TMath::Exp(-par[2]*(x[0]-par[1]))))*(1./(1.+TMath::Exp(-par[5]*(x[0]-par[4]))))+par[3];
    //double fdreturn = par[0]/(1+TMath::Exp((x[0]-par[1])*par[2]))+par[3];
    return fdreturn;
}

double fermi_dirac_sym_double_charge(double x, double *par)
{
   double fdreturn = (par[0]/(1.+TMath::Exp(-par[1]*(x-par[2])))+par[3])*(1./(1.+TMath::Exp(-par[4]*(x-par[5]))));
   //double fdreturn = par[0]/(1+TMath::Exp((x[0]-par[1])*par[2]))+par[3];
   return fdreturn;
}

double slope_at_x(double *arr, double dt, PEAKPARAM *par)
{
    int i1 = par->stime_pos;
    int i2 = par->maxtime_pos;
    for (int i =  i1; i > i2; --i)
    {
        if (arr[i] < 0.2*par->ampl)
        {
          i1 = i;
          i2 = (i+1);
          break;
        }
    }
    double x1 = i1*dt; 
    double x2 = i2*dt;
    double y1 = arr[i1];
    double y2 = arr[i2];

    double slope = (y2-y1)/(x2-x1); 
    //pass = i1;
    return slope;
}

double Xpoint_linear_interpolation(double *arr, double dt, PEAKPARAM *par)
{   
    double slope = slope_at_x(arr, dt, par ) ;
    double x1 = par->stime_pos*dt;
    double y1 = arr[par->stime_pos];

    double x = (0.2*par->ampl-y1)/slope + x1;
    //return fabs(x);  
    return x;  
}


bool TimeSigmoid(int maxpoints, double *arr, double dt, PEAKPARAM *par, int evNo, double sig_shift, int tshift)
{
        //cout<<MAGENTA<<"Preparing for Sigmoind Fit" <<endlr;
      ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
      minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
      minimizer->SetMaxIterations(10000);  // for GSL 
      minimizer->SetTolerance(1e-3);
      minimizer->SetPrintLevel(1);  
        //dt=arrt[1]-arrt[0];
        par->sig_start_pos = (int) (par->stime_pos - 1./dt ); // minus 1 ns
        while(par->sig_start_pos<tshift) par->sig_start_pos+=1;
   
        // par->sig_end_pos = (int) (par->maxtime_pos + sig_shift/dt ); // add 2 ns (this was adding too much!!)
        par->sig_end_pos = par->maxtime_pos;
        while(par->sig_end_pos>=maxpoints-50) par->sig_end_pos-=1;

        int Npoints = par->sig_end_pos - par->sig_start_pos+1; 
        if(Npoints >100 || Npoints<=1)
        {
          cout<<BLUE<<"Start time pos = "<<par->sig_start_pos *dt <<"  sig_shift = "<< sig_shift *dt <<BLUE<<"  endpoint = "<< par->sig_end_pos *dt <<endl;
          cout<<GREEN<<"Start time pos = "<<par->sig_start_pos <<"  sig_shift = "<< sig_shift  <<BLUE<<"  endpoint = "<< par->sig_end_pos  <<endl;
          cout<<MAGENTA<<"Attention : "<<" in Event " << evNo <<" Sigmoid fit has few points ==> Number of points on sig_waveform ==>"<< Npoints <<endlr;
          cout<<RED<<"END point sigmoid = "<<par->sig_end_pos<<endlr;
         //           return (kFALSE);
        }

       double x[1000], y[1000], erx[1000], ery[1000];
       int points =0;
#ifdef DEBUGMSG
  cout << "par->rms: " << par->rms << endl;
#endif

       for (int i = 0; i <= par->sig_end_pos - par->sig_start_pos; i++)
        {
            y[i] = arr[par->sig_start_pos + i];
            //x[i] = arrt[i+sig_start_pos];
            x[i] = (par->sig_start_pos + i)*dt;
            erx[i] = 0;
            ery[i] = par->rms;
            points++;
        }
        //cout<<GREEN<<"points for the fit = "<<points<<endlr;
       

        TGraphErrors* sig_waveform = new TGraphErrors(points, x, y, erx, ery);

        double fit_start_point = x[0];
        double fit_end_point = x[par->sig_end_pos - par->sig_start_pos];
        
        TF1 *sig_fit =  new TF1("sig_fit",fermi_dirac,fit_start_point,fit_end_point, 4);
        
        double sig_pars[4];
     
        double y_half_point = 0.5*par->ampl;
    
        // double  x_mid_left = (x[par->sig_end_pos - par->sig_start_pos]+ x[1])/2.;
        double  x_mid_left = (fit_end_point + x[1])/2.;

        // double steepness_left = 5.0/(x[par->sig_end_pos - par->sig_start_pos] - x[1]);
        double steepness_left = 5.0/(fit_end_point - x[1]);

        sig_pars[0] = arr[par->maxtime_pos]; // - arr[sig_start_pos];
        sig_pars[1] = x_mid_left; 
        sig_pars[2] = steepness_left;
        sig_pars[3] = 0.0;//arr[sig_start_pos];

        sig_fit->SetParameters(sig_pars[0],sig_pars[1],sig_pars[2],sig_pars[3]);

      //        double range = par->maxtime_pos-(int)(0.05/dt);
      double range = (0.5/dt);
      sig_fit->SetParError(0, 0.01*abs(sig_fit->GetParameter(0)));
      sig_fit->SetParError(1, 0.05*range);
      sig_fit->SetParError(2, 0.05);
      sig_fit->SetParError(3, 0.001*abs(sig_fit->GetParameter(0)));
#ifdef DEBUGMSG
      for (int i = 0; i < 4; i++) {
        cout << "Parameter " << i << ": " << sig_fit->GetParameter(i)
             << " (error: " << sig_fit->GetParError(i) << ")" << endlr;
      }
#endif
//Save the fit results and get the success of the fit
      TFitResultPtr r_single = sig_waveform->Fit("sig_fit", "QMR0S");
#ifdef DEBUGMSG
      cout << MAGENTA << "Final parameters for sig_fit MM:" << endlr;
      for (int i = 0; i < 4; i++) {
        cout << "Parameter " << i << ": " << sig_fit->GetParameter(i)
             << " (error: " << sig_fit->GetParError(i) << ")" << endlr;
      }
#endif       
        //sig_fit->GetParameters(&sig_pars[0]);

        for (int i=0;i<4;i++)
        {
          //par->sigmoidR[i]=sig_pars[i];
          par->sigmoidR[i] = sig_fit->GetParameter(i);

          //cout<<"SIGMOID PARAMETERS = "<< i <<" = "<<par->sigmoidR[i]<<endl;
        }
      bool SigmoidfitSuccess = isSigmoidfitSuccessful(r_single);

      if(!SigmoidfitSuccess) {
        cout<<RED<<"Attention : "<<" in Event " << evNo <<" Sigmoid fit has failed ==> Number of points on sig_waveform ==>"<< Npoints <<endlr;
        // Open a file to write the failed event number
        ofstream failedEventsFile("failed_events_sigmoid.txt", ios::app);
        if (failedEventsFile.is_open()) {
          failedEventsFile << "Event " << evNo << " failed to Sigmoid fit." << endl;
          failedEventsFile.close();
        } else {
          cerr << "Unable to open file to write failed event." << endl;
        }
      }

      par->tfit20 =  sig_pars[1] - (1./sig_pars[2])*(TMath::Log(sig_pars[0]/((0.2*par->ampl-sig_pars[3])-1.)));
      //cout<<RED<<"sigmoid timepoint ="<< par->tfit20<<endlr;
      par->chi2_sigmoid = sig_fit->GetChisquare();

#ifdef DEBUGMSG
  cout << "chi2 sigmoid " << par->chi2_sigmoid << endl;


  sig_waveform->SetMarkerSize(1.0); // Adjust the size to your preference
  sig_waveform->SetMarkerStyle(20); // Use a specific marker style
  TCanvas *c1 = new TCanvas("c-timesigmoid", "Fit Result Time Sigmoid", 2000, 1400);
  // sig_waveformd->SetMarkerColor(kRed); // Optional: Set a marker color
  sig_waveform->Draw("AP");
  sig_fit->Draw("same");

  // Plot vertical lines for the range of the fit at x[0] and x[par->sig_end_pos]
  TLine *l1 = new TLine(fit_start_point, 0, fit_start_point, par->ampl);
  TLine *l2 = new TLine(fit_end_point, 0, fit_end_point, par->ampl);
  l1->SetLineColor(kRed);
  l2->SetLineColor(kRed);
  l1->Draw("same");
  l2->Draw("same");

  // Create a TPaveText to manually display the fit parameters
  TPaveText *stats = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC"); // NDC: normalized coordinates
  stats->SetFillColor(0);  // Transparent background
  stats->SetTextAlign(12); // Align left
  stats->SetBorderSize(1);

  // Add fit parameters manually (adjust based on your parameters)
  stats->AddText("Fit Parameters (sig_fit):");
  stats->AddText(Form("Param 0: %.3f #pm %.3f", sig_fit->GetParameter(0), sig_fit->GetParError(0)));
  stats->AddText(Form("Param 1: %.3f #pm %.3f", sig_fit->GetParameter(1), sig_fit->GetParError(1)));
  stats->AddText(Form("Param 2: %.3f #pm %.3f", sig_fit->GetParameter(2), sig_fit->GetParError(2)));
  stats->AddText(Form("Param 3: %.3f #pm %.3f", sig_fit->GetParameter(3), sig_fit->GetParError(3)));
  stats->AddText(Form("Chi2: %f", sig_fit->GetChisquare()));
  stats->AddText(Form("NDF: %d", sig_fit->GetNDF()));
  stats->AddText(Form("Chi2/NDF: %f", sig_fit->GetChisquare() / sig_fit->GetNDF()));
  stats->AddText(Form("Fit Range: %.3f - %.3f", fit_start_point, fit_end_point));

  stats->Draw();

  c1->Update();
  // cin.get();
  // c1->SaveAs("fit_result_single.png");
#endif

       return SigmoidfitSuccess;


} 
bool TimeSigmoidMCP(int maxpoints, double *arr, double dt, PEAKPARAM *par, int evNo, double sig_shift, int tshift)
{
      ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
      minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
      minimizer->SetMaxIterations(10000);  // for GSL 
      minimizer->SetTolerance(1e-3);
      minimizer->SetPrintLevel(1);  
        //cout<<MAGENTA<<"Preparing for Sigmoind Fit" <<endlr;
        //dt=arrt[1]-arrt[0];
        // par->sig_start_pos = (int) (par->stime_pos - 0.1/dt ); // minus 1 ns
        par->sig_start_pos = par->stime_pos;
        while(par->sig_start_pos<tshift) par->sig_start_pos+=1;
   
        // par->sig_end_pos = (int) (par->maxtime_pos + sig_shift/dt ); // add 2 ns
        par->sig_end_pos = par->maxtime_pos;
        while(par->sig_end_pos>=maxpoints-50) par->sig_end_pos-=1;

        int Npoints = par->sig_end_pos - par->sig_start_pos+1; 
        if(Npoints >100 || Npoints<=1) {cout<<BLUE<<"Start time pos = "<<par->sig_start_pos *dt <<"  sig_shift = "<< sig_shift *dt <<BLUE<<"  endpoint = "<< par->sig_end_pos *dt <<endl;
          cout<<GREEN<<"Start time pos = "<<par->sig_start_pos <<"  sig_shift = "<< sig_shift  <<BLUE<<"  endpoint = "<< par->sig_end_pos  <<endl;
          cout<<MAGENTA<<"Attention : "<<" in Event " << evNo <<" Sigmoid fit MCP has few points ==> Number of points on sig_waveform ==>"<< Npoints <<endlr;
          cout<<RED<<"END point sigmoid = "<<par->sig_end_pos<<endlr;
        }
          

       double x[1000], y[1000], erx[1000], ery[1000];
       int points =0;
   
       for (int i = 0; i <= par->sig_end_pos - par->sig_start_pos; i++)
        {
            y[i] = arr[par->sig_start_pos + i];
            //x[i] = arrt[i+sig_start_pos];
            x[i] = (par->sig_start_pos + i)*dt;
            erx[i] = 0;
            ery[i] = par->rms;
            points++;

        }
        //cout<<GREEN<<"points for the fit = "<<points<<endlr;
       

        TGraphErrors* sig_waveform = new TGraphErrors(points, x, y, erx, ery);
        double fit_start_point = x[0];
        double fit_end_point = x[par->sig_end_pos - par->sig_start_pos];

        TF1 *sig_fit =  new TF1("sig_fit",fermi_dirac,fit_start_point,fit_end_point, 4);
        
        double sig_pars[4];
     
        double y_half_point = 0.5*par->ampl;
    
        double  x_mid_left = (fit_end_point + x[1])/2.;

        double steepness_left = 5.0/(fit_end_point - x[1]);


        sig_pars[0] = arr[par->maxtime_pos]; // - arr[sig_start_pos];
        sig_pars[1] = x_mid_left; 
        sig_pars[2] = steepness_left;
        sig_pars[3] = 0.0;//arr[sig_start_pos];

        sig_fit->SetParameters(sig_pars[0],sig_pars[1],sig_pars[2],sig_pars[3]);

//        double range = par->maxtime_pos-(int)(0.05/dt);
        double range = (0.5/dt);
        sig_fit->SetParError(0, 0.01*abs(sig_fit->GetParameter(0)));
        sig_fit->SetParError(1, 0.05*range);
        sig_fit->SetParError(2, 0.05);
        sig_fit->SetParError(3, 0.001*abs(sig_fit->GetParameter(0)));
#ifdef DEBUGMSG
        for (int i = 0; i < 4; i++) {
          cout << "Parameter " << i << ": " << sig_fit->GetParameter(i)
               << " (error: " << sig_fit->GetParError(i) << ")" << endlr;
        }
#endif
        //Save the fit results and get the success of the fit
        TFitResultPtr r_single = sig_waveform->Fit("sig_fit", "QMR0S");
#ifdef DEBUGMSG
        cout << MAGENTA << "Final parameters for sig_fit MCP:" << endlr;
        for (int i = 0; i < 4; i++) {
          cout << "Parameter " << i << ": " << sig_fit->GetParameter(i)
               << " (error: " << sig_fit->GetParError(i) << ")" << endlr;
        }
#endif
        //sig_waveform->Fit("sig_fit", "QMR0S");
        //sig_waveform->Fit("sig_fit", "QMR0S");

        //sig_fit->GetParameters(&sig_pars[0]);

        for (int i=0;i<4;i++)
        {
          //par->sigmoidR[i]=sig_pars[i];
          par->sigmoidR[i] = sig_fit->GetParameter(i);
          //cout<<"SIGMOID PARAMETERS = "<< i <<" = "<<par->sigmoidR[i]<<endl;
        }
//debugging the fit results
//cin.get(); //press enter to continue
    bool SigmoidfitSuccess = isSigmoidfitSuccessful(r_single);
    if(!SigmoidfitSuccess) {
      cout<<RED<<"Attention : "<<" in Event " << evNo <<" Sigmoid fit has failed ==> Number of points on sig_waveform ==>"<< Npoints <<endlr;
      // Open a file to write the failed event number
      ofstream failedEventsFile("failed_events_sigmoid_MCP.txt", ios::app);
      if (failedEventsFile.is_open()) {
        failedEventsFile << "Event " << evNo << " failed to Sigmoid fit." << endl;
        failedEventsFile.close();
      } else {
        cerr << "Unable to open file to write failed event." << endl;
      }
    }
#ifdef DEBUGMSG
      TCanvas *c1 = new TCanvas("c1-a", "Fit Result MCP", 1000, 800);
      // sig_waveformd->SetMarkerColor(kRed); // Optional: Set a marker color
      sig_waveform->Draw("AP");
      sig_fit->Draw("same");

      // Create a TPaveText to manually display the fit parameters
      TPaveText *stats1 = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC"); // NDC: normalized coordinates
      stats1->SetFillColor(0);  // Transparent background
      stats1->SetTextAlign(12); // Align left
      stats1->SetBorderSize(1);

      // Add fit parameters manually (adjust based on your parameters)
      stats1->AddText("Fit Parameters (sig_fit_mcp):");
      stats1->AddText(Form("Param 0: %.3f #pm %.3f", sig_fit->GetParameter(0), sig_fit->GetParError(0)));
      stats1->AddText(Form("Param 1: %.3f #pm %.3f", sig_fit->GetParameter(1), sig_fit->GetParError(1)));
      stats1->AddText(Form("Param 2: %.3f #pm %.3f", sig_fit->GetParameter(2), sig_fit->GetParError(2)));
      stats1->AddText(Form("Param 3: %.3f #pm %.3f", sig_fit->GetParameter(3), sig_fit->GetParError(3)));
      stats1->AddText(Form("Chi2/NDF: %.3f", sig_fit->GetChisquare() / sig_fit->GetNDF()));

      stats1->Draw();

      c1->Update();
      c1->SaveAs("fit_result_single_mcp.png");

#endif
      par->tfit20 =  sig_pars[1] - (1./sig_pars[2])*(TMath::Log(sig_pars[0]/((0.2*par->ampl-sig_pars[3])-1.)));
      //cout<<RED<<"sigmoid timepoint ="<< par->tfit20<<endlr;
      par->chi2_sigmoid = sig_fit->GetChisquare();

  return SigmoidfitSuccess;


} 
void TimeSigmoidDraw(int maxpoints, double *arr, double *arrt, PEAKPARAM* par, int evNo, TCanvas *sig_canvas)
{
   
    double x[1000], y[1000], erx[1000], ery[1000];
    int points =0;

    for (int i = 0; i <= par->sig_end_pos - par->sig_start_pos; i++)
    {       
        y[i] = arr[par->sig_start_pos + i];
        //x[i] = arrt[i+sig_start_pos];
        x[i] = arrt[par->sig_start_pos + i];
        erx[i] = 0;
        ery[i] = par->rms;
        points++;

    }
    //cout<<BLUE<<"draw points for the fit = "<<points<<endlr;

    TGraphErrors* sig_waveform = new TGraphErrors(points, x, y, erx, ery);

    // Create a TGraph for the fitted sigmoid curve
    TF1 *sigmoidFit = new TF1("sigmoidFit", fermi_dirac, x[0], x[par->sig_end_pos - par->sig_start_pos], 4);
    sigmoidFit->SetParameters(par->sigmoidR[0], par->sigmoidR[1], par->sigmoidR[2], par->sigmoidR[3]);

    sigmoidFit->SetParName(0,"max_amplitude");
    sigmoidFit->SetParName(1,"x_mid_point");
    sigmoidFit->SetParName(2,"steepness_left");
    sigmoidFit->SetParName(3,"baseline");
    // Overlay the original data and the fitted sigmoid curve
    //TCanvas* sig_canvas = new TCanvas("sig_canvas", "Sigmoid Fit", 800, 600);
    sig_waveform->SetMarkerStyle(20); 
    sig_waveform->Draw("AP");
    sigmoidFit->Draw("L SAME");



    
    if(sig_waveform!=0)
    {
        TH1F *h1 = (TH1F*) sig_waveform->GetHistogram();

        h1->GetXaxis()->SetTitle("Time [ns]");
        h1->GetYaxis()->SetTitle("Amplitude [V]");
    }
    // Save the canvas to an image file
    //sig_canvas->SaveAs("sigmoid_fit_visualization.png");
    
    //cout<<"Sigmoid rise has been drawn for event "<<evNo<<endlr;


}
bool FullSigmoid(int maxpoints, double *arr, double dt, PEAKPARAM *par, int evNo, double sig_shift, int tshift)
{
      // ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
      // //minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      // minimizer->SetMaxFunctionCalls(10);
      // minimizer->SetMaxIterations(10000);  // for GSL
      // minimizer->SetTolerance(5e-4);
      // minimizer->SetPrintLevel(0);
      //

      ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
      ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
      ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
      ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
      ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);

       int points;
       double x[1000], y[1000], erx[1000], ery[1000];
        //cout<<" loop boundaries "<<par->sig_end_pos - par->sig_start_pos<<endl;
        for (int i = 0; i <= par->sig_end_pos - par->sig_start_pos; i++)
        {
            y[i] = arr[par->sig_start_pos + i];
            //x[i] = arrt[i+sig_start_pos];
            x[i] = (par->sig_start_pos + i)*dt;
            erx[i] = 0;
            ery[i] = par->rms;
            points++;

        }
        //cout<<GREEN<<"points for sigmoid fit = "<<points<<endlr;

        //int start = (int) (par->stime_pos - 6./dt);  // 6 ns
        int sigstart = par->sig_start_pos;
        int sigpeak = par->maxtime_pos + (int) (3./dt);
        //while(start < 0) start++;

        int sigend = par->maxtime_pos + (int)(10./dt);
        if(sigend<=sigstart) {
          sigend = sigstart + 100;
          cout<<" sig end is less than sig start "<<endl;
        }

        par->tot_sig_end_pos =  sigend;
   
        double y_half_point = 0.5*par->ampl;
        double x_mid_left = (x[par->sig_end_pos - par->sig_start_pos]+ x[1])/2.;
        double steepness_left = 5./(x[par->sig_end_pos-par->sig_start_pos] - x[1]);
        
        double x_d[1000], y_d[1000], erx_d[1000], ery_d[1000];

        int Npointsd = 0;
        int extention =  (int) (SIGMOID_EXTENTION / dt);
        //cout<<"loooooooop boundaries II "<< sigend+extention-sigstart<<endl; 
        for (int i = 0; i < sigend+extention-sigstart &&i<1000; ++i)
        {
            x_d[i] = (i+sigstart)*dt;
            y_d[i] = arr[i+sigstart];
            erx_d[i] = 0;
            ery_d[i] =par->rms; 
            if (i> sigend-sigstart) ery_d[i]*=10.;
            Npointsd++;
            
        }

        double sig_lim_min = sigstart*dt;
        double sig_lim_max = sigend*dt;
        double sig_pars_d[4];

      TF1 *sig_fitd =  new TF1("sig_fitd",fermi_dirac,sig_lim_min,sig_lim_max, 4);

      TGraphErrors* sig_waveformd = new TGraphErrors(Npointsd, x_d, y_d, erx_d, ery_d);


      for (int i=0;i<4; i++)
        sig_pars_d[i] = par->sigmoidR[i];


      sig_fitd->SetParameters(sig_pars_d);

      double x_mid_right =  par->maxtime_pos*dt + 3. ;  // ns

      /// double sigmoid fit here
      double sig_lim_min2 = par->maxtime_pos*dt;
      // double sig_lim_max2 = par->maxtime_pos*dt+(1.5); //ns
      par->tot_sig_end_pos = par->maxtime_pos + (int) (1.5/dt) + 1;
      double sig_lim_max2 = par->tot_sig_end_pos*dt;
      double sig_pars[4];

      if (x_mid_right <= sig_lim_min2 || x_mid_right >= sig_lim_max2)
         x_mid_right = (sig_lim_min2+sig_lim_max2)/2.;

      if (steepness_left <= 0 || steepness_left >= 5)
         steepness_left = (0.+5.)/2.;
         
      sig_pars[0] = arr[par->maxtime_pos];
      sig_pars[1] = x_mid_right;
      sig_pars[2] = -steepness_left;
      sig_pars[3] = arr[par->maxtime_pos]*0.1;
      //cout<<GREEN<<"Preparing 2nd Simple sigmoid FIT - normalization"<<endlr;
      //cout<<RED<<" Initial parameters sig_fit2 = "<< sig_pars[0] <<" "<< sig_pars[3]<<" From the sig_fitd "<< sig_fitd->GetParameter(0)<< " "<< sig_fitd->GetParameter(3)<<endlr;


  //Fit implementation with Root Default minimizer
        TF1 *sig_fit2 =  new TF1("sig_fit2",fermi_dirac,sig_lim_min2,sig_lim_max2, 4);
        sig_fit2->SetParameters(sig_pars);
        sig_fit2->FixParameter(0,sig_fitd->GetParameter(0));
        sig_fit2->FixParameter(3,sig_fitd->GetParameter(3));
        //cout<<GREEN<<"HERE IS THE fitiing FIX PARAMETERS SIFIT2!" << endlr;

        //the error encountered comes from the fact that all the parameters are fixed while they are not set to be fixed
        //Allow parameters 1 and 2 to vary
        // sig_fit2->SetParameter(1, sig_fit2->GetParameter(1));
        // sig_fit2->SetParameter(2, sig_fit2->GetParameter(2));

        // Set initial step sizes for the varying parameters
        sig_fit2->SetParError(1, 0.2 * abs(sig_fit2->GetParameter(1)));
        sig_fit2->SetParError(2, 0.2 * abs(sig_fit2->GetParameter(2)));

        //parameter limits were affecting the fit and producing the error
        sig_fit2->SetParLimits(1, sig_lim_min2, sig_lim_max2); // Assuming positive midpoint
        sig_fit2->SetParLimits(2, -5, 0);  // Assuming negative steepness

#ifdef DEBUGMSG
   //Debugging the fit parameters

          // Print initial parameter values and errors, and set step sizes
          cout << GREEN << "Initial parameters BEFORE Fitting for sig_fit2:" << endlr;
          for (int i = 0; i < 4; i++) {
            double initialError = sig_fit2->GetParError(i);
            double initialValue = sig_fit2->GetParameter(i);


            // Check if parameter is fixed
            Double_t parMin, parMax;
            sig_fit2->GetParLimits(i, parMin, parMax);
            //bool isFixed = (parMin == parMax);
            bool isFixed = (initialError == 0);

            cout << "Parameter " << i << ": " << initialValue
                 << " +/- " << initialError
                 << " (fixed: " << (isFixed ? "yes" : "no") << ")" << endlr;

            // Set step size if parameter is not fixed
            if (!isFixed) {
              double stepSize = (initialError > 0) ? initialError : 0.1 * std::abs(initialValue);
              sig_fit2->SetParError(i, stepSize);
              cout << "  Setting step size for parameter " << i << " to " << stepSize << endlr;
            }
          }

#endif
        // Perform the fit
        TFitResultPtr r = sig_waveformd->Fit("sig_fit2", "QMR0S");
        //Debugging the fit
        if (r->IsValid())
          {
#ifdef DEBUGMSG
          r->Print("V"); //prints the info of the fit
          TMatrixDSym cov = r->GetCovarianceMatrix(); //get covariance matrix
          cout << MAGENTA << "Fit successful!" << endlr;
#endif
          }

        else {
          cout << RED << "Fit failed sig_fit2." << endlr;
        }

        // Print final parameter values
#ifdef DEBUGMSG
        cout << GREEN << "Final parameters for sig_fit2:" << endlr;
        for (int i = 0; i < 4; i++) {
          cout << "Parameter " << i << ": " << sig_fit2->GetParameter(i)
               << " (error: " << sig_fit2->GetParError(i) << ")" << endlr;
        }

      cout << "Fit Result:" << endl;
      cout << "Chi-square: " << r->Chi2() << endl;
      cout << "NDF: " << r->Ndf() << endl;
      cout << "Probability: " << r->Prob() << endl;
      cout << "Status: " << r->Status() << endl;
      //sig_waveformd->Fit("sig_fit2", "QMR0S");
      //sig_waveformd->Fit("sig_fit2", "QMR0S");
      cout<<MAGENTA<<"HERE IS THE fitiing SIFIT2!" << endlr;

 #endif
      sig_fit2->SetRange(sig_lim_min-3.,sig_lim_max2+10.);
      sig_fit2->SetLineColor(kCyan);


      for (int i=0;i<4; i++)
        par->sigmoidF[i] = sig_fit2->GetParameter(i);
#ifdef DEBUGMSG
 //Debugging the fit
      // TCanvas *c = new TCanvas("c", "Fit Result", 800, 600);
      // sig_waveformd->Draw("AP");
      // sig_fit2->Draw("same");
      // c->SaveAs("fit_result.png");
      //
#endif

      //cout<<RED<<"HERE IS THE ERROR SIFIT2!" << endlr;


  // TF1 *sig_fit1 =  new TF1("sig_fit1",fermi_dirac_sym_1,sig_lim_min,sig_lim_max2+10., 3);
  //     sig_fit1->SetParameter(1,sig_fit2->GetParameter(1));
  //     sig_fit1->SetParameter(2,sig_fit2->GetParameter(2));
  //
  // // Print final parameter values
  // cout << BLUE << "Final parameters for sig_fit1:" << endlr;
  // for (int i = 0; i < 4; i++) {
  //   cout << "Parameter " << i << ": " << sig_fit1->GetParameter(i)
  //        << " (error: " << sig_fit1->GetParError(i) << ")" << endlr;
  // }

  TF1 *sig_fittot =  new TF1("sig_fittot",fermi_dirac_sym_double,sig_lim_min,sig_lim_max2+SIGMOID_EXTENTION, 6);
      for (int i=0;i<4;i++)
        sig_fittot->SetParameter(i,sig_fitd->GetParameter(i));

      sig_fittot->SetParameter(3+1,sig_fit2->GetParameter(1));
      sig_fittot->SetParameter(3+2,sig_fit2->GetParameter(2));

#ifdef DEBUGMSG
  cout<<" sig_lim_min = "<<sig_lim_min<<" sig_lim_max2 = "<<" Sigmoid extension "<<SIGMOID_EXTENTION<<endl;
#endif
  double range = sig_lim_max2 + SIGMOID_EXTENTION - sig_lim_min;
    sig_fittot->SetParError(0, 0.1*abs(sig_fit2->GetParameter(0)));
    sig_fittot->SetParError(1, 0.2*range);
    sig_fittot->SetParError(2, 0.02);
    sig_fittot->SetParError(3, 0.01*abs(sig_fitd->GetParameter(3)));
    sig_fittot->SetParError(4, 0.01*range);
    sig_fittot->SetParError(5, 0.03);

  // Set parameter limits
  sig_fittot->SetParLimits(1, sig_lim_min, sig_lim_max2+SIGMOID_EXTENTION); // Assuming positive midpoint
  sig_fittot->SetParLimits(2, -5, 5);  // Assuming negative steepness
  sig_fittot->SetParLimits(4, sig_lim_min, sig_lim_max2+SIGMOID_EXTENTION); // Assuming positive midpoint
  sig_fittot->SetParLimits(5, -5, 5);  // Assuming negative steepness



  bool any_zero = false;
  for (int i = 0; i < 6; i++) {
    if (sig_fittot->GetParError(i) == 0) {
      any_zero = true;
      break;
    }
  }

  if (any_zero) {
    cout << RED << "At least one parameter has an error of 0." << endlr;
    cin.get(); //press enter to continue
  }


  // Print final parameter values
#ifdef DEBUGMSG
  cout << MAGENTA << "Final parameters for sig_fittot:" << endlr;
  for (int i = 0; i < 6; i++) {
    cout << "Parameter " << i << ": " << sig_fittot->GetParameter(i)
         << " (error: " << sig_fittot->GetParError(i) << ")" << endlr;
  }
#endif

  double fit_double_start_point = sig_lim_min;
  double fit_double_end_point = sig_lim_max2+2.6;
  gErrorIgnoreLevel = kError;
  TFitResultPtr r_tot = sig_waveformd->Fit("sig_fittot", "QMR0S");
  gErrorIgnoreLevel = kInfo;
  // TFitResultPtr r_tot = sig_waveformd->Fit("sig_fittot", "MR0S");
  sig_fittot->SetRange(fit_double_start_point,fit_double_end_point);
  sig_fittot->SetLineColor(kRed);


      for (int i=0;i<6;i++)
        par->sigmoidtot[i]=sig_fittot->GetParameter(i);

      if (r_tot->IsValid())
      {
    #ifdef DEBUGMSG
        r_tot->Print("V"); //prints the info of the fit
        TMatrixDSym cov = r_tot->GetCovarianceMatrix(); //get covariance matrix
    #endif
        cout << MAGENTA << "Fit successful!" << endlr;
      }

      else {
        cout << RED << "Fit failed sig_fittot." << endlr;
        //cin.get(); //press enter to continue

        // Open a file to write the failed event number
        ofstream failedEventsFile("failed_events_double_Sigmoid.txt", ios::app);

        if (failedEventsFile.is_open()) {
          failedEventsFile << "Event " << evNo << " failed to fit." << endl;
          failedEventsFile.close();
        } else {
          cerr << "Unable to open file to write failed event." << endl;
        }
      }


  //Debugging the fit
#ifdef DEBUGMSG
  sig_waveformd->SetMarkerSize(1.0); // Adjust the size to your preference
  sig_waveformd->SetMarkerStyle(20); // Use a specific marker style

  TCanvas *c1 = new TCanvas("c1-a", "Fit Result", 1000, 800);
  // sig_waveformd->SetMarkerColor(kRed); // Optional: Set a marker color
  sig_waveformd->Draw("AP");
  sig_fit2->Draw("same");

  // Create a TPaveText to manually display the fit parameters
  TPaveText *stats1 = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC"); // NDC: normalized coordinates
  stats1->SetFillColor(0);  // Transparent background
  stats1->SetTextAlign(12); // Align left
  stats1->SetBorderSize(1);

  // Add fit parameters manually (adjust based on your parameters)
  stats1->AddText("Fit Parameters (sig_fit2):");
  stats1->AddText(Form("Param 0: %.3f #pm %.3f", sig_fit2->GetParameter(0), sig_fit2->GetParError(0)));
  stats1->AddText(Form("Param 1: %.3f #pm %.3f", sig_fit2->GetParameter(1), sig_fit2->GetParError(1)));
  stats1->AddText(Form("Param 2: %.3f #pm %.3f", sig_fit2->GetParameter(2), sig_fit2->GetParError(2)));
  stats1->AddText(Form("Param 3: %.3f #pm %.3f", sig_fit2->GetParameter(3), sig_fit2->GetParError(3)));
  stats1->AddText(Form("Chi2/NDF: %.3f", sig_fit2->GetChisquare() / sig_fit2->GetNDF()));

  stats1->Draw();

  c1->Update();
  c1->SaveAs("fit_result_single.png");
  c1->SaveAs("fit_result_single.eps");

  TCanvas *c2 = new TCanvas("c2-a", "Fit Result", 1000, 800);
  sig_waveformd->Draw("AP");
  sig_fitd->Draw("same");

  TPaveText *stats2 = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC"); // NDC: normalized coordinates
  stats2->SetFillColor(0);  // Transparent background
  stats2->SetTextAlign(12); // Align left
  stats2->SetBorderSize(1);

  stats2->AddText("Fit Parameters (sig_fitd):");
  stats2->AddText(Form("Param 0: %.3f #pm %.3f", sig_fitd->GetParameter(0), sig_fitd->GetParError(0)));
  stats2->AddText(Form("Param 1: %.3f #pm %.3f", sig_fitd->GetParameter(1), sig_fitd->GetParError(1)));
  stats2->AddText(Form("Param 2: %.3f #pm %.3f", sig_fitd->GetParameter(2), sig_fitd->GetParError(2)));
  stats2->AddText(Form("Param 3: %.3f #pm %.3f", sig_fitd->GetParameter(3), sig_fitd->GetParError(3)));
  //stats2->AddText(Form("Chi2/NDF: %.3f", sig_fitd->GetChisquare() / sig_fitd->GetNDF()));

  stats2->Draw();

  c2->Update();
  c2->SaveAs("fit_result_single_left.png");
  c2->SaveAs("fit_result_single_left.eps");

  TCanvas *c3 = new TCanvas("c3-a", "Fit Result", 1000, 800);
  sig_waveformd->Draw("AP");
  sig_fittot->Draw("same");

  TPaveText *stats3 = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC"); // NDC: normalized coordinates
  stats3->SetFillColor(0);  // Transparent background
  stats3->SetTextAlign(12); // Align left
  stats3->SetBorderSize(1);

  stats3->AddText("Fit Parameters (sig_fittot):");
  stats3->AddText(Form("Param 0: %.3f #pm %.3f", sig_fittot->GetParameter(0), sig_fittot->GetParError(0)));
  stats3->AddText(Form("Param 1: %.3f #pm %.3f", sig_fittot->GetParameter(1), sig_fittot->GetParError(1)));
  stats3->AddText(Form("Param 2: %.3f #pm %.3f", sig_fittot->GetParameter(2), sig_fittot->GetParError(2)));
  stats3->AddText(Form("Param 3: %.3f #pm %.3f", sig_fittot->GetParameter(3), sig_fittot->GetParError(3)));
  stats3->AddText(Form("Param 4: %.3f #pm %.3f", sig_fittot->GetParameter(4), sig_fittot->GetParError(4)));
  stats3->AddText(Form("Param 5: %.3f #pm %.3f", sig_fittot->GetParameter(5), sig_fittot->GetParError(5)));
  stats3->AddText(Form("Chi2/NDF: %.3f", sig_fittot->GetChisquare() / sig_fittot->GetNDF()));
  stats3->Draw();

  c3->Update();
  c3->SaveAs("fit_result_tot.png");
  c3->SaveAs("fit_result_tot.eps");

  // Print final parameter values
  cout << BLUE << "Final parameters for sig_fittot:" << endlr;
  for (int i = 0; i < 6; i++) {
    cout << "Parameter " << i << ": " << sig_fittot->GetParameter(i)
         << " (error: " << sig_fittot->GetParError(i) << ")" << endlr;
  }
#endif
  bool doubleSigmoidfitSuccess = isdoubleSigmoidfitSuccessful(r_tot);
  par->chi2_doubleSigmoid = sig_fittot->GetChisquare();

#ifdef DEBUGMSG
  if(doubleSigmoidfitSuccess) {
    cout<<GREEN<<"Fit successful!"<<endlr;
  }
  else {
    cout<<RED<<"Fit failed sig_fittotal (second check)."<<endlr;
  }
#endif
  double fit_integral = sig_fittot->Integral(sig_lim_min, (sig_lim_max2+SIGMOID_EXTENTION));
  //par->echargefit =  abs(fit_integral) / 50;
  par->echargefit =  fit_integral;
  par->e_peak_end_pos = fit_double_end_point;



  //cout<<YELLOW<<"Double Sigmoid Processed"<<endlr;
  //cout<<BLUE<<"Epeak charge fit = " << par->echargefit <<endlr;
  //cout<<YELLOW<<"Total charge = "<< par->totchargefixed<<endl;

  return doubleSigmoidfitSuccess;

}

//This is the original Draw Function
void FullSigmoidDraw(int maxpoints, double *arr, double *arrt, double dt, PEAKPARAM* par, int evNo, TCanvas *sigcanv)
{
    int start = par->sig_start_pos;  // 6 ns
    int sigpeak = par->maxtime_pos;

    int end = par->tot_sig_end_pos;
    double x_d[1000], y_d[1000], erx_d[1000], ery_d[1000];

    int Npointsd = 0;
    int extention =  (int) (SIGMOID_EXTENTION / dt);

    for (int i = 0; i < end+extention-start &&i<1000; ++i)
    {
        x_d[i] = (i+start)*dt;
        y_d[i] = arr[i+start];
        erx_d[i] = 0;
        ery_d[i] =par->rms;
        if (i> end-start) ery_d[i]*=10.;
        Npointsd++;

    }

    TGraphErrors* sig_waveformd = new TGraphErrors(Npointsd, x_d, y_d, erx_d, ery_d);

    double sig_lim_min = par->sig_start_pos * dt;
    double sig_lim_max2 = par->tot_sig_end_pos * dt;


    TF1 *sigmoidFitTOT =  new TF1("sigmoidFitTOT",fermi_dirac_sym_double,sig_lim_min,sig_lim_max2+2.6, 6);
    for (int i=0;i<6;i++)
        sigmoidFitTOT->SetParameter(i,par->sigmoidtot[i]);


    sig_waveformd->GetHistogram()->GetXaxis()->SetTitle("Time [ns]");
    sig_waveformd->GetHistogram()->GetYaxis()->SetTitle("Amplitude [V]");
    sig_waveformd->Draw("AP");

    sigmoidFitTOT->SetLineColor(kBlue);
    sigmoidFitTOT->Draw("L SAME");
    //sig_fittot->Draw("same");

  //cout<<RED<<"HERE IS THE ERROR SIGFITot from DRAW!" << endlr;


    TF1 *sig_fitd =  new TF1("sig_fitd",fermi_dirac,sig_lim_min,par->sig_end_pos*dt, 4);

      for (int i=0;i<4; i++)
        sig_fitd->SetParameter(i,par->sigmoidR[i]);

      sig_fitd->SetLineColor(kGreen+2);
      sig_fitd->Draw("L SAME");

      TF1 *sig_fit2 =  new TF1("sig_fit2",fermi_dirac,par->maxtime_pos *dt-3,par->maxtime_pos *dt +1.5 +10, 4);
      for (int i=0;i<4; i++)
        sig_fit2->SetParameter(i,par->sigmoidF[i]);

      sig_fit2->SetLineColor(kCyan+2);
      sig_fit2->Draw("L SAME");

      TF1 *sig_fit1 =  new TF1("sig_fit1",fermi_dirac_sym_1,par->maxtime_pos *dt-3,par->maxtime_pos *dt +1.5 +10, 3);
      sig_fit1->SetParameter(0, 0);
      sig_fit1->SetParameter(1,sig_fit2->GetParameter(1));
      sig_fit1->SetParameter(2,sig_fit2->GetParameter(2));
      sig_fit1->SetLineColor(kYellow+2);
      sig_fit1->Draw("L SAME");

}


bool isdoubleSigmoidfitSuccessful(TFitResultPtr r) {
  if (r->IsValid()) {
    return true;
  } else {
    return false;
  }
}

bool isSigmoidfitSuccessful(TFitResultPtr r) {
  if (r->IsValid()) {
    return true;
  } else {
    return false;
  }
}


int AnalyseLongPulseCiv(int points,int evNo, double* data, double dt, double* drv, PEAKPARAM *par, double threshold, double sig_shift, int tshift)
{
  /// use the integrated+filtered pulse to define a region where a trigger occured. (integral above threshold) 
  ///pulses are considered negative!!!
#ifdef DEBUGMSG
  cout<<MAGENTA<<"Starting Analysis for cividec at start point " << tshift <<endlr;
#endif
  if (points - tshift < 50) return -1;
#ifdef DEBUGMSG
  cout << "Threshold = " << threshold << endlr;
  cin.get();
#endif
  int ntrig=0;
  int tpoint=0;

  double drv_start_trig = -0.0025;  /// emperical
//   drv_start_trig = -threshold;
  double drv_end_trig = 0.00002;
  double drv_second_pulse_fraction_trigger = 0.2;  // Look for second pulse, end first pulse if derivative is
                                                   // less than this fraction of the first pulse
  //double drv_end_trig = 0.002;
  if (threshold <0.0025)
    drv_end_trig = 0.000005;

  par->tot[0]=0;

  for (int i=tshift; i<points; i++)   {
      //cout << i << "  data[i] = " << data[i] << " threshold = " << threshold << " drv[i] = " << drv[i] << " drv_start_trig = " << drv_start_trig << endlr;
      if (i > 2500) break;
      if (data[i]<=threshold && drv[i]<drv_start_trig) {
      tpoint = i;
      ntrig=1;
//       par->tot[0]=1;
      break;
    }
    else
      tpoint=i;
  }

#ifdef DEBUGMSG
  cout << RED << "Trigger point = " << tpoint << endlr;
  cout << RED << "value = "<<data[tpoint]<<endlr;
#endif
  if (ntrig<=0) return (-1); // cout<<"No trigger in event!"<<endl;

//    cout<<"tpoint = "<<tpoint*dt<<endl;
  if (tpoint>=points-10) return (-1);


  double miny = data[tpoint];
  double mindy = drv[tpoint];
  bool secondary_pulse = false;
  par->charge=0.;
  par->maxtime_pos=tpoint;
  par->ampl=data[tpoint];

  for (int i=tpoint; i<points; i++)
  {
    // cout << "i = " << i << " data[i] = " << data[i] << " miny = " << miny << endlr;
    if (data[i]<miny)
    {
      par->ampl=data[i];
      par->maxtime_pos=i;
      miny=data[i];
    }
    if(drv[i]<mindy) {
      mindy = drv[i];
    }
    if (data[i]<=threshold || drv[i]<drv_start_trig)
    {
      // cout << "data[i] = " << data[i] << " threshold = " << threshold << " drv[i] = " << drv[i] << " drv_start_trig = " << drv_start_trig << endlr;
      par->tot[0]++;
      par->ftime_pos=i; /// this is added to avoid a pulse at the end of the data that does not return to 0!!!
      // cout<<RED<<"Secondary pulse detected in event "<<evNo<<" derivative start point"<<par->ftime_pos<<endlr;

    }
    else //if (data[i]>threshold)
    {
      // cout << "data[i] = " << data[i] << " threshold = " << threshold << " drv[i] = " << drv[i] << " drv_start_trig = " << drv_start_trig << endlr;
      par->ftime_pos=i;
      par->tot[0]--;
      break;
    }
    /// note down the point the signal has gone above the threshold
  }
#ifdef DEBUGMSG
  cout << BLUE << "Maxtime position = " << par->maxtime_pos << endlr;
  cout << BLUE << "End of the pulse = " << par->ftime_pos << endlr; //correct end of the pulse
  cout << BLUE << "TOT (points) = "<< par->tot[0] <<" = " <<par->tot[0]*dt<<" ns"<<endlr;
#endif
  /// fast scan for risetime, risecharge and t_start
  par->t90=tpoint;
  par->t10=tpoint;
  par->stime_pos=tpoint;
  par->ttrig=tpoint;
  par->risecharge=0.;

  //find the 90% time
  for (int i=par->maxtime_pos; i>0; i--)
  {
     if (data[i]>=par->ampl*0.9)
     {
       par->t90=i;
       break;
     }
  }
//find the 10% time
  for (int i=par->t90; i>0; i--)
  {
    par->risecharge+=data[i];
    if (data[i]>=par->ampl*0.1 || (data[i]>0.5*threshold && fabs(drv[i])<=drv_end_trig ))
    {
      par->t10=i;
      break;
    }
  }
  //find the 10% time at the falling edge
  for (int i=par->maxtime_pos; i<points; i++)
  {

    if (data[i]>=par->ampl*0.1 || (data[i]>0.5*threshold  && fabs(drv[i])<=drv_end_trig ))
    {
      par->tb10=i;
      break;
    }
  }
//find start point
  for (int i=(int) par->t10; i>0; i--)
  {
#ifdef DEBUGMSG
    cout<<MAGENTA<<i<<"  "<<i*dt<<"  "<<data[i]<<"  "<<drv[i]<<"   "<<drv[i]*drv[tpoint]<<"   "<<threshold<<endlr;
#endif
    par->stime_pos=i;
//     if (i<tpoint)
//     {
//       par->charge+=data[i];
//     }
    if (data[i]>threshold/5. || (fabs(drv[i])<=drv_end_trig || drv[i]*drv[tpoint]<0 ) ) // && data[i]>threshold*0.8 ) )
    {
      par->stime_pos=i;
      break;
    }
  }

//   cout << GREEN << "End of the pulse = " << par->ftime_pos << endlr; //correct end of the pulse
//   cout<< MAGENTA <<"pulse duration at that point = "<<(par->ftime_pos-par->stime_pos)*dt<<endlr;
// cin.get();
  //find the point of derivative array that correspond to the par->ftime_pos
//   cout<<"derivative point value at the ftime_pos = "<<drv[par->ftime_pos]<<endl;


  for (int i=par->ftime_pos; i<points; i++) {
    par->ftime_pos=i;

    if (data[i]>threshold/5. || (fabs(drv[i])<=drv_end_trig && data[i]>threshold*0.8))
    {
      //cout<<BLUE<<"end of the pulse position"<<par->ftime_pos<<endlr;

      //cout<<RED<<drv[i]<<"  "<<threshold/5.<<" data "<<data[i]<<" "<<threshold*0.8<<endlr;
      // cin.get();
      break;
    }

    if (drv[i] < mindy * drv_second_pulse_fraction_trigger)
    {
#ifdef DEBUGMSG
      cout<<RED<<"Secondary pulse detected in event "<<evNo<<" derivative start point"<<par->ftime_pos<<endlr;
#endif
      secondary_pulse = true;
      break;
    }
  }



    // if ( (i - tpoint) * dt > start_second_pulse_check_time)
    // { // Check for secondary pulse via derivative
    //   // cout<<BLUE<<"Check for derivative "<<par->ftime_pos << " " << drv[i] << " " << mindy <<endlr;
    //   if (drv[i] < mindy * 0.1)
    //   {
    //     // cout<<RED<<"Secondary pulse detected in event "<<evNo<<" derivative start point"<<par->ftime_pos<<endlr;
    //     par->tot[0]--;
    //     secondary_pulse = true;
    //     break;
    //   }
    // }


// // extend to CIVIDEC_PULSE_DURATION ns in order to get the ion tail
  // for (int i=par->ftime_pos; i<points && i<par->stime_pos + CIVIDEC_PULSE_DURATION/dt; i++)
  // {
  //   par->ftime_pos=i;
  // }

  // cout<<"End of the pulse extended to = "<<par->ftime_pos<<endl;
  // cin.get();

  par->charge=0.;
  par->echarge=0.;
  par->ioncharge=0.;
  for (int i=par->stime_pos;i<=par->ftime_pos;i++)
    par->charge+=data[i];

  /// instead of using the risecharge, we assigne to this variable the "CIVIDEC charge", i.e. an integral of fixed duration. If signals are good, it should be charge = risecharge. 
  par->risecharge=0.;
  for (int i=par->stime_pos;i<points && i<par->stime_pos + CIVIDEC_PEAK_DURATION/dt ;i++)
    par->risecharge+=data[i];

  par->e_peak_end_ampl = 11111111.;
  par->e_peak_end_pos = par->maxtime_pos; 

  int j = par->maxtime_pos +10; 
  if (j>points-1) //array points
    j = points-2; //array points
  for (int i = par->maxtime_pos + 1.5; i < j && i+1<tshift; ++i) //array points
  {
    if (data[i] < par->e_peak_end_ampl)
    {
        par->e_peak_end_ampl = data[i];
        par->e_peak_end_pos  = i;
            if(data[i+1]>0)
                break;
    }
   }
//Calculate the charge from the start to the end of the electron peak
  for (int i = par->stime_pos ; i<par->e_peak_end_pos; i++)
  {
    par->charge +=data[i];
  }

  #ifdef DEBUGMSG
   cout<<YELLOW<<"Epeak charge before double sigmoid "<< par->charge<<" at e_peak_end_point ="<< par->e_peak_end_pos*dt<<endlr;
   cin.get();
 //calculate the integral from the start point to the end point of the waveform on a constant window of 120ns
     cout<<"CIVIDEC pulse duration in points = "<<CIVIDEC_PULSE_DURATION/dt<<" or in ns = "<<CIVIDEC_PULSE_DURATION<<endl;
     cout<<"CIVIDEC epeak pulse duration in points = "<<CIVIDEC_PEAK_DURATION/dt<<" or in ns = "<<CIVIDEC_PEAK_DURATION<<endl;
     cout<<"trigger point "<<par->stime_pos<<endl;
     cin.get();
#endif
      par->totchargefixed = 0;
      double tot_charge_fixed_position = 0;
      tot_charge_fixed_position = par->stime_pos+CIVIDEC_PULSE_DURATION/dt;
      for (int i = par->stime_pos; i < tot_charge_fixed_position; i++) {
        //integrate the waveform
        par->totchargefixed += data[i];
      }
      // cout<<RED<< "Total Charge FIXED on 120ns window: " << par->totchargefixed << endlr;
      // cout<<BLUE<<"Total Charge between start-end position: "<<par->charge<<endlr;
      // cin.get();

    /// make the sig fit for sigmoind timepoint.
#ifdef DEBUGMSG
     cout<<RED<<"Starting Sigmoid fit "<<endlr;
#endif
     bool SigmoidfitSuccess = TimeSigmoid(points, data, dt,par, evNo, sig_shift, tshift);
     par->SigmoidfitSuccess = SigmoidfitSuccess;
#ifdef DEBUGMSG
     cout<<GREEN<<"Starting Sigmoid interpolation "<<endlr;
#endif
       par->tnaive20 = Xpoint_linear_interpolation(data, dt, par);
       //cout<<BLUE<<"NAIVE TIME = "<< par->tnaive20<<endlr;
       //cout<<"Sigmoid Timepoint = "<<par->tfit20<<endl;
       //double Xpoint_linear_interpolation(double *arr, double dt, PEAKPARAM *par )
       //cout<<MAGENTA<<"NaiveTime, tfit20 =  "<< par->tfit20 <<endlr;
#ifdef DEBUGMSG
     cout<<BLUE<<"Starting Double Sigmoid "<<endlr;
     cout<<RED<<"FIT TIME = "<< par->tfit20<<endlr;
#endif

    // Epeak charge calculation and Fit Success check
      bool doubleSigmoidfitSuccess =  FullSigmoid(points, data, dt, par, evNo, sig_shift, tshift);
      par->doubleSigmoidfitSuccess = doubleSigmoidfitSuccess;
//Calculate the epeak charge in a fixed window of 6ns after the peak
      par->echargefixed = 0;
      double epeak_charge_fixed_position = 0;
      epeak_charge_fixed_position = par->stime_pos+CIVIDEC_PEAK_DURATION/dt;
      for (int i = par->stime_pos; i < epeak_charge_fixed_position; i++) {
        //integrate the waveform
        par->echargefixed += data[i];
      }
//Calculate the charge of the ion tail
      for (int i = par->e_peak_end_pos/dt; i<par->ftime_pos; i++)
      {
        par->ioncharge +=data[i];
      }
      #ifdef DEBUGMSG
       cout<<MAGENTA<<" Ion charge AFTER double sigmoid "<< par->ioncharge<<" at e_peak_end_point ="<< par->e_peak_end_pos/dt<<endlr;
       cin.get();
       #endif
      // cout<<RED<< "Epeak Charge FIXED on 6ns window: " << par->totchargefixed << endlr;
      // cout<<BLUE<<"Epeak Charge fit: "<<par->echargefit<<endlr;
      // cin.get();

      par->te_peak_end = par->e_peak_end_pos * dt;
      par->risecharge *= dt;
      par->tot[0] *= dt;
      par->maxtime = par->maxtime_pos*dt;
      par->t90 *= dt;
      par->t10 *= dt;
      par->tb10 *= dt;
      //cout<<"electron peak end point @ "<< e_peak_end.x <<endl;
      par->sampl = data[par->stime_pos];
      par->fampl = data[par->ftime_pos];
      par->bslch = -0.5 * (data[par->stime_pos] + data[par->ftime_pos])*(par->ftime_pos - par->stime_pos +1.)*dt;
      par->width = (par->ftime_pos-par->stime_pos)*dt;
      par->ampl*=-1.;
      par->charge*=-1.*dt;   ///charge is calculated in V * ns.
      par->totchargefixed*=-1.*dt; ///charge is calculated in V * ns.
      par->ioncharge*=-1.*dt; ///charge is calculated in V * ns.
      par->risecharge*=-1.*dt;
      par->risetime = (par->t90-par->t10)*dt;
      par->echargefit *= -1.*dt;  ///calculated from the integral of the fit V * ns.
      par->echargefixed *= -1.*dt; ///charge is calculated in V * ns.

  //cout<<YELLOW<<"First quick scan of parameters finished "<<endlr;
    
//   cout<<"tstart = "<<par->stime_pos*dt<<endl;
//   cout<<"t10 = "<<par->t10*dt<<endl;
//   cout<<"t90 = "<<par->t90*dt<<endl;
//   cout<<"tmax = "<<par->maxtime_pos*dt<<endl;
//   cout<<"tb10 = "<<par->tb10*dt<<endl;
//   cout<<"tend = "<<par->ftime_pos*dt<<endl;
//   
//   cout<<"rt = "<<(par->t90-par->t10)*dt<<endl;
//   cout<<"tot = "<<(par->tot)*dt<<endl;
//   cout<<"DT = "<<(par->tb10-par->t10)*dt<<endl;
//   cout<<"DTall = "<<(par->ftime_pos-par->stime_pos)*dt<<endl;
//   
//   cout<<"ampl = "<<par->ampl<<endl;
//   cout<<"charge = "<<par->charge*dt/N_INTEGRATION_POINTS<<endl;
//   cout<<"risecharge = "<<par->risecharge*dt/N_INTEGRATION_POINTS<<endl;
  return (par->ftime_pos);
}

int AnalyseLongPulseMCP(int points,int evNo, double* data, double dt, double* drv, PEAKPARAM *par, double threshold, double sig_shift, int tshift)
{
  /// use the integrated+filtered pulse to define a region where a trigger occured. (integral above threshold) 
  ///pulses are considered negative!!!
#ifdef DEBUGMSG
  cout<<MAGENTA<<"Starting Analysis for cividec MCP "<<endlr;
#endif
  if (points - tshift < 50) return -1;
  
  int ntrig=0;
  int tpoint=0;
  
  // double drvtrig = 0.00002;
  // if (threshold <0.0025)
  //   drvtrig = 0.000005;


  double drv_start_trig = -0.05;
  drv_start_trig = -threshold;
  double drv_end_trig = 0.00002;
  double drv_second_pulse_fraction_trigger = 0.2;  // Look for second pulse, end first pulse if derivative is less than this fraction of the first pulse
  //double drv_end_trig = 0.002;
  if (threshold <0.0025)
    drv_end_trig = 0.000005;

  par->tot[0]=0;
  
  for (int i=tshift; i<points; i++)   {
    if (data[i]<=threshold && drv[i]<drv_start_trig) {
      tpoint = i;
      ntrig=1;
      //       par->tot[0]=1;
      break;
    }
    else
      tpoint=i;
  }
  
#ifdef DEBUGMSG
  cout << RED << "Trigger point  MCP = " << tpoint << endlr;
#endif
  if (ntrig<=0) return (-1); // cout<<"No trigger in event!"<<endl;
  
//    cout<<"tpoint = "<<tpoint*dt<<endl;
  if (tpoint>=points-10) return (-1);

  double miny = data[tpoint];
  double mindy = drv[tpoint];
  bool secondary_pulse = false;
  par->charge=0.;
  par->maxtime_pos=tpoint;
  par->ampl=data[tpoint];


  for (int i=tpoint; i<points; i++)
  {
    if (data[i]<miny)
    {
      par->ampl=data[i];
      par->maxtime_pos=i;
      miny=data[i];
    }
    if(drv[i]<mindy) {
      mindy = drv[i];
    }
    if (data[i]<=threshold || drv[i]<drv_start_trig)
    {
      par->tot[0]++;
      par->ftime_pos=i; /// this is added to avoid a pulse at the end of the data that does not return to 0!!!
      // cout<<RED<<"Secondary pulse detected in event "<<evNo<<" derivative start point"<<par->ftime_pos<<endlr;

    }
    else //if (data[i]>threshold)
    {
      par->ftime_pos=i;
      par->tot[0]--;
      break;
    }
    /// note down the point the signal has gone above the threshold
  }

#ifdef DEBUGMSG
  cout << BLUE << "Maxtime position MCP = " << par->maxtime_pos << endlr;
  cout << BLUE << "End of the pulse MCP = " << par->ftime_pos << endlr; //correct end of the pulse
#endif
  /// fast scan for risetime, risecharge and t_start
  par->t90=tpoint;
  par->t10=tpoint;
  par->stime_pos=tpoint;
  par->ttrig=tpoint;
  par->risecharge=0.;

  for (int i=par->maxtime_pos; i>0; i--)
  {
     if (data[i]>=par->ampl*0.9)
     {
       par->t90=i;
       break;
     }
  }

  for (int i=par->t90; i>0; i--)
  {
    par->risecharge+=data[i];
    if (data[i]>=par->ampl*0.1 || (data[i]>0.5*threshold && fabs(drv[i])<=drv_end_trig ))
    {
      par->t10=i;
      break;
    }
  }
  //find the 10% time at the falling edge

  for (int i=par->maxtime_pos; i<points; i++)
  {

    if (data[i]>=par->ampl*0.1 || (data[i]>0.5*threshold && fabs(drv[i])<=drv_end_trig ))
    {
      par->tb10=i;
      break;
    }
  }
  
  //find start point
  for (int i=(int) par->t10; i>0; i--)
  {
    //     cout<<data[i]<<"  "<<fabs(drv[i])<<"   "<<threshold<<endl;
    par->stime_pos=i;
    //     if (i<tpoint)
    //     {
    //       par->charge+=data[i];
    //     }
    if (data[i]>threshold/5. || fabs(drv[i])<=drv_end_trig)//&& data[i]>threshold*0.8 ) )
    {
      par->stime_pos=i;
      break;
    }
  }

  // cout << GREEN << "End of the pulse MCP = " << par->ftime_pos << endlr; //correct end of the pulse
  // cout<< MAGENTA <<"pulse duration at that point MCP = "<<(par->ftime_pos-par->stime_pos)*dt<<endlr;
  // cin.get();

  for (int i=par->ftime_pos; i<points; i++) {
    par->ftime_pos=i;

    if (data[i]>threshold/5. || (fabs(drv[i])<=drv_end_trig && data[i]>threshold*0.8))
    {
      //cout<<BLUE<<"end of the pulse position"<<par->ftime_pos<<endlr;

      //cout<<RED<<drv[i]<<"  "<<threshold/5.<<" data "<<data[i]<<" "<<threshold*0.8<<endlr;
      // cin.get();
      break;
    }

    if (drv[i] < mindy * drv_second_pulse_fraction_trigger)
    {
      //cout<<RED<<"Secondary pulse detected in event "<<evNo<<" derivative start point"<<par->ftime_pos<<endlr;
      secondary_pulse = true;
      break;
    }
  }
//   cout<<"ftime_pos = "<<par->ftime_pos<<endl;
/// extend to CIVIDEC_PULSE_DURATION ns in order to get the ion tail
///
  // cout<<"End of the pulse extended to = "<<par->ftime_pos<<endl;
  // cin.get();
  // for (int i=par->ftime_pos; i<points && i<par->stime_pos + CIVIDEC_PULSE_DURATION/dt; i++)
  // {
  //   par->ftime_pos=i;
  // }
//   cout<<"extended to = "<<par->ftime_pos<<endl;
  
  par->charge=0.;
  for (int i=par->stime_pos;i<=par->ftime_pos;i++)
    par->charge+=data[i];

  /// instead of using the risecharge, we assigne to this variable the "CIVIDEC charge", i.e. an integral of fixed duration. If signals are good, it should be charge = risecharge. 
  par->risecharge=0.;
  for (int i=par->stime_pos;i<points && i<par->stime_pos + CIVIDEC_PEAK_DURATION/dt ;i++)
    par->risecharge+=data[i];

  par->e_peak_end_ampl = 11111111.;
  par->e_peak_end_pos = par->maxtime_pos; 

  int j = par->maxtime_pos +10; 
  if (j>points-1) //array points
    j = points-2; //array points
  for (int i = par->maxtime_pos + 1.5; i < j && i+1<tshift; ++i) //array points
  {
    if (data[i] < par->e_peak_end_ampl)
    {
        par->e_peak_end_ampl = data[i];
        par->e_peak_end_pos  = i;
            if(data[i+1]>0)
                break;
    }
   }
  
   par->te_peak_end = par->e_peak_end_pos * dt;
   par->risecharge *= dt;
   par->tot[0] *= dt;
   par->maxtime = par->maxtime_pos*dt;
   par->t90 *= dt;
   par->t10 *= dt;
   par->tb10 *= dt;

   par->sampl = data[par->stime_pos];
   par->fampl = data[par->ftime_pos];
   par->bslch = -0.5 * (data[par->stime_pos] + data[par->ftime_pos])*(par->ftime_pos - par->stime_pos +1.)*dt;
   par->width = (par->ftime_pos-par->stime_pos)*dt;
   par->ampl*=-1.;
   par->charge*=-1.*dt;   ///charge is calculated in V * ns. 
   par->risetime = (par->t90-par->t10)*dt;
  //cout<<YELLOW<<"First quick scan of parameters finished "<<endlr;

/// make the sig fit for sigmoind timepoint.
   // par->tfit20 = TimeSigmoidMCP(points, data, dt,par, evNo, sig_shift, tshift);

  bool SigmoidfitSuccess = TimeSigmoidMCP(points, data, dt,par, evNo, sig_shift, tshift);
  par->SigmoidfitSuccess = SigmoidfitSuccess;

  //cin.get();
   //cout<<BLUE<<"Time Sigmoid processed, tfit20 =  "<< par->tfit20 <<endlr;

   //cout<<RED<<"FIT TIME = "<< par->tfit20<<endlr;
   par->tnaive20 = Xpoint_linear_interpolation(data, dt, par);
   //cout<<BLUE<<"NAIVE TIME = "<< par->tnaive20<<endlr;
   //cout<<"Sigmoid Timepoint = "<<par->tfit20<<endl;
   //double Xpoint_linear_interpolation(double *arr, double dt, PEAKPARAM *par )
   //cout<<MAGENTA<<"NaiveTime, tfit20 =  "<< par->tfit20 <<endlr;

///make double sigmoid fit here for charge calculation 
    //FullSigmoid(points, data, dt, par, evNo, sig_shift, tshift);
    
    //cout<<YELLOW<<"FULL SIGMOID PROCESSED"<<endlr;


//   cout<<"tstart = "<<par->stime_pos*dt<<endl;
//   cout<<"t10 = "<<par->t10*dt<<endl;
//   cout<<"t90 = "<<par->t90*dt<<endl;
//   cout<<"tmax = "<<par->maxtime_pos*dt<<endl;
//   cout<<"tb10 = "<<par->tb10*dt<<endl;
//   cout<<"tend = "<<par->ftime_pos*dt<<endl;
//   
//   cout<<"rt = "<<(par->t90-par->t10)*dt<<endl;
//   cout<<"tot = "<<(par->tot)*dt<<endl;
//   cout<<"DT = "<<(par->tb10-par->t10)*dt<<endl;
//   cout<<"DTall = "<<(par->ftime_pos-par->stime_pos)*dt<<endl;
//   
//   cout<<"ampl = "<<par->ampl<<endl;
//   cout<<"charge = "<<par->charge*dt/N_INTEGRATION_POINTS<<endl;
//   cout<<"risecharge = "<<par->risecharge*dt/N_INTEGRATION_POINTS<<endl;
  return (par->ftime_pos);
}


void AddPar(PEAKPARAM* ipar, PEAKPARAM* spar) //The function copies the values from the ipar 
                                             //object to the spar object, with some modifications based on the value of dt.
{
   spar->maxtime_pos=ipar->maxtime_pos;
   spar->stime_pos=ipar->stime_pos;
   spar->ftime_pos=ipar->ftime_pos;
   spar->e_peak_end_pos=ipar->e_peak_end_pos;
   spar->sig_start_pos = ipar->sig_start_pos;
   spar->sig_end_pos = ipar->sig_end_pos;
   spar->tot_sig_end_pos = ipar->tot_sig_end_pos;
   
   spar->maxtime = ipar->maxtime; 
   spar->ampl=ipar->ampl;
   spar->e_peak_end_ampl = ipar->e_peak_end_ampl;
   spar->sampl=ipar->sampl;
   spar->fampl=ipar->fampl;
   //spar->t20=ipar->t20;
   //spar->st20=ipar->st20;
   spar->tfit20=ipar->tfit20;
   spar->tnaive20 = ipar->tnaive20;
   spar->te_peak_end = ipar->te_peak_end;
   
   
   //spar->sechargefixed=ipar->sechargefixed;
   //spar->secharge=ipar->secharge;
   //spar->echarge=ipar->echarge;
   spar->scharge=ipar->scharge;
   spar->charge=ipar->charge;
   spar->echargefixed=ipar->echargefixed;
   spar->echargefit=ipar->echargefit;
   spar->totchargefixed=ipar->totchargefixed;
   
   spar->risetime=ipar->risetime; ///10% - 90%
   spar->risecharge=ipar->risecharge;
   
   spar->width=ipar->width;
   
   spar->tot[0]=ipar->tot[0]; //the tot[0], charge, and risecharge values in the spar object are modified 
   
    for (int i=0;i<4;i++)
   {
     spar->sigmoidR[i]=ipar->sigmoidR[i];
     spar->sigmoidF[i]=ipar->sigmoidF[i];
   }
   for (int i=0;i<6;i++)
     spar->sigmoidtot[i]=ipar->sigmoidtot[i];

   
   spar->t10=ipar->t10;
   spar->tb10=ipar->tb10;
   spar->t90=ipar->t90;
   spar->ttrig=ipar->ttrig;
   
   spar->bslch=ipar->bslch; 
   spar->bsl=ipar->bsl;
   spar->rms=ipar->rms;

  spar->SigmoidfitSuccess = ipar->SigmoidfitSuccess;
   spar->doubleSigmoidfitSuccess = ipar->doubleSigmoidfitSuccess;
   //copy the values of bsl, rms, stime_pos, ftime_pos, and maxtime_pos 
   //from the ipar object to the corresponding variables in the spar object.
   
   //based on the value of dt

}
void AddPar(IPARAM* ipar, IPARAM* spar, double dt)
{
   spar->bsl=ipar->bsl;
   spar->rms=ipar->rms;
   spar->stime=ipar->stime;
   spar->ftime=ipar->ftime;
   spar->maxtime=ipar->maxtime;
   spar->tot=ipar->tot * dt;
   spar->ampl=ipar->ampl;
   spar->charge=ipar->charge * dt;
   spar->risecharge=ipar->risecharge * dt;
   spar->t10=ipar->t10 * dt;
   spar->tb10=ipar->tb10 * dt;
   spar->t90=ipar->t90 * dt;
   spar->ttrig=ipar->ttrig * dt;
   spar->width=ipar->width * dt;
   spar->sampl=ipar->sampl;
   spar->fampl=ipar->fampl;
   spar->bslch=ipar->bslch*dt; 
  
}

void DivideH(TH1D* h1, TH1D* h2)
{
   double x, y, sx;
   int n1 = h1->GetNbinsX();
   int n2 = h2->GetNbinsX();
//    if (n1 > n2) n1 = n2;
   
   for (int i=0;i<=n1; i++)
   {
     x = h1->GetBinContent(i);
     y = h1->GetBinWidth(i);  /// in seconds!
     sx = sqrt(x);
     if (y!=0)
     {
        x/=y;
	sx/=y;
     }
     h1->SetBinContent(i,x);
     h1->SetBinError(i,sx);
   }
}

void DivideHspark(TH1D* h1, TH1D* h2, TH1D* hspark)
{
   double x, y, sx, nsparks;
   double binwidth = h1->GetBinWidth(1);
   double npulses = floor(binwidth / 1.2) ;
   if (npulses==0) 
   {
     cout<<" ATTENTION !!!! Time binwidth smaller than LINAC4 periode !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
     npulses = 1;
   }
//    npulses++;
   cout<<"****  Rate evolution BinWidth = "<<binwidth<<endl;
   cout<<"****  Corresponding to "<<npulses<<" Linac4 pulses"<<endl;
   int n1 = h1->GetNbinsX();
   int n2 = h2->GetNbinsX();
   if (n1 != n2) 
   {
     cout<<"ATTENTION !!!!!!   Plot and normalization plot do not have the same number of bins!!!!!!!"<<endl;
     if (n1>n2) n1 = n2;
   }
   for (int i=0;i<=n1; i++)
   {
     x = h1->GetBinContent(i);
     y = h2->GetBinContent(i);
     nsparks = hspark->GetBinContent(i);
     
     if (y+nsparks > npulses) {
       cout<<" ******     found bin with "<<y<<" entries and " <<nsparks<<" sparks & recoveries    *******"<<endl;
       npulses = y;
     }
     y = npulses; 
     sx = sqrt(x);
     if (y!=0)
     {
        x/=y;
	sx/=y;
     }
     h1->SetBinContent(i,x);
     h1->SetBinError(i,sx);
   }
}

void ScaleHistoErr(TH1D* h1, double scale)
{
   double x, sx;
   double binwidth = h1->GetBinWidth(1);

   int n1 = h1->GetNbinsX();
   for (int i=0;i<=n1; i++)
   {
     x = h1->GetBinContent(i);
     sx = sqrt(x);
     if (scale!=0)
     {
        x*=scale;
	    sx*=scale;
     }
     h1->SetBinContent(i,x);
     h1->SetBinError(i,sx);
   }
}

int FilterHisto(TH1* h1, double scale=0.1)
{
   double x, sx;
   scale=fabs(scale);
   if (scale>1) return 0;
   sx = h1->GetMaximum();

   int n1 = h1->GetNbinsX();
   for (int i=0;i<=n1; i++)
   {
     x = h1->GetBinContent(i);
     if (x < scale*sx)
     {
       h1->SetBinContent(i,0.);
     }
   }
   return 1;
}


TCanvas *statsEditing() {
   // Create and plot a test histogram with stats
   TCanvas *se = new TCanvas;
   TH1F *h = new TH1F("h","test",100,-3,3);
   h->FillRandom("gaus",3000);
   gStyle->SetOptStat();
   h->Draw();
   se->Update();
   // Retrieve the stat box
   TPaveStats *ps = (TPaveStats*)se->GetPrimitive("stats");
   ps->SetName("mystats");
   TList *listOfLines = ps->GetListOfLines();
   // Remove the RMS line
   TText *tconst = ps->GetLineWith("RMS");
   listOfLines->Remove(tconst);
   // Add a new line in the stat box.
   // Note that "=" is a control character
   TLatex *myt = new TLatex(0,0,"Test = 10");
   myt ->SetTextFont(42);
   myt ->SetTextSize(0.04);
   myt ->SetTextColor(kRed);
   listOfLines->Add(myt);
   // the following line is needed to avoid that the automatic redrawing of stats
   h->SetStats(0);
   se->Modified();
   return se;
}



#endif
