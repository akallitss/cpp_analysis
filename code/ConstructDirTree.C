#include "MyFunctions.h"

int ConstructDirTree()
{
  
  char dirname[500];
  sprintf(dirname,"%s/",BASEDIRNAME);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/",DATADIRNAME);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/tmp/",WORKDIR);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/logs/",WORKDIR);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/logs/analysis",WORKDIR);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/",PLOTDIR);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/moreplots/",PLOTDIR);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s/",OUTDIRNAME);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  sprintf(dirname,"%s",PARAMDIRNAME);
  cout<<"Making directory: "<<dirname<<endl;
  gSystem->mkdir(dirname,kTRUE);
  
  return 1;

}
