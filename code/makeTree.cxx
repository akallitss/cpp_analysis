#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>

#include "MyFunctions.h"

using namespace std;

int main(int argc, char **argv)
{
  int poolid=0;
  int runid=0;
  if (argc>2)  {
    int tst=sscanf(argv[1],"%d", &runid);
    int tst2=sscanf(argv[2],"%d", &poolid);
    if (tst<=0||tst2<=0) {
      cout<<"syntax: makeTree <run number> <pool number>"<<endl; 
      return -1;
    }
  }
  else {
    cout<<"syntax: makeTree <run number> <pool number>"<<endl; 
    return -2;
  }
    
  if (runid<MINRUN || runid>MAXRUN){
    cout<<"run "<<runid<<" out of bounts ("<<MINRUN<<"-"<<MAXRUN<<")"<<endl;
    return -3;
  }
  char poolname[1000];
  if (poolid==0)
      sprintf(poolname,"GDD");
  else
      sprintf(poolname,"Pool%d",poolid);
  
  char fname[100];
  sprintf(fname,"%s/logs/logRun%03d-%s.txt",WORKDIR,runid,poolname);
  FILE *ftmp = fopen(fname,"w");
  fclose(ftmp);
  char command[2000];
  sprintf(command,"root -b -x '%s/MakeTreefromRawTreePicosec.C+(%d,%d)' -q >> %s &\n",CODEDIR,runid,poolid,fname);
  cout<<"executing:\n"<<command<<endl;
  system(command);

  return (1);
}
