#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

//#include <TLatex.h>
#include <TSystem.h>
#include <TH1D.h>

#include "getChannel.h"
#include "fillMAPS.h"
#include "correction.h"

using namespace std;

static map<int,int> m1;
static map<int,int>::iterator it_m1;

static map<int,double> m2;
static map<int,double> m3;
static map<int,double>::iterator it_m2;
static map<int,double>::iterator it_m3;

static map<string,int> m4;
static map<string,int> m5;
static map<string,int>::iterator it_m4;
static map<string,int>::iterator it_m5;

static map<string,int> m6;
static map<string,int> m7;
static map<string,int>::iterator it_m6;
static map<string,int>::iterator it_m7;

static map<int,double> m8;
static map<int,double> m9;
static map<int,double>::iterator it_m8;
static map<int,double>::iterator it_m9;



void getMaps(){
  getRunNumbers(&m1,&m2,&m3);
  getMapMAPMT(&m4,&m5);
  getMapMPPC(&m6,&m7);
  getTimeCalibrationDataMAPMT(&m8);
  getTimeCalibrationDataMPPC(&m9);
}

void upstreamMaroc(int fiber, THeader *run){
  bool answer=false;
  for(int i = 0; i < 8; i+=2){
    if(fiber == run->fiberRef[i]){
      answer=true;
      break;
    }
  }
  run->upstreamBoard=answer;
}

int getMarocBoard(int fiber, THeader *run){
  int board = -1;
  for(int i = 0; i < 8; i++){
    if(fiber == run->fiberRef[i]){
      board = run->marocBoard[i];
      break;
    }
  } 
  return board;
}

int getMarocChip(int mCh){
  int ret = -1;
  if(mCh < 64) ret = 0;
  if(mCh > 63 && mCh < 128) ret = 1;
  if(mCh > 127) ret = 2;
  return ret;
}


double timeCalibrationMAPMT(double time, int channel, int pmt){
  if(m8.empty())getMaps();
  int place = pmt*256+channel;
  double correction = m8.at(place);
  double newTime = time + correction;
  //if(channel<10)cout <<"Place: "<<pmt<<place<<" Correcion:  "<<time <<" + "<<correction <<" = " <<newTime <<endl;
  //if(channel<10)printf("Place (%2d, %3d)  %4d corr %7.2f %7.2f = %7.2f \n",pmt,channel,place,time,correction,newTime);
  //cin.get();
  return newTime;
}

double timeCalibrationMPPC(double time, int channel, int pmt){
  if(m9.empty())getMaps();
  int place = pmt*256+channel;
  double correction = m9.at(place);
  double newTime = time + correction;
  return newTime;
}


int getMAPMT_anode(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream){
  if(m4.empty()){
    //cout <<"[WARNING] filling again the static maps\n";
    getMaps();
  }
  string label = "N/A";
  int anode=-1;
  if(marocBoard==2 && chip==0) label=Form("IN1_%02d",mCh);
  else if(marocBoard==3 && chip==1) label=Form("IN1_%02d",mCh-64);
  else if(chip==2) label=Form("IN3_%02d",mCh-128);
  else if(marocBoard==3 && chip==0) anode=-1; 
  else {
    cout <<Form("[ERROR] Wrong chip number, marocBoard 2, chip %d\n",chip);
    exit(EXIT_FAILURE);
  }
  if(marocBoard!=3 || chip!=0){
    if(marocUpstream) anode = m4.at(label);
    else  anode = m5.at(label);
    if(anode==-1){
      cout <<"[ERROR] Something wrong linking setup file and maps. Check it!\n";
      exit(EXIT_FAILURE);
    }
  }
  return anode;
}

int getMPPC_ch(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream){
	if(m6.empty()){
		//cout <<"[WARNING] filling again the static maps\n";
		getMaps();
	}
	string label = "N/A";
	int nCh=-1;
	if(chip==0) label=Form("IN1_%02d",mCh);
	else if(chip==1) label=Form("IN1_%02d",mCh-64);
	else if(chip==2) label=Form("IN3_%02d",mCh-128);
	else {
		cout <<Form("[ERROR] Wrong chip number, marocBoard 2, chip %d\n",chip);
		exit(EXIT_FAILURE);
	}
	if(marocUpstream) nCh=m6.at(label);
	else nCh=m7.at(label);
	if(nCh==-1){
		cout <<"[ERROR] Something wrong linking setup file and maps. Check it!\n";
		exit(EXIT_FAILURE);
	}
	return nCh; 
}
/*
   void onlineTimeCalibration(THeader *run){
   TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",runHead->suite.c_str(),runHead->runNum);
   TTree *t = (TTree*) fdRICH->Get("data");
   uint evt, nedge;
   int slot[MAXDATA], fiber[MAXDATA], chM[MAXDATA], pol[MAXDATA], time[MAXDATA];
   double trigtime;
   t->SetBranchAddress("evt",&evt);
   t->SetBranchAddress("nedge",&nedge);
   t->SetBranchAddress("fiber",&fiber);
   t->SetBranchAddress("ch",&chM);
   t->SetBranchAddress("pol",&pol);
   t->SetBranchAddress("time",&time);

   FILE *fout;
   fout = fopen("../map/time_calib_MAPMT.txt","w");

   TH2D *hCalVsCh= new TH2D("hCalVsCh","#Delta time calibrated Vs Channel; ch[#]; #Delta_t [ns]",1026,-0.5,1025.5, tstaBin,staMin,staMax);
   TH2D *hNotVsCh= new TH2D("hNotVsCh","#Delta time not calibrated Vs Channel; ch[#]; #Delta_t [ns]",1026,-0.5,1025.5, tstaBin,staMin,staMax);

   vector<vector<TH1D*>> vCh;
   vector<vector<TF1*>> vCh;
   for(int i = 0; i < 8; i++){
   vector<TH1D*> tmpV;
   vector<TF1*> tmpF;
   for(int j = 0; j < 192; j++){
   TH1D *htmp = new TH1D(Form("h_5_%d_%d",i+4,j),Form("Time 5_%d_%d",i+4,j),100,0,100);
   tmpV.push_back(htmp);
   TF1 *ftmp = new TF1(Form("f_5_%d_%d",i+4,j),"gaus(0)",0,100);
   tmpF.push_back(ftmp);
   }
   vCh.push_back(tmpV);
   vF.push_back(tmpF);
   }
   for(int i=0;i<nentries;i++){
   data->GetEntry(i);
   int upEdge=0;
   for(int j = 0; j < nedge; j++){
   if(pol[j]==0)upEdge++;
   hitFlag[j]=false;
   }
   for(int j = 0; j < nedge; j++){
   for(int k = 0; k < nedge; k++){
   if(j == k)continue;
   if((time[j] < time[k] && pol[j]==pol[k]) || pol[j] < pol[k] ){
   int tmpTime=time[j];
   int tmpPol=pol[j];
   int tmpCh=ch[j];
   int tmpSlot=slot[j];
   int tmpFiber=fiber[j];
   time[j]=time[k];
   pol[j]=pol[k];
   ch[j]=ch[k];
   slot[j]=slot[k];
   fiber[j]=fiber[k];
   time[k]=tmpTime;
   pol[k]=tmpPol;
   ch[k]=tmpCh;
   slot[k]=tmpSlot;
   fiber[k]=tmpFiber;
   }
   }
   }
   double t0=0;
   bool trigger=false;
   int ntrig=0;
   for(int j = 0; j < nedge; j++){
   if(pol[j]!=0)continue;
   for(int k = 0; k < nedge; k++){
   if(pol[k]!=0)continue;
   if(fiber[j]!=fiber[k])continue;
if(ch[j]!=ch[k])continue;
if(time[k] < time[j] && pol[j]==pol[k])continue;
if(time[k] < time[j])continue;
double start_corr = timeWalkCorrection(time[j],time[k]);
if(fiber[j]==11 && ch[j]==27){
	trigger=true;
	ntrig++;
	t0=start_corr;
	hitFlag[j]=false;
}else if(fiber[j]!=11 || ch[j]>63){
	hitFlag[j]=false;
}else{
	hitFlag[j]=false;
}
}
}
if(trigger==false || ntrig!=1)continue;
for(int j = 0; j < nedge; j++){
	if(hitFlag[j]==false)continue;
	if(fiber[j] == 11 && ch[j] < 64) continue;
	vCh[fiber[j]-4][ch[j]]->Fill(time[j]-t0);
}
}//Loop on tree.
double mean=0;
int cont=0;
TCanvas *cFit = new TCanvas("cFit","cFit",1600,900);
for(int i = 0; i < vCh.size(); i++){
	for(int j = 0; j < vCh[i].size(); j++){
		if(vCh[i][j]->GetEntries()==0)continue;
		//vCh[i][j]->GetXaxis()->SetRangeUser(40,52);
		double p0 = vCh[i][j]->GetBinContent(vCh[i][j]->GetMaximumBin());
		double p1 = vCh[i][j]->GetBinCenter(vCh[i][j]->GetMaximumBin());
		vF[i][j]->SetParameters(p0,p1,3);
		vCh[i][j]->Fit(vF[i][j],"Q","",p1-3,p1+3);
		hCorr->Fill(vF[i][j]->GetParameter(1));
		mean+=vF[i][j]->GetParameter(1);
		cont++;
	}
}
cFit->Close();
mean/=cont;
double corr[1024];
double meanCorr=0;
int contCorr=0;
for(int i = 0; i < 1024; i++)corr[i]=0;
for(int i = 0; i < vCh.size(); i++){
	for(int j = 0; j < vCh[i].size(); j++){
		if(vCh[i][j]->GetEntries()==0) continue;
		vCh[i][j]->Draw();
		c2->Update();
		c2->Print(out_pdf.c_str());
		if(abs(vF[i][j]->GetParameter(1)-mean) > 10) cout <<Form("Check fiber %d ch %d\n",i+4,j);
		int k=reference(i+4,j);
		corr[k]=-vF[i][j]->GetParameter(1);
		meanCorr-=vF[i][j]->GetParameter(1);
		contCorr++;
		vCh[i][j]->Write();
		vF[i][j]->Write();
	}
}
meanCorr/=contCorr;
for(int i = 0; i < 1024; i++){
	if(corr[i]!=0) fprintf(fout,"%d %lf\n",i+1,corr[i]-meanCorr);
	else  fprintf(fout,"%d %lf\n",i+1,meanCorr);
}
fclose(fout);
for(int i=0;i<nentries;i++){
	data->GetEntry(i);
	for(int j = 0; j < nedge; j++){
		if(hitFlag[j]==true && pol[j]==0){
			hCalVsCh->Fill(reference(fiber[j],ch[j]),time[j]+(corr[reference(fiber[j],ch[j])]-meanCorr));
			hNotVsCh->Fill(reference(fiber[j],ch[j]),time[j]);
		}
	}
}
TCanvas *c7 = new TCanvas("c7","Correction canvas",1600,900);
c7->Divide(2);
c7->Draw();
c7->cd(1);
gPad->SetGrid(1,1);
auto tp1=hCalVsCh->ProfileX();
tp1->GetYaxis()->SetRangeUser(490,530);
tp1->Draw();
c7->cd(2);
gPad->SetGrid(1,1);
auto tp2=hNotVsCh->ProfileX();
tp2->GetYaxis()->SetRangeUser(490,530);
tp2->Draw();
c7->Update();
c7->Print(out_pdf.c_str());
//c7->Print("c7.png");
c7->Print("onlineCalibration.root");
c7->Close();
}*/
