#define MAXDATA 10000
#include <iostream>
#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TError.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TSystem.h>

#include "definition.h"
#include "utility.h"
#include "computing.h"

static TH1D *hTime;
static vector <TH1D*> hHitCoinc;
static vector <TH1D*> hHitStart;
static vector <TH1D*> hHitEnd;
static vector <TH1D*> hHitDur;
static vector <TH2D*> hHitCorrDur;
static vector <TH2D*> hNotCalibVsCh;
static vector <TH2D*> hCalibVsCh;

static vector <TH1D*> vGasRadPMT;
static vector <TH1D*> vAeroRadPMT;

static TH2D *hHitCorr;
static TH2D *hMap;
static vector<TH2D*> vMap;
static TH2D *hMapNC;
static TH2D *hnMap;
static vector<TH2D*> vPiMap;
static vector<TH2D*> vKMap;
static vector<TH2D*> vPrMap;
static TH1D *hRadius;
static TH1D *hbRadius;
static vector<TH1D*> vnRadius;
static vector<TH1D*> vPiRadius;
static vector<TH1D*> vKRadius;
static vector<TH1D*> vPrRadius;
static vector<TH1D*> vCutPiRadius;
static vector<TH1D*> vCutKRadius;
static vector<TH1D*> vCutPrRadius;

static TH1D *hbHitStart;
static TH1D *hbHitEnd;
static TH2D *hbRadVsStart;
static TH2D *hbRadVsEnd;

static TH1D *hPiRad;
static TH1D *hKRad;
static TH1D *hPrRad;

static TH1D *htrig;
static TH1D *hx474;
static TH1D *hx519;
static TH1D *hx537;
static TH1D *hBeamCh;

static TH1D *hEdge;
static TH1D *hHitReco;

static TH2D *hUpGEM;
static TH2D *hDnGEM;
static TH2D *hBeam;
static TH2D *hBeamTheta;

static vector<TH1D*> vrsdRadius;
static vector<TH1D*> vrsdTime;

static vector<TH1D*> vspRadius;
static vector<TH1D*> vspTime;
static vector<TH1D*> vspPhoton;
static vector<TH1D*> vspnRadius;
static vector<TH1D*> vspnTime;
static vector<TH1D*> vspnPhoton;
static vector<TH1D*> vcutRadius;
static vector<TH1D*> vcutTime;
static vector<TH1D*> vcutPhoton;

static vector<TH1D*> vspSigPhoGas;
static vector<TH1D*> vspSigPhoAero;
static vector<TH1D*> vspnSigPhoGas;
static vector<TH1D*> vspnSigPhoAero;
static vector<TH1D*> vcutSigPhoGas;
static vector<TH1D*> vcutSigPhoAero;
static TH1D *hspSigVsPhoGas;
static TH1D *hspnSigVsPhoGas;
static TH1D *hcutSigVsPhoGas;
static TH1D *hspSigVsPhoAero;
static TH1D *hspnSigVsPhoAero;
static TH1D *hcutSigVsPhoAero;

static TH2D *hCalibVsTime;
static TH2D *hTimeVsCh;
static TH2D *hTimeWalkVsTime;
static TH2D *hStartVsDur;
static TH1D *hSingPhotRad;
static TH1D *hTimeRespectTrigger;

static vector<vector<TH1D*>> vNotTime;
static vector<vector<TH1D*>> vCalTime;

static int minPhoGas=2;
static int maxPhoGas=50;
static int minPhoAero=2;
static int maxPhoAero=15;

const char *labelCh[3]={"Pi","K","Pr"};

//---------------------------------------------------
void fillHistoMon(THeader *run){
  //---------------------------------------------------

  gErrorIgnoreLevel=kWarning;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int evt, nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10], cutPhoton[10], fiber[MAXDATA], ch[MAXDATA], time[MAXDATA], otime[MAXDATA];
  double y[MAXDATA], x[MAXDATA], r[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nt[MAXDATA], nttw[MAXDATA], dur[MAXDATA], rsdRadius[MAXDATA], rsdTime[MAXDATA], spRadius[10], spTime[10], spnRadius[10], spnTime[10], cutRadius[10], cutTime[10];
  float gx0, gy0, gx1, gy1, gxa, gya, gxtheta, gytheta;
  bool trigSig[MAXDATA], goodHit[MAXDATA], goodPhoton[MAXDATA], coincPhoton[MAXDATA],externalPhoton[MAXDATA], innerPhoton[MAXDATA], outerPhoton[MAXDATA], goodSP[10], goodSPN[10], goodCUT[10];
  double trigtime, x474time, x519time, x537time;
  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("x474time",&x474time);
  t->SetBranchAddress("x519time",&x519time);
  t->SetBranchAddress("x537time",&x537time);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("otime",&otime);
  t->SetBranchAddress("nx",&nx);
  t->SetBranchAddress("ny",&ny);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("nttw",&nttw);
  t->SetBranchAddress("dur",&dur);
  t->SetBranchAddress("trigSig",&trigSig);
  t->SetBranchAddress("goodHit",&goodHit);
  t->SetBranchAddress("goodPhoton",&goodPhoton);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("externalPhoton",&externalPhoton);
  t->SetBranchAddress("innerPhoton",&innerPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  t->SetBranchAddress("rsdRadius",&rsdRadius);
  t->SetBranchAddress("rsdTime",&rsdTime);
  t->SetBranchAddress("gx0",&gx0);
  t->SetBranchAddress("gy0",&gy0);
  t->SetBranchAddress("gx1",&gx1);
  t->SetBranchAddress("gy1",&gy1);
  t->SetBranchAddress("gxa",&gxa);
  t->SetBranchAddress("gya",&gya);
  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  t->SetBranchAddress("spRadius",&spRadius);
  t->SetBranchAddress("spTime",&spTime);
  t->SetBranchAddress("spPhoton",&spPhoton);
  t->SetBranchAddress("goodSP",&goodSP);
  t->SetBranchAddress("spnRadius",&spnRadius);
  t->SetBranchAddress("spnTime",&spnTime);
  t->SetBranchAddress("spnPhoton",&spnPhoton);
  t->SetBranchAddress("goodSPN",&goodSPN);
  t->SetBranchAddress("cutRadius",&cutRadius);
  t->SetBranchAddress("cutTime",&cutTime);
  t->SetBranchAddress("cutPhoton",&cutPhoton);
  t->SetBranchAddress("goodCUT",&goodCUT);


  double n[3] = {1+run->px474*(nN2-1), 1+run->px519*(nCO2-1), 1+run->px537*(nCO2-1)};
  double thsPi[3] = {mPi/sqrt(n[0]*n[0]-1), mPi/sqrt(n[1]*n[1]-1), mPi/sqrt(n[2]*n[2]-1)};
  double thsK[3] = {mK/sqrt(n[0]*n[0]-1), mK/sqrt(n[1]*n[1]-1), mK/sqrt(n[2]*n[2]-1)};
  double thsPr[3] = {mPr/sqrt(n[0]*n[0]-1), mPr/sqrt(n[1]*n[1]-1), mPr/sqrt(n[2]*n[2]-1)};

  bool x474Pi=false, x474K=false, x474Pr=false;
  if(run->energyGeV > thsPi[0])x474Pi=true;
  if(run->energyGeV > thsK[0])x474K=true;
  if(run->energyGeV > thsPr[0])x474Pr=true;

  bool x519Pi=false, x519K=false, x519Pr=false;
  if(run->energyGeV > thsPi[1])x519Pi=true;
  if(run->energyGeV > thsK[1])x519K=true;
  if(run->energyGeV > thsPr[1])x519Pr=true;

  bool x537Pi=false, x537K=false, x537Pr=false;
  if(run->energyGeV > thsPi[2])x537Pi=true;
  if(run->energyGeV > thsK[2])x537K=true;
  if(run->energyGeV > thsPr[2])x537Pr=true;

  cout <<Form("Pi thresholds: %lf %lf %lf Status %d %d %d\n",thsPi[0],thsPi[1],thsPi[2],x474Pi,x519Pi,x537Pi);
  cout <<Form("K thresholds: %lf %lf %lf Status %d %d %d\n",thsK[0],thsK[1],thsK[2],x474K,x519K,x537K);
  cout <<Form("Pr thresholds: %lf %lf %lf Status %d %d %d\n",thsPr[0],thsPr[1],thsPr[2],x474Pr,x519Pr,x537Pr);
  //cin.get();

  cout <<Form("Filling monitoring histograms for run %d\n",run->runNum);
  for(int i = 0; i < t->GetEntries(); i++){
    if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    //GEM CUTS
    //printf("GEM %7.2f [%7.2f] %7.2f [%7.2f] %7.3f [%7.3f] \n",abs(gxa),GEM_CUT_X,abs(gya),GEM_CUT_Y,sqrt(gxtheta*gxtheta+gytheta*gytheta),GEM_CUT_R);
    if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
    if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;

    bool isPion   = false;
    bool isKaon   = false;
    bool isProton = false;
    if(run->beamChLogic==2){
      if(x474Pi == true && x519Pi==true && x537Pi==true && x474time!= 0 && x519time!=0 && x537time!=0){
        hBeamCh->Fill(0.,1);
        isPion = true;
      }
      if(x474K == false && x519K==false && x537K==true && x474time== 0 && x519time==0 && x537time!=0){
        hBeamCh->Fill(1.,1);
        isKaon = true;
      }
      if(x474Pr == false && x519Pr==false && x537Pr==false && x474time== 0 && x519time==0 && x537time==0){
        hBeamCh->Fill(2.,1);
        isProton = true;
      }
    }else if(run->beamChLogic==3){
      if(x474Pi == true && x519Pi==true && x537Pi==true && x474time!= 0 && x519time!=0 && x537time!=0){
        hBeamCh->Fill(0.,1);
        isPion = true;
      }
      if(x474K == false && x519K==true && x537K==true && x474time== 0 && x519time!=0 && x537time!=0){
        hBeamCh->Fill(1.,1);
        isKaon = true;
      }
      if(x474Pr == false && x519Pr==false && x537Pr==true && x474time== 0 && x519time==0 && x537time!=0){
        hBeamCh->Fill(2.,1);
        isProton = true;
      }
    }else cout <<"No beam cherenkov available for this run\n";
    if(isPion==true)printf("Is_pion!\n");
    if(isKaon==true)printf("Is_kaon!\n");
    if(isProton==true)printf("Is_proton!\n");

    hUpGEM->Fill(gx0,gy0);
    hDnGEM->Fill(gx1,gy1);
    hBeam->Fill(gxa,gya);
    hBeamTheta->Fill(gxtheta,gytheta);
    if(trigtime!=0)htrig->Fill(trigtime);
    if(x474time!=0)hx474->Fill(x474time);
    if(x519time!=0)hx519->Fill(x519time);
    if(x537time!=0)hx537->Fill(x537time);
    if(goodSP[4]==true && spPhoton[4] > minPhoGas && spPhoton[4]<maxPhoGas) vspSigPhoGas[spPhoton[4]]->Fill(spRadius[4]);
    if(goodSP[9]==true && spPhoton[9] > minPhoAero && spPhoton[9]<maxPhoAero) vspSigPhoAero[spPhoton[9]]->Fill(spRadius[9]);
    if(goodSPN[4]==true && spnPhoton[4] > minPhoGas && spnPhoton[4]<maxPhoGas) vspnSigPhoGas[spnPhoton[4]]->Fill(spnRadius[4]);
    if(goodSPN[9]==true && spnPhoton[9] > minPhoAero && spnPhoton[9]<maxPhoAero) vspnSigPhoAero[spnPhoton[9]]->Fill(spnRadius[9]);
    if(goodCUT[4]==true && cutPhoton[4] > minPhoGas && cutPhoton[4]<maxPhoGas) vcutSigPhoGas[cutPhoton[4]]->Fill(cutRadius[4]);
    if(goodCUT[9]==true && cutPhoton[9] > minPhoAero && cutPhoton[9]<maxPhoAero) vcutSigPhoAero[cutPhoton[9]]->Fill(cutRadius[9]); 

   /* if(goodCUT[9] == true && cutRadius[9]<145){
      printf(" Plotta eve %3d %4d \n",i,evt);
      for(int j = 0; j < nedge; j++){
        if(goodPhoton[j] && externalPhoton[j]) printf("%3d %3d %3d %3d (%4d, %4d) %7.2f %7.2f %lf \n",evt,j,externalPhoton[j],pol[j],fiber[j],ch[j],nttw[j],dur[j],nr[j]);
      }
      //cin.get();
    }*/
    for(int j = 0; j < 10; j++){
      if(goodSP[j]==true && spPhoton[j] > 2){
        vspRadius[j]->Fill(spRadius[j]);
        vspTime[j]->Fill(spTime[j]);
        vspPhoton[j]->Fill(spPhoton[j]);
      }
      if(goodSPN[j]==true && spnPhoton[j] > 2){
        vspnRadius[j]->Fill(spnRadius[j]);
        vspnTime[j]->Fill(spnTime[j]);
        vspnPhoton[j]->Fill(spnPhoton[j]);
      }
      if(goodCUT[j]==true && cutPhoton[j] > 2){
        vcutRadius[j]->Fill(cutRadius[j]);
        vcutTime[j]->Fill(cutTime[j]);
        vcutPhoton[j]->Fill(cutPhoton[j]);
      }
    }

    if(goodCUT[4]==true){
      if(isPion==true)vCutPiRadius[1]->Fill(cutRadius[4]);
      if(isPion==true)vCutPiRadius[0]->Fill(cutRadius[4]);
      if(isKaon==true)vCutKRadius[1]->Fill(cutRadius[4]);
      if(isKaon==true)vCutKRadius[0]->Fill(cutRadius[4]);
      if(isProton==true)vCutPrRadius[1]->Fill(cutRadius[4]);
      if(isProton==true)vCutPrRadius[0]->Fill(cutRadius[4]);
    }
    if(goodCUT[9]==true){
      if(isPion==true)vCutPiRadius[2]->Fill(cutRadius[9]);
      if(isPion==true)vCutPiRadius[0]->Fill(cutRadius[9]);
      if(isKaon==true)vCutKRadius[2]->Fill(cutRadius[9]);
      if(isKaon==true)vCutKRadius[0]->Fill(cutRadius[9]);
      if(isProton==true)vCutPrRadius[2]->Fill(cutRadius[9]);
      if(isProton==true)vCutPrRadius[0]->Fill(cutRadius[9]);
    }



    if(i<10)printf(" Plotta eve %3d %4d \n",i,evt);
    int nHit = 0;
    for(int j = 0; j < nedge; j++){
      if(i<10)printf(" %3d %3d (%4d, %4d) %7.2f %7.2f %2d ",j,pol[j],fiber[j],ch[j],nttw[j],dur[j],pol[j]);
      if(dur[j]>CUT_MIN_DUR && trigSig[j]==false){
        if(pol[j]==0)hHitStart[0]->Fill(nttw[j]);
        if(pol[j]==1)hHitEnd[0]->Fill(nttw[j]);
        if(innerPhoton[j]==true){
          if(pol[j]==0)hHitStart[1]->Fill(nttw[j]);
          if(pol[j]==1)hHitEnd[1]->Fill(nttw[j]);
          if(i<10)printf(" --> in\n");
        }
        if(outerPhoton[j]==true){
          if(pol[j]==0)hHitStart[2]->Fill(nttw[j]);
          if(pol[j]==1)hHitEnd[2]->Fill(nttw[j]);
          if(i<10)printf(" --> out \n");
        }
      }else{
        if(i<10)printf("bad \n");
      }
      if(goodHit[j]==true && trigSig[j]==false && pol[j]==0){
        hHitDur[0]->Fill(dur[j]);
        hHitCorrDur[0]->Fill(dur[j],nttw[j]);
        if(innerPhoton[j]==true){
          hHitDur[1]->Fill(dur[j]);
          hHitCorrDur[1]->Fill(dur[j],nttw[j]);
        }
        if(outerPhoton[j]==true){
          hHitDur[2]->Fill(dur[j]);
          hHitCorrDur[2]->Fill(dur[j],nttw[j]);
        }
        int ref=reference(pmt[j],ch[j]);
        if(goodPhoton[j]==true){
          hCalibVsCh[0]->Fill(ref,nttw[j]);
          hNotCalibVsCh[0]->Fill(ref,time[j]);
          if(innerPhoton[j]==true){
            hCalibVsCh[1]->Fill(ref,nttw[j]);
            hNotCalibVsCh[1]->Fill(ref,time[j]);
          }
          if(outerPhoton[j]==true){
            hCalibVsCh[2]->Fill(ref,nttw[j]);
            hNotCalibVsCh[2]->Fill(ref,time[j]);
          }
        }
      }
      if(goodPhoton[j]==true){
        hHitCoinc[0]->Fill(nttw[j]);
        if(innerPhoton[j]==true){
          hHitCoinc[1]->Fill(nttw[j]);
          vnRadius[1]->Fill(nr[j]);
          if(isPion==true){
            vPiRadius[1]->Fill(nr[j]);
            vPiMap[1]->Fill(x[j],y[j]);
          }
          if(isKaon==true){
            vKRadius[1]->Fill(nr[j]);
            vKMap[1]->Fill(x[j],y[j]);
          }
          if(isProton==true){
            vPrRadius[1]->Fill(nr[j]);
            vPrMap[1]->Fill(x[j],y[j]);
          }
        }
        if(outerPhoton[j]==true){
          hHitCoinc[2]->Fill(nttw[j]);
          vnRadius[2]->Fill(nr[j]);
          if(isPion==true){
            vPiRadius[2]->Fill(nr[j]);
            vPiMap[2]->Fill(x[j],y[j]);
          }
          if(isKaon==true){
            vKRadius[2]->Fill(nr[j]);
            vKMap[2]->Fill(x[j],y[j]);
          }
          if(isProton==true){
            vPrRadius[2]->Fill(nr[j]);
            vPrMap[2]->Fill(x[j],y[j]);
          }
        }
        hMap->Fill(x[j],y[j]);
        if(dur[j] > CUT_MIN_DUR){
          vMap[0]->Fill(nx[j],ny[j]);
          if(innerPhoton[j]==true)vMap[1]->Fill(nx[j],ny[j]);
          if(outerPhoton[j]==true)vMap[2]->Fill(nx[j],ny[j]);
          if(innerPhoton[j]==false && outerPhoton[j]==false){
            vMap[3]->Fill(nx[j],ny[j]);
            hbRadius->Fill(nr[j]);
            if(pol[j]==0){
              hbHitStart->Fill(nttw[j]);
              hbRadVsStart->Fill(nttw[j],nr[j]);
            }
            if(pol[j]==1){
              hbHitEnd->Fill(nttw[j]);
              //hbRadVsEnd->Fill(nttw[j],nr[j]);
            }
          }
        }
        hnMap->Fill(nx[j],ny[j]);
        hRadius->Fill(r[j]);
        vnRadius[0]->Fill(nr[j]);
        if(innerPhoton[j]==true)vGasRadPMT[pmt[j]]->Fill(nr[j]);
        if(outerPhoton[j]==true)vAeroRadPMT[pmt[j]]->Fill(nr[j]);
        if(isPion==true){
          vPiRadius[0]->Fill(nr[j]);
          vPiMap[0]->Fill(x[j],y[j]);
        }
        if(isKaon==true){
          vKRadius[0]->Fill(nr[j]);
          vKMap[0]->Fill(x[j],y[j]);
        }
        if(isProton==true){
          vPrRadius[0]->Fill(nr[j]);
          vPrMap[0]->Fill(x[j],y[j]);
        }
      }
    }
  }
  printEnd();
  fIn->Close();
}

//---------------------------------------------------
void displayMonitor(THeader *run){
  //---------------------------------------------------

  gErrorIgnoreLevel=kWarning;
  string out_pdf0 = Form("%s/output/plot/%s/displayMonitor.pdf[",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf  = Form("%s/output/plot/%s/displayMonitor.pdf",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf1 = Form("%s/output/plot/%s/displayMonitor.pdf]",run->suite.c_str(),run->outputDir.c_str());
  string out_root = Form("%s/output/plot/%s/displayMonitor.root",run->suite.c_str(),run->outputDir.c_str());

  TList *save = new TList();

  TCanvas *c0 = new TCanvas("c0","Inner and outer",1600,900);
  c0->Divide(2,2);
  c0->Draw();
  c0->cd(1);
  hbRadius->Draw();
  c0->cd(2);
  hbHitStart->GetXaxis()->SetRangeUser(run->timeOuMin,run->timeInMax);
  hbHitStart->Draw();
  c0->cd(3);
  gStyle->SetOptStat(0);
  hbRadVsStart->GetXaxis()->SetRangeUser(run->timeOuMin,run->timeInMax);
  TH2D *htmpbRadVsStart = (TH2D*) hbRadVsStart->Clone("htmpbRadVsStart");
  hbRadVsStart->GetYaxis()->SetRangeUser(20,60);
  hbRadVsStart->Draw("colz");
  c0->cd(4);
  htmpbRadVsStart->SetTitle("External radius vs start");
  htmpbRadVsStart->GetYaxis()->SetRangeUser(130,190);
  htmpbRadVsStart->Draw("colz");
  c0->Update();
  c0->Print(out_pdf0.c_str());
  c0->Print(out_pdf.c_str());


  TCanvas *c1 = new TCanvas("c1","Timing gas",1600,900);
  c1->Divide(4,2);
  c1->Draw();
  gStyle->SetOptStat(1);
  //  c1->Print(out_pdf0.c_str());
  for (int irad=1; irad<3; irad++){

    double tmin = run->timeInMin;
    if(run->timeOuMin<tmin)tmin=run->timeOuMin;
    double tmax = run->timeInMax;
    if(run->timeOuMax>tmax)tmax=run->timeOuMax;

    double tmin_rad = run->timeInMin; 
    double tmax_rad = run->timeInMax; 
    if(irad==2){ 
      tmin_rad = run->timeOuMin;
      tmax_rad = run->timeOuMax;
    }
    double dmin = CUT_MIN_DUR;
    //printf("Time ref (%7.2f:%7.2f) gas (%7.2f:%7.2f)   aerogel (%7.2f:%7.2f)",tmin,tmax,
    //		run->timeInMin,run->timeInMax,run->timeOuMin,run->timeOuMax);

    c1->cd(1);
    //Hit start and end, showing coincidence window extremes
    THStack *hsHitTime = new THStack("hsHitTime","Hit start and End - Duration > 35;[ns]");
    hHitEnd[irad]->SetLineColor(2);
    hsHitTime->Add(hHitStart[irad]);
    hsHitTime->Add(hHitEnd[irad]);
    hsHitTime->Draw("nostack");
    hsHitTime->GetXaxis()->SetRangeUser(run->timeOuMin-0.1*run->timeOuMin,run->timeInMax+110);
    //hsHitTime->GetXaxis()->SetRangeUser(300,600);
    hsHitTime->Draw("nostack");
    TLine *l1 = new TLine(tmin_rad,0,tmin_rad,hHitStart[0]->GetBinContent(hHitStart[0]->GetMaximumBin()));
    l1->SetLineColor(3);
    l1->Draw("same");
    TLine *l2 = new TLine(tmax_rad,0,tmax_rad,hHitStart[0]->GetBinContent(hHitStart[0]->GetMaximumBin()));
    l2->SetLineColor(3);
    l2->Draw("same");
    c1->Update();

    c1->cd(5);
    //Zoom on coincidence peak, fitted. Showing coincidence window extremes.
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    hHitCoinc[irad]->Draw();
    TF1 *fcoinc = new TF1("fcoinc","gaus(0)",200,400);
    fcoinc->SetParameters(20000,355,2);
    hHitCoinc[irad]->Fit(fcoinc,"Q","",run->timeOuMin,run->timeInMax);
    l1->SetY2(hHitCoinc[irad]->GetBinContent(hHitCoinc[irad]->GetMaximumBin())/2);
    l2->SetY2(hHitCoinc[irad]->GetBinContent(hHitCoinc[irad]->GetMaximumBin())/2);
    l1->Draw("same");
    l2->Draw("same");
    c1->Update();

    c1->cd(2);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);
    hx474->SetLineColor(2);
    hx519->SetLineColor(3);
    hx537->SetLineColor(4);
    THStack *hsBeamCherenkov= new THStack("hsBeamCherenkov","Beam cherenkov;time [ns];counts [#]");
    gPad->SetLogy();	
    hsBeamCherenkov->Add(hx474);
    hsBeamCherenkov->Add(hx519);
    hsBeamCherenkov->Add(hx537);
    hsBeamCherenkov->Draw("nostack");
    //hsBeamCherenkov->GetXaxis()->SetRangeUser(0,460);
    hsBeamCherenkov->Draw("nostack");
    gPad->BuildLegend(.3,.7,.7,.9);

    c1->cd(6);
    htrig->GetXaxis()->SetRangeUser(htrig->GetBinCenter(htrig->FindFirstBinAbove(1))-10,htrig->GetBinCenter(htrig->FindLastBinAbove(1))+10);
    htrig->Draw();

    c1->cd(3);
    //Trigger time plot. If I can, start and end.
    hHitDur[irad]->GetXaxis()->SetRangeUser(0,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(1))+5);
    hHitDur[irad]->Draw();
    TLine *l3 = new TLine(dmin,0,dmin,hHitDur[irad]->GetBinContent(hHitDur[0]->GetMaximumBin()));
    l3->SetLineColor(3);
    l3->Draw("same");

    c1->cd(4);
    gPad->SetLogz();
    hHitCorrDur[irad]->GetXaxis()->SetRangeUser(0,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(3))+5);
    hHitCorrDur[irad]->GetYaxis()->SetRangeUser(run->timeOuMin-20,run->timeInMax+20);
    hHitCorrDur[irad]->Draw("colz");
    TLine *l4 = new TLine(0,tmin_rad,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(3))+5,tmin_rad);
    l4->SetLineColor(3);
    l4->Draw("same");
    TLine *l5 = new TLine(0,tmax_rad,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(3))+5,tmax_rad);
    l5->SetLineColor(3);
    l5->Draw("same");
    TLine *l6 = new TLine(dmin,run->timeOuMin-20,dmin,run->timeInMax+20);
    l6->SetLineColor(3);
    l6->Draw("same");

    c1->cd(7);
    //Time not calibrated VS Channel
    hNotCalibVsCh[irad]->Print(Form("hNotCalib_%d.root",irad));
    TProfile *tp0 = hNotCalibVsCh[irad]->ProfileX();
    tp0->GetYaxis()->SetRangeUser(run->timeOuMin,run->timeInMax+5);
    tp0->Draw();

    c1->cd(8);
    hCalibVsCh[irad]->Print(Form("hCalib_%d.root",irad));
    TProfile *tp1 = hCalibVsCh[irad]->ProfileX();
    tp1->GetYaxis()->SetRangeUser(run->timeOuMin,run->timeInMax+5);
    tp1->Draw();

    c1->Update();
    c1->Print(out_pdf.c_str());

    TCanvas *cT  = new TCanvas(Form("cT%d",irad),"c0",1600,900);
    cT->Draw();
    hNotCalibVsCh[irad]->Draw("colz");
    cT->Update();
    cT->Print(Form("hNotCalivVsCh_%d.root",irad));
    hCalibVsCh[irad]->Draw("colz");
    cT->Update();
    cT->Print(Form("hCalivVsCh_%d.root",irad));
    cT->Close();
  }
  for (int j=1; j<=1024; j++){
    TProfile *tp_gas = hNotCalibVsCh[1]->ProfileX();
    TProfile *tp_aer = hNotCalibVsCh[2]->ProfileX();
    double off_gas = tp_gas->GetBinContent(j);
    double off_aer = tp_aer->GetBinContent(j);
    double err_gas = tp_gas->GetBinError(j);
    double err_aer = tp_aer->GetBinError(j);
    double zero = 0.0;
    /*printf(" jjj %3d %7.2f %7.2f %7.2f %7.2f ",j,off_gas,off_aer,err_gas,err_aer);
      if(off_gas!=0 && off_aer==0)printf(" prof %3d %7.2f \n",j,360-off_gas);
      if(off_gas==0 && off_aer!=0)printf(" prof %3d %7.2f \n",j,355-off_aer);
      if(off_gas!=0 && off_aer!=0){
      if(err_gas<=err_aer) printf(" prof %3d %7.2f \n",j,360-off_gas);
      if(err_gas>err_aer) printf(" prof %3d %7.2f \n",j,355-off_aer);
      }
      if(off_gas==0 && off_aer==0)printf(" prof %3d %7.2f \n",j,zero);*/
  }



  TCanvas *c2 = new TCanvas("c2","Radii canvas",1600,900);
  c2->Divide(3,3);
  c2->Draw();
  c2->cd(1);
  hBeamCh->SetMinimum(0);
  hBeamCh->Draw();
  //Cherenkov combination for pion 1
  c2->cd(2);
  //Cherenkov combination for pion 2
  c2->cd(3);
  //Cherenkov combination for pion+kaon

  c2->cd(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  //Fitted GAS radius for SP.
  vspRadius[4]->SetTitle("sp radius - gas");
  TH1D *cp00 = (TH1D*)vspRadius[4]->Clone("hspRadius_fitIn");
  TF1 *fspRadius_4=new TF1();
  applyFit(cp00,fspRadius_4,"fspRadius_4",false);
  cp00->Draw(); 
  c2->cd(5);
  //Fitted GAS radius for SPN.
  vspnRadius[4]->SetTitle("Corrected sp radius - gas");
  TH1D *cp01 = (TH1D*)vspnRadius[4]->Clone("hspnRadius_fitIn");
  TF1 *fspnRadius_4=new TF1();
  applyFit(cp01,fspnRadius_4,"fspnRadius_4",false);
  cp01->Draw(); 
  c2->cd(6);
  //Fitted GAS radius for CUT.
  vcutRadius[4]->SetTitle("RMS cut radius - gas");
  TH1D *cp02 = (TH1D*)vcutRadius[4]->Clone("hcutRadius_fitIn");
  TF1 *fcutRadius_4=new TF1();
  applyFit(cp02,fcutRadius_4,"fcutRadius_4",false);
  cp02->Draw(); 
  c2->cd(7);
  //Fitted AERO radius for SP.
  vspRadius[9]->SetTitle("Sp radius - aerogel");
  TH1D *cp10 = (TH1D*)vspRadius[9]->Clone("hspRadius_fitOut");
  TF1 *fspRadius_9 = new TF1();
  applyFit(cp10,fspRadius_9,"fspRadius_4",true);
  cp10->Draw();
  c2->cd(8);
  //Fitted AERO radius for SPN.
  vspnRadius[9]->SetTitle("Corrected sp radius - aerogel");
  TH1D *cp11 = (TH1D*)vspnRadius[9]->Clone("hspnRadius_fitOut");
  TF1 *fspnRadius_9 = new TF1();
  applyFit(cp11,fspnRadius_9,"fspnRadius_4",true);
  cp11->Draw();
  c2->cd(9);
  //Fitted AERO radius for CUT.
  vcutRadius[9]->SetTitle("RMS cut radius - aerogel");
  TH1D *cp12 = (TH1D*)vcutRadius[9]->Clone("hcutRadius_fitOut");
  TF1 *fcutRadius_9 = new TF1();
  applyFit(cp12,fcutRadius_9,"fcutRadius_4",true);
  cp12->Draw();
  c2->Update();
  c2->Print(out_pdf.c_str());


  TCanvas *c3 = new TCanvas("c3","Time and photon for radius",1600,900);
  c3->Divide(3,4);
  c3->cd(1);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);
  vspTime[4]->SetTitle("Mean time of the event - gas - Single Particle");
  vspTime[4]->Draw();
  c3->cd(4);
  vspPhoton[4]->SetTitle("# photons for particle - gas - Single Particle");
  vspPhoton[4]->Draw();
  c3->cd(7);
  vspTime[9]->SetTitle("Mean time of the event - aerogel - Single Particle");
  vspTime[9]->Draw();
  c3->cd(10);
  vspPhoton[9]->SetTitle("# photons for particle - aerogel - Single Particle");
  vspPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspPhoton[9]->Draw();
  c3->cd(2);
  vspnTime[4]->SetTitle("Mean time of the event - gas - Corrected Single Particle");
  vspnTime[4]->Draw();
  c3->cd(5);
  vspnPhoton[4]->SetTitle("# photons for particle - gas - Corrected Single Particle");
  vspnPhoton[4]->Draw();
  c3->cd(8);
  vspnTime[9]->SetTitle("Mean time of the event - aerogel - Corrected Single Particle");
  vspnTime[9]->Draw();
  c3->cd(11);
  vspnPhoton[9]->SetTitle("# photons for particle - aerogel - Corrected Single Particle");
  vspnPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspnPhoton[9]->Draw();
  c3->cd(3);
  vcutTime[4]->SetTitle("Mean time of the event - gas - Cut on RMS");
  vcutTime[4]->Draw();
  c3->cd(6);
  vcutPhoton[4]->SetTitle("# photons for particle - gas - Cut on RMS");
  vcutPhoton[4]->Draw();
  c3->cd(9);
  vcutTime[9]->SetTitle("Mean time of the event - aerogel - Cut on RMS");
  vcutTime[9]->Draw();
  c3->cd(12);
  vcutPhoton[9]->SetTitle("# photons for particle - aerogel - Cut on RMS");
  vcutPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vcutPhoton[9]->Draw();
  c3->Update();
  c3->Print(out_pdf.c_str());


  TCanvas *c4 = new TCanvas("c4","GEM canvas",1600,900);
  c4->Divide(2,2);
  c4->Draw();
  c4->cd(1);
  gPad->SetGrid();
  hUpGEM->GetXaxis()->SetRangeUser(-30,30);
  hUpGEM->GetYaxis()->SetRangeUser(-30,30);
  hUpGEM->Draw("colz");
  c4->cd(2);
  gPad->SetGrid();
  hDnGEM->GetXaxis()->SetRangeUser(-30,30);
  hDnGEM->GetYaxis()->SetRangeUser(-30,30);
  hDnGEM->Draw("colz");
  c4->cd(3);
  gPad->SetGrid();
  hBeam->Draw("colz");
  c4->cd(4);
  gPad->SetGrid();
  hBeamTheta->Draw("colz");
  c4->Update();
  c4->Print(out_pdf.c_str());

  //SIGMA vs Photon Number
  for(int i = 1; i < maxPhoGas; i++){
    if(vspSigPhoGas[i]->GetEntries()==0) continue;
    TF1 *fspGas = getFun(vspSigPhoGas[i-1],false);
    hspSigVsPhoGas->SetBinContent(i,fspGas->GetParameter(2));
    hspSigVsPhoGas->SetBinError(i,fspGas->GetParError(2));
  }
  for(int i = 1; i < maxPhoGas; i++){
    if(vspnSigPhoGas[i]->GetEntries()==0) continue;
    TF1 *fspnGas = getFun(vspnSigPhoGas[i-1],false);
    hspnSigVsPhoGas->SetBinContent(i,fspnGas->GetParameter(2));
    hspnSigVsPhoGas->SetBinError(i,fspnGas->GetParError(2));
  }
  for(int i = 1; i < maxPhoGas; i++){
    if(vcutSigPhoGas[i]->GetEntries()==0) continue;
    TF1 *fcutGas = getFun(vcutSigPhoGas[i-1],false);
    hcutSigVsPhoGas->SetBinContent(i,fcutGas->GetParameter(2));
    hcutSigVsPhoGas->SetBinError(i,fcutGas->GetParError(2));
  }
  for(int i = 1; i < maxPhoAero; i++){
    if(vspSigPhoAero[i]->GetEntries()==0) continue;
    TF1 *fspAero = getFun(vspSigPhoAero[i-1],true);
    hspSigVsPhoAero->SetBinContent(i,fspAero->GetParameter(2));
    hspSigVsPhoAero->SetBinError(i,fspAero->GetParError(2));
  }
  for(int i = 1; i < maxPhoAero; i++){
    if(vspnSigPhoAero[i]->GetEntries()==0) continue;
    TF1 *fspnAero = getFun(vspnSigPhoAero[i-1],true);
    hspnSigVsPhoAero->SetBinContent(i,fspnAero->GetParameter(2));
    hspnSigVsPhoAero->SetBinError(i,fspnAero->GetParError(2));
  }
  for(int i = 1; i < maxPhoAero; i++){
    if(vcutSigPhoAero[i]->GetEntries()==0) continue;
    TF1 *fcutAero = getFun(vcutSigPhoAero[i-1],true);
    hcutSigVsPhoAero->SetBinContent(i,fcutAero->GetParameter(2));
    hcutSigVsPhoAero->SetBinError(i,fcutAero->GetParError(2));
  }


  TCanvas *c5 = new TCanvas("c5","Single photon",1600,900);
  c5->Divide(4,2);
  c5->cd(1);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0); 
  vnRadius[1]->Fit("gaus","Q","",30,50);
  vnRadius[1]->Draw();
  c5->cd(5);
  vnRadius[2]->Fit("gaus","Q","",150,170);
  vnRadius[2]->Draw(); 
  c5->cd(2);
  fitSigma(hspSigVsPhoGas,false);
  hspSigVsPhoGas->GetXaxis()->SetRangeUser(0,25);
  hspSigVsPhoGas->GetYaxis()->SetRangeUser(0,1.5);
  hspSigVsPhoGas->Draw("E");
  c5->cd(6);
  fitSigma(hspSigVsPhoAero,true);
  hspSigVsPhoAero->GetXaxis()->SetRangeUser(0,25);
  hspSigVsPhoAero->GetYaxis()->SetRangeUser(0,4);
  hspSigVsPhoAero->Draw("E");
  c5->cd(3);
  fitSigma(hspnSigVsPhoGas,false);
  hspnSigVsPhoGas->GetXaxis()->SetRangeUser(0,25);
  hspnSigVsPhoGas->GetYaxis()->SetRangeUser(0,1.5);
  hspnSigVsPhoGas->Draw("E");
  c5->cd(7);
  fitSigma(hspnSigVsPhoAero,true);
  hspnSigVsPhoAero->GetXaxis()->SetRangeUser(0,25);
  hspnSigVsPhoAero->GetYaxis()->SetRangeUser(0,4);
  hspnSigVsPhoAero->Draw("E");
  c5->cd(4);
  fitSigma(hcutSigVsPhoGas,false);
  hcutSigVsPhoGas->GetXaxis()->SetRangeUser(0,25);
  hcutSigVsPhoGas->GetYaxis()->SetRangeUser(0,1.5);
  hcutSigVsPhoGas->Draw("E");
  c5->cd(8);
  fitSigma(hcutSigVsPhoAero,true);
  hcutSigVsPhoAero->GetXaxis()->SetRangeUser(0,25);
  hcutSigVsPhoAero->GetYaxis()->SetRangeUser(0,4);
  hcutSigVsPhoAero->Draw("E");
  c5->Update();
  c5->Print(out_pdf.c_str());


  TCanvas *c6 = new TCanvas("c7","Selected particle canvas",1600,900);
  c6->Divide(3,3);
  c6->cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vPiRadius[0]->Draw();
  c6->cd(2);
  vKRadius[0]->Draw();
  c6->cd(3);
  vPrRadius[0]->Draw();
  c6->cd(4);
  vPiRadius[1]->Fit("gaus","Q","",25,50);
  vPiRadius[1]->Draw();
  c6->cd(5);
  vKRadius[1]->Fit("gaus","Q","",25,50);
  vKRadius[1]->Draw();
  c6->cd(6);
  vPrRadius[1]->Fit("gaus","Q","",25,50);
  vPrRadius[1]->Draw();
  c6->cd(7);
  vPiRadius[2]->Fit("gaus","Q","",150,170);
  vPiRadius[2]->Draw();
  c6->cd(8);
  vKRadius[2]->Fit("gaus","Q","",150,170);
  vKRadius[2]->Draw();
  c6->cd(9);
  vPrRadius[2]->Fit("gaus","Q","",150,170);
  vPrRadius[2]->Draw();
  c6->Update();
  c6->Print(out_pdf.c_str());

  TCanvas *c8 = new TCanvas("c8","Selected particle THStack",1600,900);
  c8->Draw();
  c8->Divide(2);
  c8->cd(1);
  gPad->SetLogy();
  THStack *hsCut1 = new THStack("hsCut1","Single particle radius - Gas;Radius [mRad];Counts [#]");
  vCutKRadius[1]->SetLineColor(2);
  vCutPrRadius[1]->SetLineColor(3);
  hsCut1->Add(vCutPiRadius[1]);
  hsCut1->Add(vCutKRadius[1]);
  hsCut1->Add(vCutPrRadius[1]);
  hsCut1->Draw("nostack");
  TLegend *lhsCut1 = new TLegend(0.7,0.7,0.9,0.9);
  lhsCut1->AddEntry(vCutPiRadius[1],"Pion","lp");
  lhsCut1->AddEntry(vCutKRadius[1],"Kaon","lp");
  lhsCut1->AddEntry(vCutPrRadius[1],"Proton","lp");
  lhsCut1->Draw("same");
  c8->cd(2);
  gPad->SetLogy();
  THStack *hsCut2 = new THStack("hsCut2","Single particle radius - Aerogel;Radius [mRad];Counts [#]");
  vCutKRadius[2]->SetLineColor(2);
  vCutPrRadius[2]->SetLineColor(3);
  hsCut2->Add(vCutPiRadius[2]);
  hsCut2->Add(vCutKRadius[2]);
  hsCut2->Add(vCutPrRadius[2]);
  hsCut2->Draw("nostack");
  TLegend *lhsCut2 = new TLegend(0.1,0.7,0.1,0.9);
  lhsCut2->AddEntry(vCutPiRadius[2],"Pion","lp");
  lhsCut2->AddEntry(vCutKRadius[2],"Kaon","lp");
  lhsCut2->AddEntry(vCutPrRadius[2],"Proton","lp");
  lhsCut2->Draw("same");
  c8->Update();
  c8->Print(out_pdf.c_str());





  TCanvas *c7 = new TCanvas("c7","Selected particle single radius",1600,900);
  c7->Draw();
  c7->Divide(3,3);
  c7->cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vCutPiRadius[0]->Draw();
  c7->cd(2);
  vCutKRadius[0]->Draw();
  c7->cd(3);
  vCutPrRadius[0]->Draw();
  c7->cd(4);
  vCutPiRadius[1]->Fit("gaus","Q","",25,50);
  vCutPiRadius[1]->Draw();
  c7->cd(5);
  vCutKRadius[1]->SetLineColor(4);;
  vCutKRadius[1]->Fit("gaus","Q","",25,50);
  vCutKRadius[1]->Draw();
  c7->cd(6);
  vCutPrRadius[1]->SetLineColor(4);;
  vCutPrRadius[1]->Fit("gaus","Q","",25,50);
  vCutPrRadius[1]->Draw();
  c7->cd(7);
  vCutPiRadius[2]->Fit("gaus","Q","",150,170);
  vCutPiRadius[2]->Draw();
  c7->cd(8);
  vCutKRadius[2]->SetLineColor(4);;
  vCutKRadius[2]->Fit("gaus","Q","",150,170);
  vCutKRadius[2]->Draw();
  c7->cd(9);
  vCutPrRadius[2]->SetLineColor(4);;
  vCutPrRadius[2]->Fit("gaus","Q","",150,170);
  vCutPrRadius[2]->Draw();
  c7->Update();
  c7->Print(out_pdf.c_str());


  TCanvas *c9 = new TCanvas("c9","Selected particle maps",1600,900);
  c9->Divide(3,3);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1); 
  c9->cd(1);
  vPiMap[0]->SetTitle("Pions - Full map");
  vPiMap[0]->Draw("colz");
  c9->cd(2);
  vKMap[0]->SetTitle("Kaons - Full map");
  vKMap[0]->Draw("colz");
  c9->cd(3);
  vPrMap[0]->SetTitle("Protons - Full map");
  vPrMap[0]->Draw("colz");

  c9->cd(4);
  vPiMap[1]->SetTitle("Pions - Gas");
  vPiMap[1]->Draw("colz");
  c9->cd(5);
  vKMap[1]->SetTitle("Kaons - Gas");
  vKMap[1]->Draw("colz");
  c9->cd(6);
  vPrMap[1]->SetTitle("Protons - Gas");
  vPrMap[1]->Draw("colz");

  c9->cd(7);
  vPiMap[2]->SetTitle("Pions - Aerogel");
  vPiMap[2]->Draw("colz");
  c9->cd(8);
  vKMap[2]->SetTitle("Kaons - Aerogel");
  vKMap[2]->Draw("colz");
  c9->cd(9);
  vPrMap[2]->SetTitle("Protons - Aerogel");
  vPrMap[2]->Draw("colz");
  c9->Update();
  c9->Print(out_pdf.c_str());


  TCanvas *c10 = new TCanvas("c10","Map canvas",1600,900);
  c10->Divide(2);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1); 
  c10->Draw();
  c10->cd(1);
  gPad->SetLogz();
  hMap->Draw("colz");
  TEllipse *geoCutRing = new TEllipse(0,0,run->geoCut);
  geoCutRing->SetLineColor(2);
  geoCutRing->SetFillStyle(0);
  geoCutRing->Draw("same");
  c10->cd(2);
  hnMap->Draw("colz");
  c10->Update();
  c10->Print(out_pdf.c_str());

  TCanvas *c11 = new TCanvas("c11","Differential map canvas",1600,900);
  c11->Divide(2,2);
  c11->Draw();
  c11->cd(1);
  vMap[0]->SetTitle("All good photons map");
  vMap[0]->Draw("colz");
  geoCutRing->Draw("same");
  c11->cd(2);
  vMap[1]->SetTitle("Inner photons map");
  vMap[1]->Draw("colz");
  geoCutRing->Draw("same");
  c11->cd(3);
  vMap[2]->SetTitle("Outer photons map");
  vMap[2]->Draw("colz");
  geoCutRing->Draw("same");
  c11->cd(4);
  vMap[3]->SetTitle("Background photons map");
  vMap[3]->Draw("colz");
  geoCutRing->Draw("same");
  c11->Update();
  c11->Print(out_pdf.c_str());


  TCanvas *c12 = new TCanvas("c12","Differential map canvas",1600,900);
  c12->Divide(4,2);
  c12->Draw();
  c12->cd(1);
  vGasRadPMT[0]->SetTitle("PMT north - Gas");
  vGasRadPMT[0]->Draw();
  c12->cd(2);
  vGasRadPMT[2]->SetTitle("PMT est - Gas");
  vGasRadPMT[2]->Draw();
  c12->cd(3);
  vGasRadPMT[1]->SetTitle("PMT south - Gas");
  vGasRadPMT[1]->Draw();
  c12->cd(4);
  vGasRadPMT[3]->SetTitle("PMT west - Gas");
  vGasRadPMT[3]->Draw();
  c12->cd(5);
  vAeroRadPMT[0]->SetTitle("PMT north - Aero");
  vAeroRadPMT[0]->Draw();
  c12->cd(6);
  vAeroRadPMT[2]->SetTitle("PMT est - Aero");
  vAeroRadPMT[2]->Draw();
  c12->cd(7);
  vAeroRadPMT[1]->SetTitle("PMT south - Aero");
  vAeroRadPMT[1]->Draw();
  c12->cd(8);
  vAeroRadPMT[3]->SetTitle("PMT west - Aero");
  vAeroRadPMT[3]->Draw();
  c12->Update();
  c12->Print(out_pdf.c_str());

  c1->Print(out_pdf1.c_str());
  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();

  cout <<Form("You can open the monitoring plots typing: evince %s\n",out_pdf.c_str());
}


/*
   void displayMonitor2(THeader *run){
   gErrorIgnoreLevel=kWarning;
//Time distibution and coincidence peak zoom, rings before and after correction
string out_pdf0 = Form("%s/output/plot/%s/displayMonitor2.pdf[",run->suite.c_str(),run->outputDir.c_str());
string out_pdf = Form("%s/output/plot/%s/displayMonitor2.pdf",run->suite.c_str(),run->outputDir.c_str());
string out_pdf1 = Form("%s/output/plot/%s/displayMonitor2.pdf]",run->suite.c_str(),run->outputDir.c_str());
string out_root = Form("%s/output/plot/%s/displayMonitor2.root",run->suite.c_str(),run->outputDir.c_str());

cout <<out_root.c_str() <<endl;
TList *save = new TList();
save->Add(hCalibVsTime);
save->Add(hCalibVsCh);
save->Add(hNotCalibVsCh);
save->Add(hTimeVsCh);
save->Add(hTimeWalkVsTime);
save->Add(hDur);
save->Add(hStartVsDur);
save->Add(hSingPhotRad);
save->Add(hTimeRespectTrigger);
save->Add(hHitEnd);
save->Add(hHitStart[0]);


TCanvas *c4 = new TCanvas("c4","c4",1600,900);
c4->Divide(4,2);
c4->Draw();
c4->Print(out_pdf0.c_str());
c4->cd(1);
hCalibVsTime->Draw("colz");
c4->cd(5);
hTimeWalkVsTime->Draw("colz");
c4->cd(2);
THStack *hsHitTime = new THStack("hsHitTime","Hit start and End - Duration > 35;[ns]");
hHitEnd->SetLineColor(2);
hsHitTime->Add(hHitStart[0]);
hsHitTime->Add(hHitEnd);
hsHitTime->Draw("nostack");
hsHitTime->GetXaxis()->SetRangeUser(300,600);
hsHitTime->Draw("nostack");
TLine *l1 = new TLine(run->timeInMin,0,run->timeInMin,hTime->GetBinContent(hTime->GetMaximumBin()));
l1->SetLineColor(3);
l1->Draw("same");
TLine *l2 = new TLine(run->timeOuMax,0,run->timeOuMax,hTime->GetBinContent(hTime->GetMaximumBin()));
l2->SetLineColor(3);
l2->Draw("same");
c4->cd(6);
hDur->Draw();
c4->cd(3);
hStartVsDur->Draw("colz");
c4->cd(7);
gPad->SetLogz();
//hHitReco->Draw();
auto tp0 = hNotCalibVsCh->ProfileX();
tp0->GetYaxis()->SetRangeUser(350,380);
tp0->Draw();
//hNotCalibVsCh->Draw("colz");//hSingPhotRad->Draw();
c4->cd(4);
hMap->Draw("colz");
c4->cd(8);
auto tp1 = hCalibVsCh->ProfileX();
tp1->GetYaxis()->SetRangeUser(350,380);
tp1->Draw();
//  hCalibVsCh->Draw("colz");//hSingPhotRad->Draw();

c4->Update();
c4->Print(out_pdf.c_str());
*/

/*  TCanvas *c5 = new TCanvas("c5","c5",1600,900);
    for(int i = 0; i < 8; i++){
    for(int j = 0; j < 192; j++){
    if(vNotTime[i][j]->GetEntries()==0 && vCalTime[i][j]->GetEntries()==0)continue;
    THStack *hs = new THStack(Form("hs_%d_%d",i+4,j),Form("Calibrated and not - h_%d_%d; Time [ns]",i+4,j));
    vNotTime[i][j]->SetLineColor(2);
    hs->Add(vNotTime[i][j]);
    hs->Add(vCalTime[i][j]);
    hs->Draw("nostack");
    c5->Update();
    c5->Print(out_pdf.c_str());
    }
    }*/
/*c4->Print(out_pdf1.c_str());

  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();
////c1->Close();
}*/




//---------------------------------------------------
void fillHisto(THeader *run){
}
//---------------------------------------------------
/*
   gErrorIgnoreLevel=kWarning;
   TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
   TFile *fIn = new TFile (fName,"READ");
   TTree *t = (TTree*) fIn->Get("dRICH");

   int evt, nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10], cutPhoton[10], fiber[MAXDATA], ch[MAXDATA], slot[MAXDATA];
   double y[MAXDATA], x[MAXDATA], r[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nttw[MAXDATA], dur[MAXDATA], rsdRadius[MAXDATA]; 
   double rsdTime[MAXDATA], spRadius[10], spTime[10], spnRadius[10], spnTime[10], cutRadius[10], cutTime[10];
   float gx0, gy0, gx1, gy1, gxa, gya, gxtheta, gytheta;
   bool trigSig[MAXDATA], goodHit[MAXDATA], coincPhoton[MAXDATA],externalPhoton[MAXDATA], goodSP[10], goodSPN[10], goodCUT[10];
   t->SetBranchAddress("evt",&evt);
   t->SetBranchAddress("nedge",&nedge);
   t->SetBranchAddress("pol",&pol);
   t->SetBranchAddress("fiber",&fiber);
   t->SetBranchAddress("slot",&slot);
   t->SetBranchAddress("ch",&ch);
   t->SetBranchAddress("pmt",&pmt);
   t->SetBranchAddress("x",&x);
   t->SetBranchAddress("y",&y);
   t->SetBranchAddress("r",&r);
   t->SetBranchAddress("nx",&nx);
   t->SetBranchAddress("ny",&ny);
   t->SetBranchAddress("nr",&nr);
   t->SetBranchAddress("nttw",&nttw);
   t->SetBranchAddress("dur",&dur);
   t->SetBranchAddress("trigSig",&trigSig);
   t->SetBranchAddress("goodHit",&goodHit);
   t->SetBranchAddress("coincPhoton",&coincPhoton);
   t->SetBranchAddress("externalPhoton",&externalPhoton);
   t->SetBranchAddress("rsdRadius",&rsdRadius);
   t->SetBranchAddress("rsdTime",&rsdTime);

   t->SetBranchAddress("gx0",&gx0);
   t->SetBranchAddress("gy0",&gy0);
   t->SetBranchAddress("gx1",&gx1);
   t->SetBranchAddress("gy1",&gy1);
   t->SetBranchAddress("gxa",&gxa);
   t->SetBranchAddress("gya",&gya);
   t->SetBranchAddress("gxtheta",&gxtheta);
   t->SetBranchAddress("gytheta",&gytheta);

   t->SetBranchAddress("spRadius",&spRadius);
   t->SetBranchAddress("spTime",&spTime);
   t->SetBranchAddress("spPhoton",&spPhoton);
   t->SetBranchAddress("goodSP",&goodSP);

   t->SetBranchAddress("spnRadius",&spnRadius);
   t->SetBranchAddress("spnTime",&spnTime);
   t->SetBranchAddress("spnPhoton",&spnPhoton);
   t->SetBranchAddress("goodSPN",&goodSPN);

   t->SetBranchAddress("cutRadius",&cutRadius);
   t->SetBranchAddress("cutTime",&cutTime);
   t->SetBranchAddress("cutPhoton",&cutPhoton);
   t->SetBranchAddress("goodCUT",&goodCUT);


   cout <<Form("Filling histograms for run %d\n",run->runNum);;
   for(int i = 0; i < t->GetEntries(); i++){
   if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
   t->GetEntry(i);
   hUpGEM->Fill(gx0,gy0);
   hDnGEM->Fill(gx1,gy1);
   hBeam->Fill(gxa,gya);
   hBeamTheta->Fill(gxtheta,gytheta);
   int nHit = 0;
   hEdge->Fill(nedge);

   if(i<10)printf(" EDGE in mon event %3d %6d \n",i,evt);
   for(int j = 0; j < nedge; j++){
//if(goodSP[4]!=true && spPhoton[4] <= 2)continue;
//if(spRadius[9]>150)continue;
//cout <<evt <<" "<<spRadius[4] <<" " <<spRadius[9] <<" " <<j <<" " <<x[j] <<" " <<y[j] <<endl; 
if(pol[j]==0){
  hTime->Fill(nttw[j]);
  hMapNC->Fill(x[j],y[j]);
}
if(goodHit[j]==1){
  if(pol[j]==0){
    if(dur[j]>CUT_MIN_DUR && trigSig[j]==false)hHitStart[0]->Fill(nttw[j]);
    nHit++;
    //if(nttw[j] > run->timeInMin && nttw[j] < run->timeOuMax) hHitDur->Fill(dur[j]);
    //if(nttw[j] > 360 && nttw[j] < 370)hHitDur->Fill(dur[j]);
    if(trigSig[j]==false)hHitDur[0]->Fill(dur[j]);
    hHitCorrDur->Fill(dur[j],nttw[j]);
    if(i<10)printf(" h %3d (%4d %4d %4d) %7.2f %7.2f \n",j,slot[j],fiber[j],ch[j],nttw[j],dur[j]);
  }
  else{
    if(dur[j]>CUT_MIN_DUR && trigSig[j]==false)hHitEnd->Fill(nttw[j]);
  }
}

if(coincPhoton[j]==true){
  if(dur[j]>35)hCoinc->Fill(nttw[j]);
  hMap->Fill(x[j],y[j]);
  hRadius->Fill(r[j]);
  hnMap->Fill(nx[j],ny[j]);
  hnRadius->Fill(nr[j]);
  if(externalPhoton[j]==false){
    vrsdRadius[0]->Fill(rsdRadius[j]);
    vrsdTime[0]->Fill(rsdTime[j]);
  }else{
    vrsdRadius[1]->Fill(rsdRadius[j]);
    vrsdTime[1]->Fill(rsdTime[j]);
  }
}
}
hHitReco->Fill(nHit);
for(int j = 0; j < 10; j++){
  if(goodSP[j]==true && spPhoton[j] > 2){
    //if(spRadius[9]>150)continue;
    vspRadius[j]->Fill(spRadius[j]);
    vspTime[j]->Fill(spTime[j]);
    vspPhoton[j]->Fill(spPhoton[j]);
  }
  if(goodSPN[j]==true && spnPhoton[j] > 2){
    vspnRadius[j]->Fill(spnRadius[j]);
    vspnTime[j]->Fill(spnTime[j]);
    vspnPhoton[j]->Fill(spnPhoton[j]);
  }
  if(goodCUT[j]==true && cutPhoton[j] > 2){
    vcutRadius[j]->Fill(cutRadius[j]);
    vcutTime[j]->Fill(cutTime[j]);
    vcutPhoton[j]->Fill(cutPhoton[j]);
  }
}
for(int j = 0; j < nedge; j++){
  if(goodHit[j]==0 || pol[j]!=0)continue;
  if(trigSig[j]==true)continue;
  //for(int k = j; k < nedge; k++){
  //if(goodHit[k]==0 || pol[k]!=1)continue;
  //if(fiber[j]!=fiber[k] || ch[j]!=ch[k]) continue;
  //hHitCorr->Fill(nttw[j],nttw[k]);
  hHitCorr->Fill(nttw[j],nttw[j]+dur[j]);
  //break;
  //}
}
if(goodSP[4]==true && spPhoton[4] > minPhoGas && spPhoton[4]<maxPhoGas) vspSigPhoGas[spPhoton[4]]->Fill(spRadius[4]);
if(goodSP[9]==true && spPhoton[9] > minPhoAero && spPhoton[9]<maxPhoAero) vspSigPhoAero[spPhoton[9]]->Fill(spRadius[9]);
if(goodSPN[4]==true && spnPhoton[4] > minPhoGas && spnPhoton[4]<maxPhoGas) vspnSigPhoGas[spnPhoton[4]]->Fill(spnRadius[4]);
if(goodSPN[9]==true && spnPhoton[9] > minPhoAero && spnPhoton[9]<maxPhoAero) vspnSigPhoAero[spnPhoton[9]]->Fill(spnRadius[9]);
if(goodCUT[4]==true && cutPhoton[4] > minPhoGas && cutPhoton[4]<maxPhoGas) vcutSigPhoGas[cutPhoton[4]]->Fill(cutRadius[4]);
if(goodCUT[9]==true && cutPhoton[9] > minPhoAero && cutPhoton[9]<maxPhoAero) vcutSigPhoAero[cutPhoton[9]]->Fill(cutRadius[9]); 
}
printEnd();
fIn->Close();
}*/


void displayPhotonAnalysis(THeader *run){}
/*
   gErrorIgnoreLevel=kWarning;
//Time distibution and coincidence peak zoom, rings before and after correction
string out_pdf0 = Form("%s/output/plot/%s/displayPhotonAnalysis.pdf[",run->suite.c_str(),run->outputDir.c_str());
string out_pdf = Form("%s/output/plot/%s/displayPhotonAnalysis.pdf",run->suite.c_str(),run->outputDir.c_str());
string out_pdf1 = Form("%s/output/plot/%s/displayPhotonAnalysis.pdf]",run->suite.c_str(),run->outputDir.c_str());
string out_root = Form("%s/output/plot/%s/displayPhotonAnalysis.root",run->suite.c_str(),run->outputDir.c_str());


TList *save = new TList();
save->Add(hspSigVsPhoGas);
save->Add(hspnSigVsPhoGas);
save->Add(hcutSigVsPhoGas);
save->Add(hspSigVsPhoAero);
save->Add(hspnSigVsPhoAero);
save->Add(hcutSigVsPhoAero);
for(int i = 0; i < minPhoGas; i++){
save->Add(vspSigPhoGas[i]);
save->Add(vspnSigPhoGas[i]);
save->Add(vcutSigPhoGas[i]);
}
for(int i = 0; i < minPhoAero; i++){
save->Add(vspSigPhoAero[i]);
save->Add(vspnSigPhoAero[i]);
save->Add(vcutSigPhoAero[i]);
}

for(int i = 0; i < maxPhoGas; i++){
if(vspSigPhoGas[i]->GetEntries()==0) continue;
TF1 *fspGas = getFun(vspSigPhoGas[i],false);
hspSigVsPhoGas->SetBinContent(i,fspGas->GetParameter(2));
hspSigVsPhoGas->SetBinError(i,fspGas->GetParError(2));
}
for(int i = 0; i < maxPhoGas; i++){
if(vspnSigPhoGas[i]->GetEntries()==0) continue;
TF1 *fspnGas = getFun(vspnSigPhoGas[i],false);
hspnSigVsPhoGas->SetBinContent(i,fspnGas->GetParameter(2));
hspnSigVsPhoGas->SetBinError(i,fspnGas->GetParError(2));
}
for(int i = 0; i < maxPhoGas; i++){
if(vcutSigPhoGas[i]->GetEntries()==0) continue;
TF1 *fcutGas = getFun(vcutSigPhoGas[i],false);
hcutSigVsPhoGas->SetBinContent(i,fcutGas->GetParameter(2));
hcutSigVsPhoGas->SetBinError(i,fcutGas->GetParError(2));
}
for(int i = 0; i < maxPhoAero; i++){
if(vspSigPhoAero[i]->GetEntries()==0) continue;
TF1 *fspAero = getFun(vspSigPhoAero[i],true);
hspSigVsPhoAero->SetBinContent(i,fspAero->GetParameter(2));
hspSigVsPhoAero->SetBinError(i,fspAero->GetParError(2));
}
for(int i = 0; i < maxPhoAero; i++){
if(vspnSigPhoAero[i]->GetEntries()==0) continue;
TF1 *fspnAero = getFun(vspnSigPhoAero[i],true);
hspnSigVsPhoAero->SetBinContent(i,fspnAero->GetParameter(2));
hspnSigVsPhoAero->SetBinError(i,fspnAero->GetParError(2));
}
for(int i = 0; i < maxPhoAero; i++){
if(vcutSigPhoAero[i]->GetEntries()==0) continue;
TF1 *fcutAero = getFun(vcutSigPhoAero[i],true);
hcutSigVsPhoAero->SetBinContent(i,fcutAero->GetParameter(2));
hcutSigVsPhoAero->SetBinError(i,fcutAero->GetParError(2));
}

gStyle->SetOptStat(0);
gStyle->SetOptFit(1);

TCanvas *c1 = new TCanvas("c1","c1",1600,900);
c1->Draw();
c1->Print(out_pdf0.c_str());

fitSigma(hspSigVsPhoGas,false);
hspSigVsPhoGas->GetXaxis()->SetRangeUser(0,25);
hspSigVsPhoGas->GetYaxis()->SetRangeUser(0,2);
hspSigVsPhoGas->Draw("E");
c1->Update();
c1->Print(out_pdf.c_str());
fitSigma(hspnSigVsPhoGas,false);
hspnSigVsPhoGas->GetXaxis()->SetRangeUser(0,25);
hspnSigVsPhoGas->GetYaxis()->SetRangeUser(0,2);
hspnSigVsPhoGas->Draw("E");
c1->Update();
c1->Print(out_pdf.c_str());
fitSigma(hcutSigVsPhoGas,false);
hcutSigVsPhoGas->GetXaxis()->SetRangeUser(0,25);
hcutSigVsPhoGas->GetYaxis()->SetRangeUser(0,2);
hcutSigVsPhoGas->Draw("E");
c1->Update();
c1->Print(out_pdf.c_str());

fitSigma(hspSigVsPhoAero,true);
hspSigVsPhoAero->GetXaxis()->SetRangeUser(0,10);
hspSigVsPhoAero->GetYaxis()->SetRangeUser(0,5);
hspSigVsPhoAero->Draw("E");
c1->Update();
c1->Print(out_pdf.c_str());
fitSigma(hspnSigVsPhoAero,true);
hspnSigVsPhoAero->GetXaxis()->SetRangeUser(0,10);
hspnSigVsPhoAero->GetYaxis()->SetRangeUser(0,5);
hspnSigVsPhoAero->Draw("E");
c1->Update();
c1->Print(out_pdf.c_str());
fitSigma(hcutSigVsPhoAero,true);
hcutSigVsPhoAero->GetXaxis()->SetRangeUser(0,10);
hcutSigVsPhoAero->GetYaxis()->SetRangeUser(0,5);
hcutSigVsPhoAero->Draw("E");
c1->Update();
c1->Print(out_pdf.c_str());  
c1->Print(out_pdf1.c_str());

TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
save->Write();
fOut->Close();
c1->Close();
}*/


void displayBase(THeader *run){}
/*
   gErrorIgnoreLevel=kWarning;
//Time distibution and coincidence peak zoom, rings before and after correction
string out_pdf0 = Form("%s/output/plot/%s/displayBase.pdf[",run->suite.c_str(),run->outputDir.c_str());
string out_pdf = Form("%s/output/plot/%s/displayBase.pdf",run->suite.c_str(),run->outputDir.c_str());
string out_pdf1 = Form("%s/output/plot/%s/displayBase.pdf]",run->suite.c_str(),run->outputDir.c_str());
string out_root = Form("%s/output/plot/%s/displayBase.root",run->suite.c_str(),run->outputDir.c_str());

TList *save = new TList();
save->Add(hTime);
save->Add(hCoinc);
save->Add(hMap);
save->Add(hnMap);
save->Add(hRadius);
save->Add(hnRadius);

TCanvas *c1 = new TCanvas("c1","c1",900,900);
c1->Draw();
c1->Print(out_pdf0.c_str());

hTime->Draw();
TLine *l1 = new TLine(run->timeInMin,0,run->timeInMin,hTime->GetBinContent(hTime->GetMaximumBin()));
l1->SetLineColor(3);
l1->Draw("same");
TLine *l2 = new TLine(run->timeOuMax,0,run->timeOuMax,hTime->GetBinContent(hTime->GetMaximumBin()));
l2->SetLineColor(3);
l2->Draw("same");
c1->Update();
c1->Print(out_pdf.c_str());

gStyle->SetOptFit(1);
TF1 *f = new TF1("fCoinc","gaus(0)",0,1000);
f->SetParameters(25000,330,4);
hCoinc->Fit("fCoinc","Q","",run->timeInMin,run->timeOuMax);
hCoinc->Draw();
l1->SetY2(hCoinc->GetBinContent(hCoinc->GetMaximumBin()));
l2->SetY2(hCoinc->GetBinContent(hCoinc->GetMaximumBin()));
l1->Draw("same");
l2->Draw("same");
c1->Update();

//  gStyle->SetOptStat(1);
//  gStyle->SetOptFit(0);

c1->Print(out_pdf.c_str());
hMap->Draw("colz");
c1->Update();
c1->Print(out_pdf.c_str());
hnMap->Draw("colz");
c1->Print(out_pdf.c_str());
hRadius->Draw();
c1->Print(out_pdf.c_str());
hnRadius->Draw();
c1->Print(out_pdf.c_str());

gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
TH1D *cp0 = (TH1D*)hnRadius->Clone("hnRadius_fitIn");
TF1 *fInRadius = new TF1("fInRadius","gaus(0)",20,60);
save->Add(fInRadius);
fInRadius->SetParameters(5000,38,3);
cp0->Fit("fInRadius","Q","",30,45);
cp0->Draw();
c1->Print(out_pdf.c_str());


TH1D *cp1 = (TH1D*)hnRadius->Clone("hnRadius_fitOut");
TF1 *fOutRadius = new TF1("fOutRadius","gaus(0)",60,90);
save->Add(fOutRadius);
fOutRadius->SetParameters(5000,163,3);
cp1->Fit("fOutRadius","Q","",155,170);
cp1->Draw();
c1->Print(out_pdf.c_str());
c1->Print(out_pdf1.c_str());

TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
save->Write();
fOut->Close();
c1->Close();
}*/

void displayRSD(THeader *run){}
/*
   gErrorIgnoreLevel=kWarning;
   string out_pdf0 = Form("%s/output/plot/%s/displayRSD.pdf[",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf = Form("%s/output/plot/%s/displayRSD.pdf",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf1 = Form("%s/output/plot/%s/displayRSD.pdf]",run->suite.c_str(),run->outputDir.c_str());
   string out_root = Form("%s/output/plot/%s/displayRSD.root",run->suite.c_str(),run->outputDir.c_str());
   TList *save = new TList();
   save->Add(vrsdRadius[0]);
   save->Add(vrsdRadius[1]);
   save->Add(vrsdTime[0]);
   save->Add(vrsdTime[1]);
   TCanvas *c1 = new TCanvas("c1","c1",1600,900);
   c1->Divide(2);
   c1->Draw();
   c1->Print(out_pdf0.c_str());

   c1->cd(1);
   vrsdRadius[0]->SetTitle("Radius residui - inner ring"); 
   TLine *l1 = new TLine(run->cutRadiusInRMS,0,run->cutRadiusInRMS,vrsdRadius[0]->GetBinContent(1));
   l1->SetLineColor(3);
   vrsdRadius[0]->Draw();
   l1->Draw("same");
   c1->Draw();
   c1->Update();
   c1->cd(2);
   vrsdTime[0]->SetTitle("Time residui - inner ring"); 
   TLine *l2 = new TLine(run->cutTimeInRMS,0,run->cutTimeInRMS,vrsdTime[0]->GetBinContent(1));
   l2->SetLineColor(3);
   vrsdTime[0]->Draw();
   l2->Draw("same");
   c1->Update();
   c1->Print(out_pdf.c_str());

   c1->cd(1);
   vrsdRadius[1]->SetTitle("Radius residui - external ring"); 
   TLine *l3 = new TLine(run->cutRadiusOutRMS,0,run->cutRadiusOutRMS,vrsdRadius[1]->GetBinContent(1)); 
   l3->SetLineColor(3);
   vrsdRadius[1]->Draw();
   l3->Draw("same");
   c1->Update();
   c1->cd(2);
   vrsdTime[1]->SetTitle("Time residui - external ring"); 
   TLine *l4 = new TLine(run->cutTimeOutRMS,0,run->cutTimeOutRMS,vrsdTime[1]->GetBinContent(1)); 
   l4->SetLineColor(3);
   vrsdTime[1]->Draw();
   l4->Draw("same");
   c1->Update();
   c1->Print(out_pdf.c_str());

   c1->Print(out_pdf1.c_str());

   TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
   save->Write();
   fOut->Close();
   c1->Close();
   }*/

void displaySP(THeader *run){}
/*
   gErrorIgnoreLevel=kWarning;
   string out_pdf0 = Form("%s/output/plot/%s/displaySP.pdf[",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf = Form("%s/output/plot/%s/displaySP.pdf",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf1 = Form("%s/output/plot/%s/displaySP.pdf]",run->suite.c_str(),run->outputDir.c_str());
   string out_root = Form("%s/output/plot/%s/displaySP.root",run->suite.c_str(),run->outputDir.c_str());
   TList *save = new TList();
   save->Add(vspRadius[4]);
   save->Add(vspTime[4]);
   save->Add(vspPhoton[4]);
   save->Add(vspRadius[9]);
   save->Add(vspTime[9]);
   save->Add(vspPhoton[9]);

   TCanvas *c1 = new TCanvas("c1","c1",1600,900);
   c1->Draw();
   c1->Print(out_pdf0.c_str());
   c1->Divide(3);
   c1->cd(1);
   vspRadius[4]->SetTitle("Single particle radius - inner ring");
   vspRadius[4]->GetXaxis()->SetRangeUser(25,50);
   vspRadius[4]->Draw();
   c1->cd(2);
   vspTime[4]->SetTitle("Mean time of the event - inner ring");
   vspTime[4]->Draw();
   c1->cd(3);
   vspPhoton[4]->SetTitle("# photons for particle - inner ring");
   vspPhoton[4]->Draw();


   c1->Print(out_pdf.c_str());
   c1->cd(1); 
   vspRadius[9]->SetTitle("Single particle radius - external ring");
   vspRadius[9]->GetXaxis()->SetRangeUser(130,150);
   vspRadius[9]->Draw();
   c1->cd(2);
   vspTime[9]->SetTitle("Mean time of the event - external ring");
   vspTime[9]->Draw();
   c1->cd(3);
   vspPhoton[9]->SetTitle("# photons for particle - external ring");
   vspPhoton[9]->GetXaxis()->SetRangeUser(0,25);
   vspPhoton[9]->Draw();
   c1->Print(out_pdf.c_str());

   c1->Divide(1);
   c1->cd(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);

   TH1D *cp0 = (TH1D*)vspRadius[4]->Clone("hspRadius_fitIn");
   TF1 *fspRadius_4=new TF1();
   save->Add(fspRadius_4);
   applyFit(cp0,fspRadius_4,"fspRadius_4",false);
   cp0->Draw();
   c1->Print(out_pdf.c_str());

   TH1D *cp1 = (TH1D*)vspRadius[9]->Clone("hspRadius_fitOut");
   TF1 *fspRadius_9 = new TF1();
   save->Add(fspRadius_9);
   applyFit(cp1,fspRadius_9,"fspRadius_4",true);
   cp1->Draw();
   c1->Print(out_pdf.c_str());

   c1->Print(out_pdf1.c_str());

   TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
   save->Write();
   fOut->Close();
   c1->Close();
   }*/

void displaySPN(THeader *run){}
/*
   gErrorIgnoreLevel=kWarning;
   string out_pdf0 = Form("%s/output/plot/%s/displaySPN.pdf[",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf = Form("%s/output/plot/%s/displaySPN.pdf",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf1 = Form("%s/output/plot/%s/displaySPN.pdf]",run->suite.c_str(),run->outputDir.c_str());
   string out_root = Form("%s/output/plot/%s/displaySPN.root",run->suite.c_str(),run->outputDir.c_str());
   TList *save = new TList();
   save->Add(vspnRadius[4]);
   save->Add(vspnTime[4]);
   save->Add(vspnRadius[9]);
   save->Add(vspnTime[9]);

   TCanvas *c1 = new TCanvas("c1","c1",1600,900);
   c1->Draw();
   c1->Print(out_pdf0.c_str());

   c1->Divide(3);
   c1->cd(1);
   vspnRadius[4]->SetTitle("Single particle radius - inner ring - corrected");
   vspnRadius[4]->GetXaxis()->SetRangeUser(25,50);
   vspnRadius[4]->Draw();
   c1->cd(2);
   vspnTime[4]->SetTitle("Mean time of the event - inner ring - corrected");
   vspnTime[4]->Draw();
   c1->cd(3);
   vspnPhoton[4]->SetTitle("# photons for particle - inner ring - corrected");
   vspnPhoton[4]->Draw();
   c1->Print(out_pdf.c_str());


   c1->cd(1); 
   vspnRadius[9]->SetTitle("Single particle radius - external ring - corrected");
   vspnRadius[9]->GetXaxis()->SetRangeUser(130,180);
   vspnRadius[9]->Draw();
   c1->cd(2);
   vspnTime[9]->SetTitle("Mean time of the event - external ring - corrected");
   vspnTime[9]->Draw();
   c1->cd(3);
   vspnPhoton[9]->SetTitle("# photons for particle - external ring - corrected");
   vspnPhoton[9]->GetXaxis()->SetRangeUser(0,25);
   vspnPhoton[9]->Draw();
   c1->Print(out_pdf.c_str());

   c1->Divide(1);
   c1->cd(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);

   TH1D *cp0 = (TH1D*)vspnRadius[4]->Clone("hspnRadius_fitIn");
   TF1 *fspnRadius_4=new TF1();
   save->Add(fspnRadius_4);
   applyFit(cp0,fspnRadius_4,"fspnRadius_4",false);
   cp0->Draw();
   c1->Print(out_pdf.c_str());

   TH1D *cp1 = (TH1D*)vspnRadius[9]->Clone("hspnRadius_fitOut");
   TF1 *fspnRadius_9 = new TF1();
   save->Add(fspnRadius_9);
   applyFit(cp1,fspnRadius_9,"fspnRadius_4",true);
   cp1->Draw();
   c1->Print(out_pdf.c_str());



   c1->Print(out_pdf1.c_str());
   TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
   save->Write();
   fOut->Close();
   c1->Close();
   }*/

void displayCUT(THeader *run){}
/*
   gErrorIgnoreLevel=kWarning;
   string out_pdf0 = Form("%s/output/plot/%s/displayCUT.pdf[",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf = Form("%s/output/plot/%s/displayCUT.pdf",run->suite.c_str(),run->outputDir.c_str());
   string out_pdf1 = Form("%s/output/plot/%s/displayCUT.pdf]",run->suite.c_str(),run->outputDir.c_str());
   string out_root = Form("%s/output/plot/%s/displayCUT.root",run->suite.c_str(),run->outputDir.c_str());
   TList *save = new TList();
   save->Add(vcutRadius[4]);
   save->Add(vcutTime[4]);
   save->Add(vcutRadius[9]);
   save->Add(vcutTime[9]);

   TCanvas *c1 = new TCanvas("c1","c1",1600,900);
   c1->Draw();
   c1->Print(out_pdf0.c_str());


   c1->Divide(3);
   c1->cd(1);
   vcutRadius[4]->SetTitle("Single particle radius - inner ring - after rms cuts");
   vcutRadius[4]->GetXaxis()->SetRangeUser(25,50);
   vcutRadius[4]->Draw();
   c1->cd(2);
   vcutTime[4]->SetTitle("Mean time of the event - inner ring - after rms cuts");
   vcutTime[4]->Draw();
   c1->cd(3);
   vcutPhoton[4]->SetTitle("# photons for particle - inner ring - after rms cuts");
   vcutPhoton[4]->Draw();
   c1->Print(out_pdf.c_str());


   c1->cd(1); 
   vcutRadius[9]->SetTitle("Single particle radius - external ring - after rms cuts");
   vcutRadius[9]->GetXaxis()->SetRangeUser(130,180);
   vcutRadius[9]->Draw();
   c1->cd(2);
   vcutTime[9]->SetTitle("Mean time of the event - external ring - after rms cuts");
   vcutTime[9]->Draw();
   c1->cd(3);
   vcutPhoton[9]->SetTitle("# photons for particle - external ring - after rms cuts");
   vcutPhoton[9]->GetXaxis()->SetRangeUser(0,25);
   vcutPhoton[9]->Draw();
   c1->Print(out_pdf.c_str());


   c1->Divide(1);
   c1->cd(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);

   TH1D *cp0 = (TH1D*)vcutRadius[4]->Clone("hcutRadius_fitIn");
   TF1 *fcutRadius_4=new TF1();
   save->Add(fcutRadius_4);
   applyFit(cp0,fcutRadius_4,"fcutRadius_4",false);
   cp0->SetTitle("Single particle radius - Gas");
   cp0->Draw();
   c1->Print(out_pdf.c_str());

   TH1D *cp1 = (TH1D*)vcutRadius[9]->Clone("hcutRadius_fitOut");
   TF1 *fcutRadius_9 = new TF1();
   save->Add(fcutRadius_9);
   applyFit(cp1,fcutRadius_9,"fcutRadius_4",true);
   cp1->SetTitle("Single particle radius - Aerogel");
   cp1->Draw();
   c1->Print(out_pdf.c_str());

   c1->Print(out_pdf1.c_str());


   TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
   save->Write();
fOut->Close();
c1->Close();
}*/


void inizializePlot(THeader *run){
  //Create /output directory
  string out_dir = Form("%s/output/plot/%s",run->suite.c_str(),run->outputDir.c_str());
  cout <<out_dir.c_str() <<endl;
  if(gSystem->AccessPathName(out_dir.c_str())){
    if(std::system(Form("mkdir -p %s",out_dir.c_str()))){
      cout <<"[ERROR] Output directory not created\n";
      exit(EXIT_FAILURE);
    }
  }

  hTime = new TH1D("hTime all","Time distribution for all hits",1000,0,1000);

  for (int i=0; i<3; i++){
    TH1D *hStart = new TH1D(Form("hHitStart_%d",i),Form("Hit start - %d",i),2000,0,2000);
    hHitStart.push_back(hStart);
    TH1D *hEnd = new TH1D(Form("hHitEnd_%d",i),Form("Hit end - %d",i),2000,0,2000);
    hHitEnd.push_back(hEnd);
    TH1D *hDur = new TH1D(Form("hHitDur_%d",i),Form("Hit duration - %d",i),101,-.50,100.5);
    hHitDur.push_back(hDur);
    TH2D *hCorrDur = new TH2D(Form("hHitCorrDur_%d",i),Form("Hit correlation; Dur [ns]; Start [ns] - %d",i),71,-.5,70.5,2000,0,2000); 
    hHitCorrDur.push_back(hCorrDur);
    TH1D *hCoinc = new TH1D(Form("hHitCoinc_%d",i),Form("Coincidence peak - %d",i),5*(int)(run->timeInMax-run->timeOuMin)+10,run->timeOuMin-5,run->timeInMax+5);
    hHitCoinc.push_back(hCoinc);
    TH2D *hCal = new TH2D(Form("hCalibVsCh_%d",i),Form("#Delta time calibrated - trigger time Vs Channel - %d",i),1026,-0.5,1025.5,run->timeInMax-run->timeOuMin+60,run->timeOuMin-30,run->timeInMax+30);
    TH2D *hNoCal = new TH2D(Form("hNotCalibVsCh_%d",i),Form("#Delta time not calibrated - trigger time Vs Channel - %d",i),1026,-0.5,1025.5,run->timeInMax-run->timeOuMin+60,run->timeOuMin-30,run->timeInMax+30);
    hCalibVsCh.push_back(hCal);
    hNotCalibVsCh.push_back(hNoCal);
  }

  hHitCorr = new TH2D("hHitCorr","Hit correlation; Start [ns]; End[ns]",1000,0,1000,1000,0,1000); 
  hHitReco = new TH1D("hHitReco","Number of reconstructed hit for event",151,-0.5,150.5);
  hEdge = new TH1D("hEdge","Number of edge for event",151,-0.5,150.5);

  for(int i = 0; i < 4; i++){
    TH1D *hGasRadPMT = new TH1D(Form("hGasRadPMT%d",i),"hGasRadPMT;r [mRad]",160,20,60);
    vGasRadPMT.push_back(hGasRadPMT);
    TH1D *hAeroRadPMT = new TH1D(Form("hAeroRadPMT%d",i),"hAeroRadPMT,r [mRad]",320,120,200);
    vAeroRadPMT.push_back(hAeroRadPMT);
  }

  if(run->sensor=="MPPC") hMap = new TH2D("hMap","Hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMap = new TH2D("hMap",Form("Hit position MAPMT - run %d;x [mm];y [mm]",run->runNum),sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMap = new TH2D("hMap","Hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);

  for(int i = 0; i < 3; i++){
    TH2D *hPiMap = (TH2D*) hMap->Clone(Form("hPiMap%d",i));
    TH2D *hKMap = (TH2D*) hMap->Clone(Form("hKMap%d",i));
    TH2D *hPrMap = (TH2D*) hMap->Clone(Form("hPrMap%d",i));
    vPiMap.push_back(hPiMap);
    vKMap.push_back(hKMap);
    vPrMap.push_back(hPrMap);
  }
  for(int i = 0; i < 4; i++){
    TH2D *hMapD = (TH2D*) hMap->Clone(Form("hMap%d",i));
    vMap.push_back(hMapD);
  }
  hbRadius = new TH1D("hbRadius","Background radius - after corrections;r [mRad]",200,0,200);
  hbHitStart=new TH1D("hbHitStart","Background start time;Time [ns]",2000,0,2000);
  hbHitEnd=new TH1D("hbHitEnd","Background end time;Time [ns]",2000,0,2000);
  hbRadVsStart = new TH2D("hbRadVsStart","Internal radius Vs start; Start [ns]; Radius [mRad]",2000,0,2000,200,0,200); 
  //hbRadVsEnd;

  if(run->sensor=="MPPC") hMapNC = new TH2D("hMapNC","All hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMapNC = new TH2D("hMapNC","All hit position MAPMT;x [mm];y [mm]",sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMapNC = new TH2D("hMapNC","All hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);

  hnMap = new TH2D("hnMap","Corrected positions of hit;x [mm];y [mm]",180,-90,90,180,-90,90);

  //STD GEM HISTO
  /*hUpGEM = new TH2D("hUpGEM","Upstream GEM; x_0[mm];y_0[mm]",300,-60,60,300,-60,60);
    hDnGEM = new TH2D("hDnGEM","Dnstream GEM; x_0[mm];y_0[mm]",300,-60,60,300,-60,60);
    hBeam = new TH2D("hBeam","Beam profile at aerogel; x_0[mm];y_0[mm]",300,-60,60,300,-60,60);
    hBeamTheta = new TH2D("hBeamTheta","Beam divergence; x_0[mm];y_0[mm]",100,-.002,.002,100,-.002,.002);*/
  hUpGEM = new TH2D("hUpGEM","Upstream GEM; x_0[mm];y_0[mm]",3000,-600,600,3000,-600,600);
  hDnGEM = new TH2D("hDnGEM","Dnstream GEM; x_0[mm];y_0[mm]",3000,-600,600,3000,-600,600);
  hBeam = new TH2D("hBeam","Beam profile at aerogel; x_0[mm];y_0[mm]",150,-30,30,150,-30,30);
  hBeamTheta = new TH2D("hBeamTheta","Beam divergence; x_0[mm];y_0[mm]",400,-.002,.002,400,-.002,.002);


  hRadius = new TH1D("hRadius","Single photon radius - before corrections;r [mRad]",400,0,200);
  //	hnRadius = new TH1D("hnRadius","Single photon radius - after corrections;r [mRad]",400,0,200);
  TH1D *hnRadius0 = new TH1D("hnRadius0","Single photon radius - after corrections;r [mRad]",400,0,200);
  TH1D *hnRadius1 = new TH1D("hnRadius1","Single photon radius - Gas;r [mRad]",400,0,100);
  TH1D *hnRadius2 = new TH1D("hnRadius2","Single photon radius - Aerogel;r [mRad]",400,100,200);
  vnRadius.push_back(hnRadius0);
  vnRadius.push_back(hnRadius1);
  vnRadius.push_back(hnRadius2);
  TH1D *hPiRadius0 = new TH1D("hPiRadius0","Pi - Single photon radius - after corrections;r [mRad]",400,0,200);
  TH1D *hPiRadius1 = new TH1D("hPiRadius1","Pi - Single photon radius - Gas;r [mRad]",160,20,60);
  TH1D *hPiRadius2 = new TH1D("hPiRadius2","Pi - Single photon radius - Aerogel;r [mRad]",400,100,200);
  vPiRadius.push_back(hPiRadius0);
  vPiRadius.push_back(hPiRadius1);
  vPiRadius.push_back(hPiRadius2);
  TH1D *hKRadius0 = new TH1D("hKRadius0","K - Single photon radius - after corrections;r [mRad]",400,0,200);
  TH1D *hKRadius1 = new TH1D("hKRadius1","K - Single photon radius - Gas;r [mRad]",160,20,60);
  TH1D *hKRadius2 = new TH1D("hKRadius2","K - Single photon radius - Aerogel;r [mRad]",400,100,200);
  vKRadius.push_back(hKRadius0);
  vKRadius.push_back(hKRadius1);
  vKRadius.push_back(hKRadius2);
  TH1D *hPrRadius0 = new TH1D("hPrRadius0","Pr - Single photon radius - after corrections;r [mRad]",400,0,200);
  TH1D *hPrRadius1 = new TH1D("hPrRadius1","Pr - Single photon radius - Gas;r [mRad]",160,20,60);
  TH1D *hPrRadius2 = new TH1D("hPrRadius2","Pr - Single photon radius - Aerogel;r [mRad]",400,100,200);
  vPrRadius.push_back(hPrRadius0);
  vPrRadius.push_back(hPrRadius1);
  vPrRadius.push_back(hPrRadius2);


  TH1D *hCutPiRadius0 = new TH1D("hCutPiRadius0","Pi - Single particle radius - after corrections;r [mRad]",400,0,200);
  TH1D *hCutPiRadius1 = new TH1D("hCutPiRadius1","Pi - Single particle radius - Gas;r [mRad]",100,25,50);
  TH1D *hCutPiRadius2 = new TH1D("hCutPiRadius2","Pi - Single particle radius - Aerogel;r [mRad]",220,120,180);
  vCutPiRadius.push_back(hCutPiRadius0);
  vCutPiRadius.push_back(hCutPiRadius1);
  vCutPiRadius.push_back(hCutPiRadius2);
  TH1D *hCutKRadius0 = new TH1D("hCutKRadius0","K - Single particle radius - after corrections;r [mRad]",400,0,200);
  TH1D *hCutKRadius1 = new TH1D("hCutKRadius1","K - Single particle radius - Gas;r [mRad]",100,25,50);
  TH1D *hCutKRadius2 = new TH1D("hCutKRadius2","K - Single particle radius - Aerogel;r [mRad]",220,120,180);
  vCutKRadius.push_back(hCutKRadius0);
  vCutKRadius.push_back(hCutKRadius1);
  vCutKRadius.push_back(hCutKRadius2);
  TH1D *hCutPrRadius0 = new TH1D("hCutPrRadius0","Pr - Single particle radius - after corrections;r [mRad]",400,0,200);
  TH1D *hCutPrRadius1 = new TH1D("hCutPrRadius1","Pr - Single particle radius - Gas;r [mRad]",100,25,50);
  TH1D *hCutPrRadius2 = new TH1D("hCutPrRadius2","Pr - Single particle radius - Aerogel;r [mRad]",220,120,180);
  vCutPrRadius.push_back(hCutPrRadius0);
  vCutPrRadius.push_back(hCutPrRadius1);
  vCutPrRadius.push_back(hCutPrRadius2);

  for(int i = 0; i < 10; i++){
    TH1D *hspRadius = new TH1D(Form("hspRadius_%d",i),Form("Single particle radius - %d - before corrections;radius [mRad]",i),800,0,200);
    vspRadius.push_back(hspRadius);
    TH1D *hspTime = new TH1D(Form("hspTime_%d",i),Form("Single particle radius - %d - before corrections;time [ns]",i),5*(int)(run->timeInMax-run->timeOuMin),run->timeOuMin,run->timeInMax);
    vspTime.push_back(hspTime);
    TH1D *hspPhoton = new TH1D(Form("hspPhoton_%d",i),Form("Single particle radius - %d - before corrections;photon [#]",i),50,0,50);
    vspPhoton.push_back(hspPhoton);

    TH1D *hspnRadius = new TH1D(Form("hspnRadius_%d",i),Form("Single particle radius - %d - after corrections;radius [mRad]",i),800,0,200);
    vspnRadius.push_back(hspnRadius);
    TH1D *hspnTime = new TH1D(Form("hspnTime_%d",i),Form("Single particle radius - %d - after corrections;time [ns]",i),5*(int)(run->timeInMax-run->timeOuMin),run->timeOuMin,run->timeInMax);
    vspnTime.push_back(hspnTime);
    TH1D *hspnPhoton = new TH1D(Form("hspnPhoton_%d",i),Form("Single particle radius - %d - after corrections;photon [#]",i),50,0,50);
    vspnPhoton.push_back(hspnPhoton);

    TH1D *hcutRadius = new TH1D(Form("hcutRadius_%d",i),Form("Single particle radius - %d - after rms cuts;radius [mRad]",i),800,0,200);
    vcutRadius.push_back(hcutRadius);
    TH1D *hcutTime = new TH1D(Form("hcutTime_%d",i),Form("Single particle radius - %d - after rms cuts;time [ns]",i),5*(int)(run->timeInMax-run->timeOuMin),run->timeOuMin,run->timeInMax);
    vcutTime.push_back(hcutTime);
    TH1D *hcutPhoton = new TH1D(Form("hcutPhoton_%d",i),Form("Single particle radius - %d - after rms cuts;photon [#]",i),50,0,50);
    vcutPhoton.push_back(hcutPhoton);
  }

  for(int i = 0; i < 2; i++){
    TH1D *hrsdRadius = new TH1D(Form("hrsdRadius_%d",i),Form("Radius residui - %d;rsd_r [mRad];counts [#]",i),120,0,12);
    vrsdRadius.push_back(hrsdRadius);
    TH1D *hrsdTime = new TH1D(Form("hrsdTime_%d",i),Form("Time residui - %d;rsd_t [ns];counts [#]",i),120,0,12);
    vrsdTime.push_back(hrsdTime);
  }

  for(int i = 0; i < maxPhoGas; i++){
    TH1D *hspSigPhoGas= new TH1D(Form("hspSigPhoGas_%02d",i),Form("hspSigPhoGas_%02d",i),400,0,200);
    vspSigPhoGas.push_back(hspSigPhoGas);
    TH1D *hspnSigPhoGas= new TH1D(Form("hspnSigPhoGas_%02d",i),Form("hspnSigPhoGas_%02d",i),400,0,200);
    vspnSigPhoGas.push_back(hspnSigPhoGas);
    TH1D *hcutSigPhoGas= new TH1D(Form("hcutSigPhoGas_%02d",i),Form("hcutSigPhoGas_%02d",i),400,0,200);
    vcutSigPhoGas.push_back(hcutSigPhoGas);
  }
  for(int i = 0; i < maxPhoAero; i++){
    TH1D *hspSigPhoAero= new TH1D(Form("hspSigPhoAero_%02d",i),Form("hspSigPhoAero_%02d",i),400,0,200);
    vspSigPhoAero.push_back(hspSigPhoAero);
    TH1D *hspnSigPhoAero= new TH1D(Form("hspnSigPhoAero_%02d",i),Form("hspnSigPhoAero_%02d",i),400,0,200);
    vspnSigPhoAero.push_back(hspnSigPhoAero);
    TH1D *hcutSigPhoAero= new TH1D(Form("hcutSigPhoAero_%02d",i),Form("hcutSigPhoAero_%02d",i),400,0,200);
    vcutSigPhoAero.push_back(hcutSigPhoAero);
  }
  hspSigVsPhoGas = new TH1D("hspSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hspnSigVsPhoGas= new TH1D("hsnpSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hcutSigVsPhoGas = new TH1D("hcutSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hspSigVsPhoAero = new TH1D("hspSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
  hspnSigVsPhoAero= new TH1D("hsnpSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
  hcutSigVsPhoAero = new TH1D("hcutSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);



  //hCoinc = new TH1D("hCoinc","Coincidence peak",5*(int)(run->timeOuMax-run->timeInMin)+10,run->timeInMin-5,run->timeOuMax+5);

  ////// MONITORING PLOT DEFINITION
  //hCalibVsTime = new TH2D("hCalibVsTime","Calibrated time vs acquired time",(int)(run->timeOuMax-run->timeInMin)+11,run->timeInMin-5.5,run->timeOuMax+5.5,(int)(run->timeOuMax-run->timeInMin)+11,run->timeInMin-5.5,run->timeOuMax+5.5);
  //hTimeWalkVsTime = new TH2D("hTimeWalkVsTime","TimeWalk corrected time vs acquired time",(int)(run->timeOuMax-run->timeInMin)+11,run->timeInMin-5.5,run->timeOuMax+5.5,5*(int)(run->timeOuMax-run->timeInMin)+50,run->timeInMin-5,run->timeOuMax+5);
  hCalibVsTime = new TH2D("hCalibVsTime","Calibrated time vs acquired time",51,349.5,400.1,(int)(run->timeInMax-run->timeOuMin)+11,run->timeOuMin-5.5,run->timeInMax+5.5);
  hTimeVsCh = new TH2D("hTimeVsCh","#Delta time not calibrated - trigger time Vs Channel",1026,-0.5,1025.5,70,330,400);
  hTimeWalkVsTime = new TH2D("hTimeWalkVsTime","TimeWalk corrected time vs acquired time",51,349.5,401.5,5*(int)(run->timeInMax-run->timeOuMin)+50,run->timeOuMin-5,run->timeInMax+5);
  //hDur = new TH1D("hDur","Hit duration",101,-.50,100.5);
  hStartVsDur = new TH2D("hStartVsDur","DuratioVsStart; Start [ns]; Duration[ns]",(int)(run->timeInMax-run->timeOuMin)+11,run->timeOuMin-5.5,run->timeInMax+5.5,101,-.5,100.5); 
  hSingPhotRad = new TH1D("hSingPhotRad","Single photon radius; Radius [mRad]",400,0,200);

  for(int i = 0; i < 32; i++){
    vector<TH1D*> tmpNot;
    vector<TH1D*> tmpCal;
    for(int j = 0; j < 1025; j++){
      TH1D* hNot = new TH1D(Form("hNot_%d_%d",i+4,j),"Time not calibrated",1001,-.5,1000.5);
      tmpNot.push_back(hNot);
      TH1D* hCal = new TH1D(Form("hCal_%d_%d",i+4,j),"Time calibrated",1001,-.5,1000.5);
      tmpCal.push_back(hCal);
    }
    vNotTime.push_back(tmpNot);
    vCalTime.push_back(tmpCal);
  } 
  hx474=new TH1D("hx474","Cherenkov x474",2000,0,2000);
  hx519=new TH1D("hx519","Cherenkov x519",2000,0,2000);
  hx537=new TH1D("hx537","Cherenkov x537",2000,0,2000);
  htrig=new TH1D("htrig","Trigger;time [ns]",2000,0,2000);


  hPiRad = new TH1D("hPiRad","#pi radius - gas; radius [mRad]",200,0,100);
  hKRad = new TH1D("hKRad","K radius - gas; radius [mRad]",200,0,100);
  hPrRad = new TH1D("hPrRad","p radius - gas; radius [mRad]",200,0,100);

  hBeamCh = new TH1D("hBeamCh","Beam Cherenkov",3,-.5,2.5);
}



