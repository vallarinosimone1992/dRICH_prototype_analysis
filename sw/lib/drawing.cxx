#define MAXDATA 10000
#include <iostream>
#include <vector>
//#include <filesystem>
//namespace fs = std::filesystem;

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
#include <TSystem.h>

#include "definition.h"
#include "utility.h"
#include "computing.h"

static TH1D *hTime;
static TH1D *hCoinc;
static TH1D *hHitStart;
static TH1D *hHitEnd;
static TH1D *hHitDur;
static TH2D *hHitCorr;
static TH2D *hHitCorrDur;
static TH2D *hMap;
static TH2D *hMapNC;
static TH2D *hnMap;
static TH1D *hRadius;
static TH1D *hnRadius;

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

static int minPhoGas=2;
static int maxPhoGas=50;
static int minPhoAero=2;
static int maxPhoAero=15;



void fillHisto(THeader *run){
  gErrorIgnoreLevel=kWarning;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int evt, nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10], cutPhoton[10], fiber[MAXDATA], ch[MAXDATA];
  double y[MAXDATA], x[MAXDATA], r[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nttw[MAXDATA], dur[MAXDATA], rsdRadius[MAXDATA], rsdTime[MAXDATA], spRadius[10], spTime[10], spnRadius[10], spnTime[10], cutRadius[10], cutTime[10];
  float gx0, gy0, gx1, gy1, gxa, gya, gxtheta, gytheta;
  bool trigSig[MAXDATA], goodHit[MAXDATA], coincPhoton[MAXDATA],outerPhoton[MAXDATA], goodSP[10], goodSPN[10], goodCUT[10];
  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("fiber",&fiber);
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


  cout <<Form("Filling histograms for run %d\n",run->runNum);;
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    hUpGEM->Fill(gx0,gy0);
    hDnGEM->Fill(gx1,gy1);
    hBeam->Fill(gxa,gya);
    hBeamTheta->Fill(gxtheta,gytheta);
    int nHit = 0;
    hEdge->Fill(nedge);
    for(int j = 0; j < nedge; j++){
      //if(goodSP[4]!=true && spPhoton[4] <= 2)continue;
      if(spRadius[9]>150)continue;
      cout <<evt <<" "<<spRadius[4] <<" " <<spRadius[9] <<" " <<j <<" " <<x[j] <<" " <<y[j] <<endl; 
      if(pol[j]==0){
        hTime->Fill(nttw[j]);
        hMapNC->Fill(x[j],y[j]);
      }
      if(goodHit[j]==1){
        if(pol[j]==0){
          if(dur[j]>35 && trigSig[j]==false)hHitStart->Fill(nttw[j]);
          nHit++;
	        //if(nttw[j] > run->timeMin && nttw[j] < run->timeMax) hHitDur->Fill(dur[j]);
	        //if(nttw[j] > 360 && nttw[j] < 370)hHitDur->Fill(dur[j]);
	        if(trigSig[j]==false)hHitDur->Fill(dur[j]);
	        hHitCorrDur->Fill(dur[j],nttw[j]);
          }
        else{
          if(dur[j]>35 && trigSig[j]==false)hHitEnd->Fill(nttw[j]);
          }
      }

      if(coincPhoton[j]==true){
        if(dur[j]>35)hCoinc->Fill(nttw[j]);
        hMap->Fill(x[j],y[j]);
        hRadius->Fill(r[j]);
        hnMap->Fill(nx[j],ny[j]);
        hnRadius->Fill(nr[j]);
        if(outerPhoton[j]==false){
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
      for(int k = j; k < nedge; k++){
        if(goodHit[k]==0 || pol[k]!=1)continue;
        if(fiber[j]!=fiber[k] || ch[j]!=ch[k]) continue;
        hHitCorr->Fill(nttw[j],nttw[k]);
        break;
      }
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
}

void displayMonitor(THeader *run){
  gErrorIgnoreLevel=kWarning;
  //Time distibution and coincidence peak zoom, rings before and after correction
  string out_pdf0 = Form("%s/output/plot/%s/displayMonitor.pdf[",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf = Form("%s/output/plot/%s/displayMonitor.pdf",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf1 = Form("%s/output/plot/%s/displayMonitor.pdf]",run->suite.c_str(),run->outputDir.c_str());
  string out_root = Form("%s/output/plot/%s/displayMonitor.root",run->suite.c_str(),run->outputDir.c_str());
  cout <<"Check: " <<out_pdf.c_str() << endl;

  TList *save = new TList();
  save->Add(hTime);
  save->Add(hCoinc);
  save->Add(hHitStart);
  save->Add(hHitEnd);
  save->Add(hHitDur);
  save->Add(hHitCorr);
  save->Add(hHitCorrDur);
  save->Add(hHitReco);
  save->Add(hEdge);
  save->Add(hMap);
  save->Add(hnMap);
  save->Add(hRadius);
  save->Add(hnRadius);

  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  c1->Divide(4,2);
  c1->Draw();
  c1->Print(out_pdf0.c_str());
  c1->cd(1);
  THStack *hsHitTime = new THStack("hsHitTime","Hit start and End - Duration > 35;[ns]");
  hHitEnd->SetLineColor(2);
  hsHitTime->Add(hHitStart);
  hsHitTime->Add(hHitEnd);
  hsHitTime->Draw("nostack");
  hsHitTime->GetXaxis()->SetRangeUser(300,600);
  hsHitTime->Draw("nostack");
  //hTime->Draw();
  TLine *l1 = new TLine(run->timeMin,0,run->timeMin,hTime->GetBinContent(hTime->GetMaximumBin()));
  l1->SetLineColor(3);
  l1->Draw("same");
  TLine *l2 = new TLine(run->timeMax,0,run->timeMax,hTime->GetBinContent(hTime->GetMaximumBin()));
  l2->SetLineColor(3);
  l2->Draw("same");
  c1->Update();
  c1->cd(5);
  hCoinc->Draw();
  TF1 *fcoinc = new TF1("fcoinc","gaus(0)",200,400);
  fcoinc->SetParameters(20000,355,2);
  hCoinc->Fit(fcoinc,"Q","",run->timeMin,run->timeMax);
  l1->SetY2(hCoinc->GetBinContent(hCoinc->GetMaximumBin()));
  l2->SetY2(hCoinc->GetBinContent(hCoinc->GetMaximumBin()));
  l1->Draw("same");
  l2->Draw("same");
  c1->cd(2);
  hHitDur->Draw();
  c1->cd(3);
  hHitCorr->Draw("colz");
  c1->cd(6);
  hHitReco->Draw();
  c1->cd(7);
  hHitCorrDur->Draw("colz");
  c1->cd(4);
  hMap->Draw("colz");
  c1->cd(8);
  hMapNC->Draw("colz");

  c1->Update();
  c1->Print(out_pdf.c_str());

  TCanvas *c2 = new TCanvas("c2","c2",1600,900);
  c2->Divide(3,2);
  c2->Draw();
  c2->cd(2);
  vspTime[4]->SetTitle("Mean time of the event - gas");
  vspTime[4]->Draw();
  c2->cd(3);
  vspPhoton[4]->SetTitle("# photons for particle - gas");
  vspPhoton[4]->Draw();
  c2->cd(5);
  vspTime[9]->SetTitle("Mean time of the event - aerogel");
  vspTime[9]->Draw();
  c2->cd(6);
  vspPhoton[9]->SetTitle("# photons for particle - aerogel");
  vspPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspPhoton[9]->Draw();
  c2->Update();


  c2->cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vspRadius[4]->SetTitle("Single particle radius - gas");
  TH1D *cp0 = (TH1D*)vspRadius[4]->Clone("hspRadius_fitIn");
  TF1 *fspRadius_4=new TF1();
  save->Add(fspRadius_4);
  applyFit(cp0,fspRadius_4,"fspRadius_4",false);
  cp0->Draw(); 
  c2->cd(4); 
  vspRadius[9]->SetTitle("Single particle radius - aerogel");
  TH1D *cp1 = (TH1D*)vspRadius[9]->Clone("hspRadius_fitOut");
  TF1 *fspRadius_9 = new TF1();
  save->Add(fspRadius_9);
  applyFit(cp1,fspRadius_9,"fspRadius_4",true);
  cp1->Draw();

  c2->Update();
  c2->Print(out_pdf.c_str());

  TCanvas *c3 = new TCanvas("c3","GEM canvas",1600,900);
  c3->Draw();
  c3->Divide(2,2);
  c3->cd(1);
  hUpGEM->Draw("colz");
  c3->cd(2);
  hDnGEM->Draw("colz");
  c3->cd(3);
  hBeam->Draw("colz");
  c3->cd(4);
  hBeamTheta->Draw("colz");
  c3->Update();
  c3->Print(out_pdf.c_str());

  c2->Print(out_pdf1.c_str());

  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();
  ////c1->Close();
}

void displayPhotonAnalysis(THeader *run){
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
}


void displayBase(THeader *run){
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
  TLine *l1 = new TLine(run->timeMin,0,run->timeMin,hTime->GetBinContent(hTime->GetMaximumBin()));
  l1->SetLineColor(3);
  l1->Draw("same");
  TLine *l2 = new TLine(run->timeMax,0,run->timeMax,hTime->GetBinContent(hTime->GetMaximumBin()));
  l2->SetLineColor(3);
  l2->Draw("same");
  c1->Update();
  c1->Print(out_pdf.c_str());

  gStyle->SetOptFit(1);
  TF1 *f = new TF1("fCoinc","gaus(0)",0,1000);
  f->SetParameters(25000,330,4);
  hCoinc->Fit("fCoinc","Q","",run->timeMin,run->timeMax);
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
}

void displayRSD(THeader *run){
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
  vrsdRadius[1]->SetTitle("Radius residui - outer ring"); 
  TLine *l3 = new TLine(run->cutRadiusOutRMS,0,run->cutRadiusOutRMS,vrsdRadius[1]->GetBinContent(1)); 
  l3->SetLineColor(3);
  vrsdRadius[1]->Draw();
  l3->Draw("same");
  c1->Update();
  c1->cd(2);
  vrsdTime[1]->SetTitle("Time residui - outer ring"); 
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
}

void displaySP(THeader *run){
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
  vspRadius[9]->SetTitle("Single particle radius - outer ring");
  vspRadius[9]->GetXaxis()->SetRangeUser(130,180);
  vspRadius[9]->Draw();
  c1->cd(2);
  vspTime[9]->SetTitle("Mean time of the event - outer ring");
  vspTime[9]->Draw();
  c1->cd(3);
  vspPhoton[9]->SetTitle("# photons for particle - outer ring");
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
}

void displaySPN(THeader *run){
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
  vspnRadius[9]->SetTitle("Single particle radius - outer ring - corrected");
  vspnRadius[9]->GetXaxis()->SetRangeUser(130,180);
  vspnRadius[9]->Draw();
  c1->cd(2);
  vspnTime[9]->SetTitle("Mean time of the event - outer ring - corrected");
  vspnTime[9]->Draw();
  c1->cd(3);
  vspnPhoton[9]->SetTitle("# photons for particle - outer ring - corrected");
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
}

void displayCUT(THeader *run){
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
  vcutRadius[9]->SetTitle("Single particle radius - outer ring - after rms cuts");
  vcutRadius[9]->GetXaxis()->SetRangeUser(130,180);
  vcutRadius[9]->Draw();
  c1->cd(2);
  vcutTime[9]->SetTitle("Mean time of the event - outer ring - after rms cuts");
  vcutTime[9]->Draw();
  c1->cd(3);
  vcutPhoton[9]->SetTitle("# photons for particle - outer ring - after rms cuts");
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
}


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

  hTime = new TH1D("hTime","Time distribution for all hits",1000,0,1000);
  hCoinc = new TH1D("hCoinc","Coincidence peak",5*(int)(run->timeMax-run->timeMin)+10,run->timeMin-5,run->timeMax+5);
  hHitStart = new TH1D("hHitStart","Hit start",1000,0,1000);
  hHitEnd = new TH1D("hHitEnd","Hit end",1000,0,1000);
  hHitDur = new TH1D("hHitDur","Hit duration",101,-.50,100.5);
  hHitCorr = new TH2D("hHitCorr","Hit correlation; Start [ns]; End[ns]",1000,0,1000,1000,0,1000); 
  hHitCorrDur = new TH2D("hHitCorrDur","Hit correlation; Dur [ns]; Start [ns]",71,-.5,70.5,100,330,430); 
  hHitReco = new TH1D("hHitReco","Number of reconstructed hit for event",151,-0.5,150.5);
  hEdge = new TH1D("hEdge","Number of edge for event",151,-0.5,150.5);


  if(run->sensor=="MPPC") hMap = new TH2D("hMap","Hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMap = new TH2D("hMap","Hit position MAPMT;x [mm];y [mm]",sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMap = new TH2D("hMap","Hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);

  if(run->sensor=="MPPC") hMapNC = new TH2D("hMapNC","All hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMapNC = new TH2D("hMapNC","All hit position MAPMT;x [mm];y [mm]",sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMapNC = new TH2D("hMapNC","All hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);

  hnMap = new TH2D("hnMap","Corrected positions of hit;x [mm];y [mm]",180,-90,90,180,-90,90);

  hUpGEM = new TH2D("hUpGEM","Upstream GEM; x_0[mm];y_0[mm]",300,-60,60,300,-60,60);
  hDnGEM = new TH2D("hDnGEM","Dnstream GEM; x_0[mm];y_0[mm]",300,-60,60,300,-60,60);
  hBeam = new TH2D("hBeam","Beam profile at aerogel; x_0[mm];y_0[mm]",300,-60,60,300,-60,60);
  hBeamTheta = new TH2D("hBeamTheta","Beam divergence; x_0[mm];y_0[mm]",100,-.002,.002,100,-.002,.002);


  hRadius = new TH1D("hRadius","Single photon radius - before corrections;r [mRad]",400,0,200);
  hnRadius = new TH1D("hnRadius","Single photon radius - after corrections;r [mRad]",400,0,200);

  for(int i = 0; i < 10; i++){
    TH1D *hspRadius = new TH1D(Form("hspRadius_%d",i),Form("Single particle radius - %d - before corrections;radius [mRad]",i),400,0,200);
    vspRadius.push_back(hspRadius);
    TH1D *hspTime = new TH1D(Form("hspTime_%d",i),Form("Single particle radius - %d - before corrections;time [ns]",i),5*(int)(run->timeMax-run->timeMin),run->timeMin,run->timeMax);
    vspTime.push_back(hspTime);
    TH1D *hspPhoton = new TH1D(Form("hspPhoton_%d",i),Form("Single particle radius - %d - before corrections;photon [#]",i),50,0,50);
    vspPhoton.push_back(hspPhoton);

    TH1D *hspnRadius = new TH1D(Form("hspnRadius_%d",i),Form("Single particle radius - %d - after corrections;radius [mRad]",i),400,0,200);
    vspnRadius.push_back(hspnRadius);
    TH1D *hspnTime = new TH1D(Form("hspnTime_%d",i),Form("Single particle radius - %d - after corrections;time [ns]",i),5*(int)(run->timeMax-run->timeMin),run->timeMin,run->timeMax);
    vspnTime.push_back(hspnTime);
    TH1D *hspnPhoton = new TH1D(Form("hspnPhoton_%d",i),Form("Single particle radius - %d - after corrections;photon [#]",i),50,0,50);
    vspnPhoton.push_back(hspnPhoton);

    TH1D *hcutRadius = new TH1D(Form("hcutRadius_%d",i),Form("Single particle radius - %d - after rms cuts;radius [mRad]",i),400,0,200);
    vcutRadius.push_back(hcutRadius);
    TH1D *hcutTime = new TH1D(Form("hcutTime_%d",i),Form("Single particle radius - %d - after rms cuts;time [ns]",i),5*(int)(run->timeMax-run->timeMin),run->timeMin,run->timeMax);
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
}
