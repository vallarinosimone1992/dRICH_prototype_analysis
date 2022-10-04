#define MAXDATA 10000
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>

#include "definition.h"
#include "utility.h"
#include "correction.h"

using namespace std;

//---------------------------------------------------
void recoHit(THeader *run)
//---------------------------------------------------
{
    TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
    TFile *fIn = new TFile (fName,"UPDATE");
    TTree *t = (TTree*) fIn->Get("dRICH");

    int evt,nedge, pmt[MAXDATA], pol[MAXDATA], slot[MAXDATA], fiber[MAXDATA], ch[MAXDATA];
    double x[MAXDATA], y[MAXDATA], r[MAXDATA], nt[MAXDATA], dur[MAXDATA], nttw[MAXDATA];
    t->SetBranchAddress("evt",&evt);
    t->SetBranchAddress("nedge",&nedge);
    t->SetBranchAddress("pol",&pol);
    t->SetBranchAddress("slot",&slot);
    t->SetBranchAddress("fiber",&fiber);
    t->SetBranchAddress("ch",&ch);
    t->SetBranchAddress("pmt",&pmt);
    t->SetBranchAddress("x",&x);
    t->SetBranchAddress("y",&y);
    t->SetBranchAddress("r",&r);
    t->SetBranchAddress("nt",&nt);

    bool goodHit[MAXDATA], externalPhoton[MAXDATA];;
    auto tgoodHit= t->Branch("goodHit",&goodHit,"goodHit[nedge]/O");
    auto tdur = t->Branch("dur",&dur,"dur[nedge]/D");
    auto tnttw = t->Branch("nttw",&nttw,"nttw[nedge]/D");
    auto tExternalPhoton= t->Branch("externalPhoton",&externalPhoton,"externalPhoton[nedge]/O");

    cout <<"Applying selection based on time and radius RMS\n";
    for(int i = 0; i < t->GetEntries(); i++){
        if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
        t->GetEntry(i);
        for(int j = 0; j < nedge; j++){
            goodHit[j]=false;
            nttw[j]=nt[j];
            dur[j]=0;
        }
        if(i<10){
            printf(" EDGE event %3d %6d\n",i,evt);
	for(int     j=0; j<nedge; j++)printf("e  %3d %3d (%4d %4d %4d) %7.2f \n",j,pol[j],slot[j],fiber[j],ch[j],nt[j]);
            printf(" HIT event %3d \n",i);
        }
        for(int j = 0; j < nedge; j++){
	    externalPhoton[j]=false;
            if(r[j] > run->radCut)externalPhoton[j]=true;
            //cout <<pol[j] <<endl;
            if(pol[j]!=0)continue;
            for(int k = 0; k < nedge; k++){
                if(pol[k]!=1)continue;
                //cout <<"Pol " <<pol[j] <<" " <<pol[k]<<endl;
                if(slot[j]!=slot[k])continue;
                if(fiber[j]!=fiber[k])continue;
                if(ch[j]!=ch[k])continue;
                //cout <<"Pmt " <<pmt[j] <<" " <<pmt[k]<<endl;
                if(nt[k] < nt[j])continue;
                //cout <<"Time " <<nt[j] <<" " <<nt[k] <<endl;
                if(nt[j] - nt[j] > run->MaxHitLength) continue;
                //cout <<"X " <<x[j] <<" " <<x[k] <<endl;
                //cout <<"Y " <<y[j] <<" " <<y[k] <<endl;
                if(x[j]==x[k] && y[j]==y[k]){
                    goodHit[j]=true;
                    goodHit[k]=true;
                    dur[j] = nt[k]-nt[j];
                    dur[k] = nt[k]-nt[j];
                    nttw[j] = timeWalkCorrection(nt[j],nt[k]);
		    externalPhoton[j]=false;
                    //if(r[j] > run->geoCut)externalPhoton[j]=true;
                    if(r[j] > run->radCut)externalPhoton[j]=true;
	            if(i<10)printf(" h %3d %3d  %7.2f %7.2f %7.2f \n",j,k,nt[j],nttw[j],dur[j]);
                    break;
                }
            }
        }
        tgoodHit->Fill();
        tdur->Fill();
        tnttw->Fill();
        tExternalPhoton->Fill();
    }
    printEnd();
    t->Write("",TObject::kOverwrite);
    fIn->Close();
}


//---------------------------------------------------
void rmsCutSelection(THeader *run){
//---------------------------------------------------
//
    TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
    TFile *fIn = new TFile (fName,"UPDATE");
    TTree *t = (TTree*) fIn->Get("dRICH");

    int nedge, pmt[MAXDATA];
    bool goodPhoton[MAXDATA],externalPhoton[MAXDATA];
    double rsdRadius[MAXDATA], rsdTime[MAXDATA];
    bool goodRMS[10];
    bool goodHit[MAXDATA];
    t->SetBranchAddress("nedge",&nedge);
    t->SetBranchAddress("pmt",&pmt);
    t->SetBranchAddress("goodPhoton",&goodPhoton);
    t->SetBranchAddress("externalPhoton",&externalPhoton);
    t->SetBranchAddress("rsdRadius",&rsdRadius);
    t->SetBranchAddress("rsdTime",&rsdTime);
    t->SetBranchAddress("goodRMS",&goodRMS);
    t->SetBranchAddress("goodHit",&goodHit);
    //t->SetBranchAddress("goodHit",&goodHit);

      float gxa, gya, gxtheta, gytheta;
      t->SetBranchAddress("gxa",&gxa);
      t->SetBranchAddress("gya",&gya);
      t->SetBranchAddress("gxtheta",&gxtheta);
      t->SetBranchAddress("gytheta",&gytheta);
    bool cutPhotonFlag[MAXDATA];
    auto tcutPhotonFlag= t->Branch("cutPhotonFlag",&cutPhotonFlag,"cutPhotonFlag[nedge]/O");

    cout <<"Applying selection based on time and radius RMS\n";
    for(int i = 0; i < t->GetEntries(); i++){
      if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
      t->GetEntry(i);
      //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
      //if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
      for(int j = 0; j < nedge; j++){
        cutPhotonFlag[j]=false;
        if(goodHit[j]==false)continue;
        int k=0;
        double cmpRadius = run->cutRadiusInRMS;
        double cmpTime = run->cutTimeInRMS;
        if(externalPhoton[j]==true) {
          k = 1;
          cmpRadius = run->cutRadiusOutRMS;
          cmpTime = run->cutTimeOutRMS;
        }
        int refTOT = 4+5*k;
        if(goodPhoton[j]==true && goodRMS[refTOT]==true){
          if(abs(rsdRadius[j]) < cmpRadius &&  abs(rsdTime[j]) < cmpTime)cutPhotonFlag[j] = true; 
        }
      } 
      tcutPhotonFlag->Fill();
    }
    printEnd();
    t->Write("",TObject::kOverwrite);
    fIn->Close();
}


//---------------------------------------------------
void selectGoodPhotons(THeader *run)
//---------------------------------------------------
{

      TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
      TFile *fIn = new TFile (fName,"UPDATE");
      TTree *t = (TTree*) fIn->Get("dRICH");

      int nedge, pol[MAXDATA],evt;
      double x[MAXDATA],y[MAXDATA],r[MAXDATA],nttw[MAXDATA], dur[MAXDATA];
      bool goodHit[MAXDATA],trigSig[MAXDATA],externalPhoton[MAXDATA];
      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("nedge",&nedge);
      t->SetBranchAddress("pol",&pol);
      t->SetBranchAddress("nttw",&nttw);
      t->SetBranchAddress("x",&x);
      t->SetBranchAddress("y",&y);
      t->SetBranchAddress("r",&r);
      t->SetBranchAddress("dur",&dur);
      t->SetBranchAddress("goodHit",&goodHit);
      t->SetBranchAddress("trigSig",&trigSig);
      t->SetBranchAddress("externalPhoton",&externalPhoton);

      float gxa, gya, gxtheta, gytheta;
      t->SetBranchAddress("gxa",&gxa);
      t->SetBranchAddress("gya",&gya);
      t->SetBranchAddress("gxtheta",&gxtheta);
      t->SetBranchAddress("gytheta",&gytheta);
      
      bool coincPhoton[MAXDATA],goodPhoton[MAXDATA];
      auto tCoincPhoton= t->Branch("coincPhoton",&coincPhoton,"coincPhoton[nedge]/O");
      auto tGoodPhoton= t->Branch("goodPhoton",&goodPhoton,"goodPhoton[nedge]/O");

      cout <<"Selecting the photon in the time coincidence window and duration\n";
      for(int i = 0; i < t->GetEntries(); i++){
          if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
          t->GetEntry(i);
	  if(i<10)printf("EVENT %3d %6d \n",i,evt);
	  //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
	  //if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
          for(int j = 0; j < nedge; j++){
              goodPhoton[j]=false;
              coincPhoton[j]=false;
              if(goodHit[j]==false)continue;

	      double timeMin = run->timeInMin;
	      double timeMax = run->timeInMax;
	      if(externalPhoton[j]==true){
	          timeMin = run->timeOuMin;
	          timeMax = run->timeOuMax;
	      }
              if(pol[j]==0 && trigSig[j]==false && nttw[j] > timeMin && nttw[j] < timeMax){ 
                  coincPhoton[j]=true;
                  if(dur[j] > run->durMin) goodPhoton[j]=true;
	      }
	      if(i<10 && pol[j]==0)printf(" sele pho %3d  time %7.2f (%7.2f : %7.2f) dur %7.2f  r %7.2f  %3d --> %3d %3d \n",
			     j,nttw[j],timeMin,timeMax,dur[j],r[j],externalPhoton[j],coincPhoton[j],goodPhoton[j]);
          }
          tCoincPhoton->Fill();
          tGoodPhoton->Fill();
      }
      printEnd();
      t->Write("",TObject::kOverwrite);
      fIn->Close();

}


//---------------------------------------------------
//void selectGoodPhotons(THeader *run)
//---------------------------------------------------
/*{

      TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
      TFile *fIn = new TFile (fName,"UPDATE");
      TTree *t = (TTree*) fIn->Get("dRICH");

      int nedge, pol[MAXDATA];
      double x[MAXDATA],y[MAXDATA],r[MAXDATA],nttw[MAXDATA], dur[MAXDATA];
      bool coincPhoton[MAXDATA],trigSig[MAXDATA];

      t->SetBranchAddress("nedge",&nedge);
      t->SetBranchAddress("pol",&pol);
      t->SetBranchAddress("nttw",&nttw);
      t->SetBranchAddress("x",&x);
      t->SetBranchAddress("y",&y);
      t->SetBranchAddress("r",&r);
      t->SetBranchAddress("dur",&dur);
      t->SetBranchAddress("coincPhoton",&coincPhoton);
      t->SetBranchAddress("trigSig",&trigSig);

      bool goodPhoton[MAXDATA], externalPhoton[MAXDATA];
      auto tGoodPhoton= t->Branch("goodPhoton",&goodPhoton,"goodPhoton[nedge]/O");
      //auto tExternalPhoton= t->Branch("externalPhoton",&externalPhoton,"externalPhoton[nedge]/O");

      cout <<"Selecting the photon in the time goodidence window and dividing rings\n";
      for(int i = 0; i < t->GetEntries(); i++){
          if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
          t->GetEntry(i);
          for(int j = 0; j < nedge; j++){
              goodPhoton[j]=false;
              //externalPhoton[j]=false;
              if(coincPhoton[j]==false)continue;
              if(dur[j] > run->durMin) goodPhoton[j]=true;
              //if(r[j] > run->geoCut)externalPhoton[j]=true;
              //if(r[j] > run->radCut)externalPhoton[j]=true;
              //Attention, if you refine the cut you should check the conversion to milliradiant.
          }
          tGoodPhoton->Fill();
          //tExternalPhoton->Fill();
      }
      printEnd();
      t->Write("",TObject::kOverwrite);
      fIn->Close();

}*/


//---------------------------------------------------
void findTimeCoincidence(THeader *run)
//---------------------------------------------------
{

    TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
    TFile *fIn = new TFile (fName,"READ");
    TTree *t = (TTree*) fIn->Get("dRICH");

    int nedge, pol[MAXDATA];
    bool externalPhoton[MAXDATA], goodHit[MAXDATA], trigSig[MAXDATA];
    double nttw[MAXDATA], dur[MAXDATA];
    t->SetBranchAddress("nedge",&nedge);
    t->SetBranchAddress("pol",&pol);
    t->SetBranchAddress("dur",&dur);
    t->SetBranchAddress("goodHit",&goodHit);
    t->SetBranchAddress("externalPhoton",&externalPhoton);
    t->SetBranchAddress("trigSig",&trigSig);
    t->SetBranchAddress("nttw",&nttw);

    TH1D *hint = new TH1D("hint","hint",2001,-0.5,2000.5);
    TH1D *hout = new TH1D("hout","hout",2001,-0.5,2000.5);

    cout <<"Finding the time-coincidence window\n";
    for(int i = 0; i < t->GetEntries(); i++){
        if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
        t->GetEntry(i);
        for(int j = 0; j < nedge; j++){
            if(pol[j]==1) continue;
            if(goodHit[j]==false)continue;
            if(trigSig[j]==true)continue;
            if(dur[j]<CUT_MIN_DUR)continue;
            if(externalPhoton[j]==false) hint->Fill(nttw[j]);
            if(externalPhoton[j]==true) hout->Fill(nttw[j]);
        }
    }
    if(SHOW_PROGRESS==true)printEnd();
    
    TCanvas *c0 = new TCanvas("c0","c0",1600,900);
    c0->Draw();
    TF1 *fint0 = new TF1("fint0","pol0",0,2000);
    hint->Fit(fint0,"Q","",1800,2000);
    TF1 *fint = new TF1("fint","gaus(0)+pol0(3)",0,2000);
    fint->SetParameters(10000,hint->GetBinCenter(hint->GetMaximumBin()),5,fint0->GetParameter(0));
    fint->SetParLimits(2,0.5,30);
    fint->FixParameter(3,fint0->GetParameter(0));
    double timeFitInMin=hint->GetBinCenter(hint->GetMaximumBin())-CUT_FIT_COINC;
    double timeFitInMax=hint->GetBinCenter(hint->GetMaximumBin())+CUT_FIT_COINC;
    hint->Fit("fint","Q","",timeFitInMin,timeFitInMax);
    hint->Draw();
    run->timeInMin=fint->GetParameter(1)-CUT_NSIGMA_TIME*fint->GetParameter(2); 
    run->timeInMax=fint->GetParameter(1)+CUT_NSIGMA_TIME*fint->GetParameter(2);
    c0->Print("cIn.root");

    TF1 *fout0 = new TF1("fout0","pol0",0,2000);
    hout->Fit(fout0,"Q","",1800,2000);
    TF1 *fout = new TF1("fout","gaus(0)+pol0(3)",0,2000);
    fout->SetParameters(10000,hout->GetBinCenter(hout->GetMaximumBin()),5,fout0->GetParameter(0));
    fout->SetParLimits(2,0.5,30);
    fout->FixParameter(3,fout0->GetParameter(0));
    double timeFitOuMin=hout->GetBinCenter(hout->GetMaximumBin())-CUT_FIT_COINC;
    double timeFitOuMax=hout->GetBinCenter(hout->GetMaximumBin())+CUT_FIT_COINC;
    hout->Fit("fout","Q","",timeFitOuMin,timeFitOuMax);
    hout->Draw();
    run->timeOuMin=fout->GetParameter(1)-CUT_NSIGMA_TIME*fout->GetParameter(2); 
    run->timeOuMax=fout->GetParameter(1)+CUT_NSIGMA_TIME*fout->GetParameter(2);
    c0->Print("cOut.root");

    printf(" Event Time window In: %7.2f %7.2f  Ou: %7.2f %7.2f \n",
            run->timeInMin,run->timeInMax,run->timeOuMin,run->timeOuMax);

    fIn->Close();

}
