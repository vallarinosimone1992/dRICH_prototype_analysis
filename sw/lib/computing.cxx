#define MAXDATA 10000
#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>

#include "computing.h"
#include "definition.h"
#include "utility.h"

using namespace std;

//---------------------------------------------------
double mmTomRad(double r, double inPath, double zMir){
//---------------------------------------------------
  double path=0;
  if(inPath > 1000){
    path=inPath-zMir;
  }else{
    path=inPath+zMir;
  }
  return 1000*atan(r/path);
}


//---------------------------------------------------
//void convertToRadiant(THeader *run)
//---------------------------------------------------
/*{
    TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
    TFile *fIn = new TFile (fName,"UPDATE");
    TTree *t = (TTree*) fIn->Get("dRICH");
    int nedge;
    double nr[MAXDATA];
    bool cutPhotonFlag[MAXDATA],externalPhoton[MAXDATA];
    t->SetBranchAddress("nedge",&nedge);
    t->SetBranchAddress("nr",&nr);
    t->SetBranchAddress("externalPhoton",&externalPhoton);

    double rRad[MAXDATA];
    auto trRad=t->Branch("rRad",&rRad,"&rRad[nedge]/D");

    cout <<"Convert radius to radiant\n";
    for(int i = 0; i < t->GetEntries(); i++){
        if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
        t->GetEntry(i);
        for(int j = 0; j < nedge; j++){
            double path=0;
            if(externalPhoton[j]==true){
                path = run->firstPath - run->firstMirrorPosition;
            }else{
                path = run->secondPath - run->secondMirrorPosition;
            }
            rRad[j] = atan(nr[j]/path);
        }
        trRad->Fill();
    }
    printEnd();
    t->Write("",TObject::kOverwrite);
    fIn->Close();

}*/

void fitSigma(TH1D *h, bool out){
  TF1 *f = new TF1("f","sqrt([0]*[0]/x + [1]*[1])",0,100);
  f->SetParameters(1.1,0.1);
  int min=0, max=0;
  if(out == false){
    min=6;
    max=22;
  }else{
    min=2;
    max=10;
  }
  h->Fit("f","Q","",min,max);
}

TF1* getFun(TH1D *h, bool out){
  TF1 *f;
  if(out == false){
    h->GetXaxis()->SetRangeUser(25,50);
  }else{
    h->GetXaxis()->SetRangeUser(130,220);
  }
  f = new TF1("f","gaus(0)",0,250);
  double p1= h->GetBinCenter(h->GetMaximumBin());
  f->SetParameters(5000,p1,2);
  h->Fit(f->GetName(),"Q","",p1-2,p1+3);
  return  f;
}



double getSigma(TH1D *h, TF1 *f,string fname, bool out){
  if(out == false){
    h->GetXaxis()->SetRangeUser(10,60);
  }else{
    h->GetXaxis()->SetRangeUser(150,220);
  }
  f = new TF1(fname.c_str(),"gaus(0)",0,250);
  double p1= h->GetBinCenter(h->GetMaximumBin());
  f->SetParameters(5000,p1,2);
  //h->Fit(f->GetName(),"Q","",p1-1,p1+1.5);
  h->Fit(f->GetName(),"Q","",p1-2,p1+4);
  return  f->GetParameter(2);
}


void applyFit(TH1D *h, TF1 *f,string fname, bool out){
  //double fitCenter, fitMin, fitMax;
  if(out == false){
    double min = h->GetBinCenter(h->GetMaximumBin())-4;
    double max = h->GetBinCenter(h->GetMaximumBin())+4;
    h->GetXaxis()->SetRangeUser(min,max);
    //h->GetXaxis()->SetRangeUser(10,55);
  }else{
    h->GetXaxis()->SetRangeUser(150,220);
  }
  f = new TF1(fname.c_str(),"gaus(0)",0,220);
  double p1= h->GetBinCenter(h->GetMaximumBin());
  f->SetParameters(5000,p1,2);
  //h->Fit(f->GetName(),"Q","",p1-1,p1+1.5);
  h->Fit(f->GetName(),"Q","",p1-2,p1+4);
  h->Draw();
}


//---------------------------------------------------
void computeRMS(THeader *run, int recompute=-1)
//---------------------------------------------------
{
    TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
    TFile *fIn = new TFile (fName,"UPDATE");
    TTree *t = (TTree*) fIn->Get("dRICH");

    int nedge, pmt[MAXDATA], time[MAXDATA], pol[MAXDATA], evt;
    double r[MAXDATA], spRadius[10], spTime[10];
    double nr[MAXDATA], nttw[MAXDATA], spnRadius[10], spnTime[10];
    bool coincPhoton[MAXDATA],externalPhoton[MAXDATA], goodPhoton[MAXDATA];
    t->SetBranchAddress("evt",&evt);
    t->SetBranchAddress("nedge",&nedge);
    t->SetBranchAddress("pmt",&pmt);
    t->SetBranchAddress("pol",&pol);
    t->SetBranchAddress("r",&r);
    t->SetBranchAddress("time",&time);
    t->SetBranchAddress("nr",&nr);
    t->SetBranchAddress("nttw",&nttw);
    t->SetBranchAddress("coincPhoton",&coincPhoton);
    t->SetBranchAddress("externalPhoton",&externalPhoton);
    t->SetBranchAddress("goodPhoton",&goodPhoton);
    t->SetBranchAddress("spRadius",&spRadius);
    t->SetBranchAddress("spTime",&spTime);
    t->SetBranchAddress("spnRadius",&spnRadius);
    t->SetBranchAddress("spnTime",&spnTime);

    float gxa, gya, gxtheta, gytheta;
    t->SetBranchAddress("gxa",&gxa);
    t->SetBranchAddress("gya",&gya);
    t->SetBranchAddress("gxtheta",&gxtheta);
    t->SetBranchAddress("gytheta",&gytheta);
    
    double rmsRadius[10], rmsTime[10], rmsPhoton[10];
    double rsdRadius[MAXDATA], rsdTime[MAXDATA];
    bool goodRMS[10];
    auto trmsRadius=t->Branch("rmsRadius",&rmsRadius,"rmsRadius[10]/D");
    auto trmsPhoton=t->Branch("rmsPhoton",&rmsPhoton,"rmsPhoton[10]/I");
    auto trmsTime=t->Branch("rmsTime",&rmsTime,"rmsTime[10]/D");
    auto tgoodRMS=t->Branch("goodRMS",&goodRMS,"goodRMS[10]/O");
    auto trsdRadius=t->Branch("rsdRadius",&rsdRadius,"rsdRadius[nedge]/D");
    auto trsdTime=t->Branch("rsdTime",&rsdTime,"rsdTime[nedge]/D");

    cout <<"Computing RMS and residui.\n";
    for(int i = 0; i < t->GetEntries(); i++){
        if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
        t->GetEntry(i);
        for(int j = 0; j < 10; j++){
            rmsRadius[j]=0;
            rmsPhoton[j]=0;
            rmsTime[j]=0;
            goodRMS[j]=false;
        }
    	//if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
    	//if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
        if(i<10)printf(" Computing eve %3d %6d \n",i,evt);
        bool countIt[MAXDATA]={false};
        if(recompute == -1){
            for(int j = 0; j < nedge; j++){
		if(goodPhoton[j]==false)continue;

                double   tmpRMSrIn=0, tmpRMStIn=0, rMeanIn=0, tMeanIn=0;
                double   tmpRMSrOut=0, tmpRMStOut=0, rMeanOut=0, tMeanOut=0;
                int   tmpCountIn=0, tmpCountOut=0;
                int   countIn=0, countOut=0;
		for(int k = 0; k < nedge; k++){
                    if(k==j)continue;
		    if(goodPhoton[k]==false)continue;
		    if(externalPhoton[k]==false){
		    	rMeanIn+=nr[k];
			tMeanIn+=nttw[k];
			countIn++;
		    }else{
		    	rMeanOut+=nr[k];
			tMeanOut+=nttw[k];
			countOut++;
		    }
		}
		if(countIn>0){
		    rMeanIn/=countIn;
		    tMeanIn/=countIn;
		}
		if(countOut>0){
		    rMeanOut/=countOut;
		    tMeanOut/=countOut;
		}
                for(int k = 0; k < nedge; k++){
                    if(k==j)continue;
                    if(goodPhoton[k]==true && externalPhoton[k]==false){
                        tmpRMSrIn+=pow(nr[k]-rMeanIn,2);
                        tmpRMStIn+=pow(nttw[k]-tMeanIn,2);
                        tmpCountIn++;
                    }
                    if(goodPhoton[k]==true && externalPhoton[k]==true){
                        tmpRMSrOut+=pow(nr[k]-rMeanOut,2);
                        tmpRMStOut+=pow(nttw[k]-tMeanOut,2);
                        tmpCountOut++;
                    }
                }
		if(tmpCountIn>0){
                    tmpRMSrIn=sqrt(tmpRMSrIn/tmpCountIn);
                    tmpRMStIn=sqrt(tmpRMStIn/tmpCountIn);
		    if(tmpRMSrIn<RMS_ANGLE_GAS)tmpRMSrIn=RMS_ANGLE_GAS;
		    if(tmpRMSrIn<RMS_TIME_GAS)tmpRMSrIn=RMS_TIME_GAS;
		}
		if(tmpCountOut>0){
                    tmpRMSrOut=sqrt(tmpRMSrOut/tmpCountOut);
                    tmpRMStOut=sqrt(tmpRMStOut/tmpCountOut);
		    if(tmpRMSrIn<RMS_ANGLE_AER)tmpRMSrIn=RMS_ANGLE_AER;
		    if(tmpRMSrIn<RMS_TIME_AER)tmpRMSrIn=RMS_TIME_AER;
		}
                if(externalPhoton[j]==true   && tmpCountOut >= 2){
                    double rsdR=nr[j]-rMeanOut;
                    double rsdT=nttw[j]-tMeanOut;
                    if(abs(rsdR) < 3*tmpRMSrOut && abs(rsdT) < 3*tmpRMStOut) countIt[j]=true;
                }else   if(externalPhoton[j]==false && tmpCountIn >= 2){
                    double rsdR=nr[j]-rMeanIn;
                    double rsdT=nttw[j]-tMeanIn;
                    if(abs(rsdR) < 3*tmpRMSrIn && abs(rsdT) < 3*tmpRMStIn) countIt[j]=true;
                }
		if(i<10 && externalPhoton[j]==true)
                    printf(" e %3d r  %7.2f vs %7.2f %7.2f  (%7.2f x3) t %7.2f vs %7.2f %7.2f  (%7.2f x3) --> %3d\n",
                        j,nr[j],rMeanOut,nr[j]-rMeanOut,tmpRMSrOut,nttw[j],tMeanOut,nttw[j]-tMeanOut,tmpRMStOut,countIt[j]);
            }
      }else{
          for(int j = 0; j < nedge; j++){
          countIt[j]=true;
        }
      }

      for(int j = 0; j < nedge; j++){
        //cout <<goodPhoton[j] <<" " <<countIt[j];
        //cin.get();
        if(countIt[j]==false)continue;
        int k=0;
        if(externalPhoton[j]==true) k = 1;
        int refPMT = pmt[j]+5*k;
        int refTOT = 4+5*k;
        rsdRadius[j]=-1000;
        rsdTime[j]=-1000;
        if(goodPhoton[j]==true){
          rsdRadius[j] = nr[j]-spnRadius[refTOT];
          rsdTime[j]   = nttw[j]-spnTime[refTOT];
          rmsRadius[refTOT]+=pow(rsdRadius[j],2);
          rmsTime[refTOT]+=pow(rsdTime[j],2);
          rmsPhoton[refTOT]+=1;

          rmsRadius[refPMT]+=pow(spnRadius[refPMT]-nr[j],2);
          rmsTime[refPMT]+=pow(spnTime[refPMT]-nttw[j],2);
          rmsPhoton[refPMT]+=1;

        }
      }
      for(int j = 0; j < 10; j++){
          if(rmsPhoton[j]>1){
              rmsRadius[j]=sqrt(rmsRadius[j]/rmsPhoton[j]);
              rmsTime[j]=sqrt(rmsTime[j]/rmsPhoton[j]);
              goodRMS[j]=true;
          }else{
              rmsRadius[j]=-10;
              rmsTime[j]=-10;
          }
      }
      trmsRadius->Fill();
      trmsPhoton->Fill();
      trmsTime->Fill();
      tgoodRMS->Fill();
      trsdRadius->Fill();
      trsdTime->Fill();
    }
    printEnd();
    t->Write("",TObject::kOverwrite);
    fIn->Close();
}


void computeCutSingleParticle(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nttw[MAXDATA];
  bool cutPhotonFlag[MAXDATA],externalPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nttw",&nttw);
  t->SetBranchAddress("externalPhoton",&externalPhoton);
  t->SetBranchAddress("cutPhotonFlag",&cutPhotonFlag);

  float gxa, gya, gxtheta, gytheta;
  t->SetBranchAddress("gxa",&gxa);
  t->SetBranchAddress("gya",&gya);
  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);

  double cutRadius[10], cutTime[10];
  int cutPhoton[10];
  bool goodCUT[10];
  auto tcutRadius= t->Branch("cutRadius",&cutRadius,"cutRadius[10]/D");
  auto tcutTime= t->Branch("cutTime",&cutTime,"cutTime[10]/D");
  auto tcutPhoton= t->Branch("cutPhoton",&cutPhoton,"cutPhoton[10]/I");
  auto tgoodCUT=t->Branch("goodCUT",&goodCUT,"goodCUT[10]/O");

  cout <<"Computing single particles property after cut based on rms.\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      cutRadius[j]=0;
      cutTime[j]=0;
      cutPhoton[j]=0;
    }
    //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
    //if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
    for(int j = 0; j < nedge; j++){
      int k=0;
      if(externalPhoton[j]==true) k = 1;
      int refPMT = pmt[j]+5*k;
      int refTOT = 4+5*k;
      if(cutPhotonFlag[j]==true){
        cutRadius[refPMT]+=nr[j];
        cutPhoton[refPMT]+=1;
        cutTime[refPMT]+=nttw[j];
        cutRadius[refTOT]+=nr[j];
        cutPhoton[refTOT]+=1;
        cutTime[refTOT]+=nttw[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(cutPhoton[j]!=0){
        cutRadius[j]/=cutPhoton[j];
        cutTime[j]/=cutPhoton[j];
        goodCUT[j]=true;
        //cout <<j<<" " <<cutRadius[j]<<" " <<cutTime[j] <<" "<<cutPhoton[j] <<endl;
      }else{
        cutRadius[j]=0;
        cutTime[j]=0;
      }
    }
    tcutRadius->Fill();
    tcutTime->Fill();
    tcutPhoton->Fill();
    tgoodCUT->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

//---------------------------------------------------
void newSingleParticle(THeader *run){
//---------------------------------------------------
//
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
  TFile *fIn = new TFile(fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nttw[MAXDATA];
  bool goodPhoton[MAXDATA],externalPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nttw",&nttw);
  t->SetBranchAddress("goodPhoton",&goodPhoton);
  t->SetBranchAddress("externalPhoton",&externalPhoton);

  float gxa, gya, gxtheta, gytheta;
  t->SetBranchAddress("gxa",&gxa);
  t->SetBranchAddress("gya",&gya);
  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  
  double spnRadius[10], spnTime[10];
  int spnPhoton[10];
  bool goodSPN[10];

  auto tSpnRadius=t->Branch("spnRadius",&spnRadius,"spnRadius[10]/D");
  auto tSpnPhoton=t->Branch("spnPhoton",&spnPhoton,"spnPhoton[10]/I");
  auto tSpnTime=t->Branch("spnTime",&spnTime,"spnTime[10]/D");
  auto tgoodSPN=t->Branch("goodSPN",&goodSPN,"goodSPN[10]/O");

  cout <<"Computing mean quantities for each single particle\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      spnRadius[j]=0;
      spnPhoton[j]=0;
      spnTime[j]=0;
      goodSPN[j]=false;
    }
    //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
    //if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
    for(int j = 0; j < nedge; j++){
      int k=0;
      if(externalPhoton[j]==true) k = 1;
      int refPMT = pmt[j]+5*k;
      int refTOT = 4+5*k;
      if(goodPhoton[j]==true){
        spnRadius[refPMT]+=nr[j];
        spnPhoton[refPMT]+=1;
        spnTime[refPMT]+=nttw[j];
        spnRadius[refTOT]+=nr[j];
        spnPhoton[refTOT]+=1;
        spnTime[refTOT]+=nttw[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(spnPhoton[j]!=0){
        spnRadius[j]/=spnPhoton[j];
        spnTime[j]/=spnPhoton[j];
        goodSPN[j]=true;
        //cout <<j<<" " <<spnRadius[j]<<" " <<spnTime[j] <<" "<<spnPhoton[j] <<endl;
      }else{
        spnRadius[j]=0;
        spnTime[j]=0;
      }
    }
    tSpnRadius->Fill();
    tSpnPhoton->Fill();
    tSpnTime->Fill();
    tgoodSPN->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}


//---------------------------------------------------
void singleParticle(THeader *run)
//---------------------------------------------------
{
    TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
    TFile *fIn = new TFile (fName,"UPDATE");
    TTree *t = (TTree*) fIn->Get("dRICH");

    int nedge, pol[MAXDATA], pmt[MAXDATA];
    double x[MAXDATA],y[MAXDATA],r[MAXDATA], nttw[MAXDATA], rmm[MAXDATA];
    bool goodPhoton[MAXDATA], externalPhoton[MAXDATA];

    t->SetBranchAddress("nedge",&nedge);
    t->SetBranchAddress("pol",&pol);
    t->SetBranchAddress("nttw",&nttw);
    t->SetBranchAddress("pmt",&pmt);
    t->SetBranchAddress("x",&x);
    t->SetBranchAddress("y",&y);
    t->SetBranchAddress("r",&r);
    t->SetBranchAddress("rmm",&rmm);
    t->SetBranchAddress("goodPhoton",&goodPhoton);
    t->SetBranchAddress("externalPhoton",&externalPhoton);

    float gxa, gya, gxtheta, gytheta;
    t->SetBranchAddress("gxa",&gxa);
    t->SetBranchAddress("gya",&gya);
    t->SetBranchAddress("gxtheta",&gxtheta);
    t->SetBranchAddress("gytheta",&gytheta);

    double spRadius[10], spTime[10], spRadiusmm[10];
    int spPhoton[10];
    bool goodSP[10];
    
    auto tSPRadius=t->Branch("spRadius",&spRadius,"spRadius[10]/D");
    auto tSPRadiusmm=t->Branch("spRadiusmm",&spRadiusmm,"spRadiusmm[10]/D");
    auto tSPPhoton=t->Branch("spPhoton",&spPhoton,"spPhoton[10]/I");
    auto tSPTime=t->Branch("spTime",&spTime,"spTime[10]/D");
    auto tgoodSP=t->Branch("goodSP",&goodSP,"goodSP[10]/O");

    cout <<"Selecting the photon\n";
    for(int i = 0; i < t->GetEntries(); i++){
        if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
        t->GetEntry(i);
        for(int j = 0; j < 10; j++){
            spRadius[j]=0;
            spRadiusmm[j]=0;
            spPhoton[j]=0;
            spTime[j]=0;
            goodSP[j]=false;
        }
        //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
        //if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
        for(int j = 0; j < nedge; j++){
            if(goodPhoton[j]==true){
                int k=0;
                if(externalPhoton[j]==true) k = 1;
                int refPMT = pmt[j]+5*k;
                int refTOT = 4+5*k;
                spRadius[refPMT]+=r[j];
                spRadiusmm[refPMT]+=rmm[j];
                spPhoton[refPMT]+=1;
                spTime[refPMT]+=nttw[j];
                spRadius[refTOT]+=r[j];
                spRadiusmm[refTOT]+=rmm[j];
                spPhoton[refTOT]+=1;
                spTime[refTOT]+=nttw[j];
            }
        }
        for(int j = 0; j < 10; j++){
            if(spPhoton[j]!=0){
                spRadius[j]/=spPhoton[j];
                spRadiusmm[j]/=spPhoton[j];
                spTime[j]/=spPhoton[j];
                goodSP[j]=true;
                //cout <<j<<" " <<spRadius[j]<<" " <<spTime[j] <<" "<<spPhoton[j] <<endl;
            }else{
                spRadius[j]=0;
                spRadiusmm[j]=0;
                spTime[j]=0;
            }
        }
        tSPRadius->Fill();
        tSPRadiusmm->Fill();
        tSPPhoton->Fill();
        tSPTime->Fill();
        tgoodSP->Fill();
    }
    printEnd();
    t->Write("",TObject::kOverwrite);
    fIn->Close();
}

