
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TChain.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooBreitWigner.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooBifurGauss.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include <fstream>
#define NCPU 4
using namespace RooFit;

Double_t thefit_real_scale(TString plotpath, Double_t *vect,  int doall = 3, Double_t LB =2.2, Double_t RB = 2.4, Double_t lm12=1.5, Double_t rm12=5, Double_t lm23=0, Double_t rm23=1.4, TString BINX = "1", TString BINY = "1"){
//Double_t thefit_real_scale(int doall = 3, Double_t LB =2.2, Double_t RB = 2.4, Double_t lm12=1.5, Double_t rm12=5, Double_t lm23=0, Double_t rm23=1.4, TString BINX = "1", TString BINY = "1"){
  //TString STR;
  //STR.Form("%i",STR_i);
  /* ////debug
  TString plotpath ="main";
  Double_t vect[8];
  lm12=3;
  rm12=3.5;
  lm23=0.5;
  rm23=1;
  //end debug*/
  TString MCPATH ="/home/mberger/hbk_work/pkpi/";
  TCanvas *c = new TCanvas("name", "name");
  RooRealVar m23("4m23","M(K#pi)", lm23, rm23);
  RooRealVar m12("4m12","M(pK)", lm12, rm12);
  RooRealVar m("m","M(p^{+}K^{-}#pi^{+})",2.2,2.4);
  TString PATH ="/home/mberger/hbk/tofit/pkpi/real";
  TString FILE = PATH  +"/calcmass.root";
  TString STR2;
   TString PLOT1 = "plots/"+plotpath+"/signal_"+BINX+"_"+BINY+".pdf";
   TString PLOT2 = "plots/"+plotpath+"/bkg_"+BINX+"_"+BINY+".pdf";
   TString PLOT3 = "plots/"+plotpath+"/total_"+BINX+"_"+BINY+".pdf";
  TString PARA1;
  TString PARA2;
  TString PARA3;
  TString PARA5;

  
  if (doall) PARA1 = "txt/"+plotpath+"/signal"+BINX+"_"+BINY+".txt";
  else PARA1 = "txt/signal.txt";
  //PARA2 = "txt/"+plotpath+"/sigcor"+BINX+"_"+BINY+".txt";
  
  if(doall>2) PARA2 = "txt/"+plotpath+"/bkg"+BINX+"_"+BINY+".txt";
  else PARA2 = "txt/bkg.txt";
    PARA3 = "txt/"+plotpath+"/yield_"+BINX+"_"+BINY+".txt";
    // PARA1 = "txt/signal.txt";
    // PARA2 = "txt/sigcor.txt";
  
  // PARA2 = "txt/bkg.txt";
  // PARA3 = "txt/yield.txt";

 
  TChain chain("hm");
  chain.Add(FILE);

   
  
  
    
    
    
    ///////////
  RooRealVar meansig("meansig","mean of signal",2.28594,2.265 ,2.35) ;
  RooRealVar sigmagau("sigmagau","width of gaussian signal",0.00345503,0.,0.2);
  RooRealVar sigmabwig("sigmabwig","width of bwign",0.0058709,0.,0.2) ;
  RooRealVar sigmabwig2("sigmabwig2","width of bwign2",345,0.,0.2) ;
  RooRealVar sigmagau2("sigmagau2","width of gaussian signal2",0.00562663,0.,0.2);
  RooRealVar sigmagau3("sigmagau3","width of gaussian signal3",0.0175976,0.,0.2);
  RooRealVar fsig ("fsig","gauss/bwign",0.6285,0.,1.);
  RooRealVar fsig2 ("fsig2","2",0.257487,0.,1.);
  //RooRealVar m("m","M(#Sigma#pi#pi)",2.2,2.4);
  RooGaussian gausssig("gausssig","gaussian PDF signal",m,meansig,sigmagau) ;
  RooGaussian gausssig2("gausssig2","gaussian PDF signal2",m,meansig,sigmagau2) ;
  //RooGaussian gausssig3("gausssig3","gaussian PDF signal3",m,meansig,sigmagau3) ;
  //    RooBifurGauss bigau("bigau","bigau",m,meansig,sigmagau2,sigmagau3) ;
  RooBreitWigner bwign("bwign","bwign",m,meansig,sigmabwig) ;
  //RooBreitWigner bwign2("bwign2","bwign2",m,meansig,sigmabwig2) ;
  //RooBifurGauss bigau("bigau","bigau",m,meansig,sigmagau2,sigmagau3) ;
  //RooAddPdf signal("signal","adpdf",RooArgList(gausssig,bwign,gausssig2),RooArgList(fsig,fsig2));
  RooAddPdf signal("signal","adpdf",RooArgList(gausssig,gausssig2,bwign),RooArgList(fsig,fsig2));
        
      
  RooRealVar coef1("coef1","coef1",-0.0987995,-1.,1.);
  RooRealVar coef2("coef2","coef2",-0.000122793,-1.,1.);
  RooRealVar coef3("coef3","coef3",0.00546435,-1.,1.);
   
        
  RooChebychev bkg("bkg","Background",m,RooArgSet(coef1,coef2,coef3)) ;
    //  m.setRange(2.26,2.31);
    RooPlot* frame = m.frame(Bins(300));
     if (doall){
    TChain chain1("hm");
    for (Int_t ll = 1; ll<3; ll++){
    
	STR2.Form("%i",ll);
	FILE = MCPATH + STR2 +"/signal_mfit.root";
	chain1.Add(FILE);
      }
    
 RooDataSet ds1("ds1","ds1",RooArgSet(m,m23,m12),Import(chain1)) ;
    
    ds1.plotOn(frame);
  
    signal.fitTo(ds1,NumCPU(NCPU));
    signal.plotOn(frame);
     frame->Draw();
    c->Print(PLOT1);
    RooArgSet(meansig,fsig,fsig2,sigmagau,sigmagau2,sigmabwig).writeToFile(PARA1);
     }
     if (!doall) RooArgSet(meansig,fsig,fsig2,sigmagau,sigmagau2,sigmabwig).readFromFile(PARA1);

       if (doall>2){
	m.setRange(LB,RB);  
    TChain chain3("hm");
    for (Int_t ll = 1; ll<3; ll++){

	STR2.Form("%i",ll);
	FILE = MCPATH + STR2 +"/bkg_mfit.root";
	chain3.Add(FILE);
      
    }   
    RooDataSet ds3("ds3","ds3",RooArgSet(m,m23,m12),Import(chain3)) ;
    frame = m.frame(Bins(300));
    ds3.plotOn(frame);

    bkg.fitTo(ds3,NumCPU(NCPU));
    bkg.plotOn(frame);
    // chi3=p->chiSquare(3);
    RooArgSet(coef1,coef2,coef3).writeToFile(PARA2);
    frame->Draw();
    c->Print(PLOT2);
    }
       if (doall<3)  RooArgSet(coef1,coef2,coef3).readFromFile(PARA2);
     
    RooRealVar nsig("numsig","number signal", 5000, 0, 2000000 );
    RooRealVar nbkg("numbkg","number background",10000 , 0,4000000 );
    ////set range
    // m.setRange("range",2.266,2.306);
    // RooExtendPdf esig("esig","extended signal p.d.f",signal,nsig,"range");
    //RooExtendPdf ebkg("ebkg","extended background p.d.f",bkg,nbkg,"range") ;

    RooRealVar scale("scalefactor","scalefactor", 1.07, 0.,2. );
   
    RooGenericPdf sigmagau1_real("sigma1_real","sigmag1","scalefactor*sigmagau",RooArgSet(scale,sigmagau));
    RooGenericPdf sigmagau2_real("sigma2_real","sigmag2","scalefactor*sigmagau2",RooArgSet(scale,sigmagau2));
    RooGenericPdf sigmabw_real("sigmabw_real","sigmabw","scalefactor*sigmabwig",RooArgSet(scale,sigmabwig));
    RooGaussian gausssig_real("gausssig_real","gaussian PDF signal",m,meansig,sigmagau1_real) ;
    RooGaussian gausssig2_real("gausssig2_real","gaussian PDF signal2",m,meansig,sigmagau2_real) ;
    RooBreitWigner bwign_real("bwign_real","bwign",m,meansig,sigmabw_real) ;

    RooAddPdf signal_real("signal_real","adpdf",RooArgList(gausssig_real,gausssig2_real,bwign_real),RooArgList(fsig,fsig2));
    
    RooAddPdf etotal("etotal","etotal",RooArgSet(signal_real,bkg),RooArgSet(nsig,nbkg));
    //RooAddPdf etotal("etotal","etotal",RooArgSet(signal,bkg),RooArgSet(nsig,nbkg));
    //RooAddPdf etotal("etotal","etotal",RooArgList(esig,ebkg));
    //meansig.setConstant(kTRUE);
    sigmabwig.setConstant(kTRUE);
    sigmagau.setConstant(kTRUE);
    sigmagau2.setConstant(kTRUE);
    sigmagau3.setConstant(kTRUE);
    fsig.setConstant(kTRUE);
    fsig2.setConstant(kTRUE);
    coef1.setConstant(kTRUE);
    coef2.setConstant(kTRUE);
    coef3.setConstant(kTRUE);
  

    //sigmabwig_cor.setConstant(kTRUE);
    //m.setRange("sigdal",2.2806,2.2926);
    // m.setRange(2.2,2.4); 
    
    RooDataSet ds("ds","ds",RooArgSet(m,m23,m12),Import(chain)) ;
    frame = m.frame(Bins(300));
    ds.plotOn(frame);
    frame->Draw();

    RooFitResult* res = etotal.fitTo(ds,Save(),NumCPU(NCPU));
    //fitresult->Print();
    RooPlot *p = etotal.plotOn(frame);
    Double_t chi2 = p->chiSquare(11);
    etotal.plotOn(frame,Components("signal_real"),LineStyle(kDashed),LineColor(kRed));
   
    etotal.plotOn(frame,Components("bkg"),LineStyle(kDashed),LineColor(kGreen));
    
    std::cout<<chi2<<std::endl;
    //cb.plotOn(frame,Components("gauss2"),LineStyle(kDashed),LineColor(kGreen));
    //testg.plotON(frame,LineStyle(kDashed),LineColor);
       frame->Draw();
       /*ds.Print();
    
    nsig.Print();
    nbkg.Print();
       */
       RooArgSet(nsig,nbkg,scale).writeToFile(PARA3);
       vect[0]=nsig.getValV();
       vect[1]=nsig.getError();
       vect[2]=nbkg.getValV();
       vect[3]=nbkg.getError();
       vect[4]=scale.getValV();
       vect[5]=scale.getError();
       vect[6]=ds.sumEntries();
       c->Print(PLOT3);
    res->Print();
    return chi2;
}
