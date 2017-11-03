// combined fit of two histogram with separate functions but a common parameter
// N.B. this macro must be compiled with ACliC 

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"


Double_t funcBg(Double_t *x, Double_t *par){
    
    // fit of the region of roi, nickel, and cu ka kb
    // this is for fitting an already scaled histogram
    
    
    Double_t xx = x[0];
    
    //par[0] = background constant
    //par[10] = background slope
    
    Double_t back = par[0] + (xx - 7000) * par[10];
    
    //par[1] = cu ka1 gain
    //par[2] = cu ka1 mean
    //par[3] = cu ka1 sigma
    
    Double_t cuKa1 = par[1]/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-par[2])*(xx-par[2]))/(2*par[3]*par[3]));
    
    Double_t cuKa2Gain = par[1] * 0.51;
    Double_t cuKa2Mean = par[2] - 19.95;
    
    Double_t cuKa2 = cuKa2Gain/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-cuKa2Mean)*(xx-cuKa2Mean))/(2*par[3]*par[3]));
    
    //par[4] = cu kb gain
    //par[5] = cu kb mean
    //par[6] = cu kb sigma
    
    Double_t cuKb = par[4]/(sqrt(2*TMath::Pi())*par[6])*TMath::Exp(-((xx-par[5])*(xx-par[5]))/(2*par[6]*par[6]));
    
    //par[7] = ni ka1 gain
    //par[8] = ni ka1 mean
    //par[9] = ni ka1 sigma
  
    
    Double_t niKa1 = par[7]/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-par[8])*(xx-par[8]))/(2*par[9]*par[9]));
    
    Double_t niKa2Gain = par[7] * 0.51;
    Double_t niKa2Mean = par[8] - 17.3;
    
    Double_t niKa2 = niKa2Gain/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-niKa2Mean)*(xx-niKa2Mean))/(2*par[9]*par[9]));

    
    Double_t roiCuFunc = back + cuKa1 + cuKa2 + cuKb + niKa1 + niKa2;
    
   
    return roiCuFunc;    
}

Double_t funcSigBg(Double_t *x, Double_t *par){
    
    // fit of the region of roi, nickel, and cu ka kb
    // this is for fitting an already scaled histogram
    
    
    Double_t xx = x[0];
    
    //par[0] = background constant
    //par[10] = background slope
    
    Double_t back = par[0] + (xx - 7000) * par[10];
    
    //par[1] = cu ka1 gain
    //par[2] = cu ka1 mean
    //par[3] = cu ka1 sigma
    
    Double_t cuKa1 = par[1]/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-par[2])*(xx-par[2]))/(2*par[3]*par[3]));
    
    Double_t cuKa2Gain = par[1] * 0.51;
    Double_t cuKa2Mean = par[2] - 19.95;
    
    Double_t cuKa2 = cuKa2Gain/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-cuKa2Mean)*(xx-cuKa2Mean))/(2*par[3]*par[3]));
    
    //par[4] = cu kb gain
    //par[5] = cu kb mean
    //par[6] = cu kb sigma
    
    Double_t cuKb = par[4]/(sqrt(2*TMath::Pi())*par[6])*TMath::Exp(-((xx-par[5])*(xx-par[5]))/(2*par[6]*par[6]));
    
    //par[7] = ni ka1 gain
    //par[8] = ni ka1 mean
    //par[9] = ni ka1 sigma
  
    
    Double_t niKa1 = par[7]/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-par[8])*(xx-par[8]))/(2*par[9]*par[9]));
    
    Double_t niKa2Gain = par[7] * 0.51;
    Double_t niKa2Mean = par[8] - 17.3;
    
    Double_t niKa2 = niKa2Gain/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-niKa2Mean)*(xx-niKa2Mean))/(2*par[9]*par[9]));

    //par[11] = forbidden gauss gain
    // mean of forbidden transition fixed at 7729 eV
    
    Double_t forbGauss = par[11]/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-7729)*(xx-7729))/(2*par[3]*par[3]));
 
    Double_t roiCuFunc = back + cuKa1 + cuKa2 + cuKb + niKa1 + niKa2 + forbGauss;
    
   
    return roiCuFunc;    
}

struct GlobalChi2 { 
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,  ROOT::Math::IMultiGenFunction & f2) : 
      fChi2_1(&f1), fChi2_2(&f2) {}

   // parameter vector is first background (in common 1 and 2) and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[11]; // p1 is for background = without current histogram
      p1[0] = par[0]; // bg constant ..common parameter
      p1[1] = par[1]; // cu ka1 gain ... free 
      p1[2] = par[2]; // cu ka1 mean ... fixed 
      p1[3] = par[3]; // cu ka1 sigma ... free 
      p1[4] = par[4]; // cu kb gain ... free
      p1[5] = par[5]; // cu kb mean ... fixed
      p1[6] = par[6]; // cu kb sigma ... free
      p1[7] = par[7]; // ni ka1 gain ... free
      p1[8] = par[8]; // ni ka1 mean ... fixed
      p1[9] = par[9]; // ni ka1 sigma ... free
      p1[10] = par[10]; // background slope ... common

      double p2[12]; // parameters for the fit with signal with current 
      p2[0] = par[0]; // bg constant ..common parameter
      p2[1] = par[11]; // cu ka1 gain ... free 
      p2[2] = par[12]; // cu ka1 mean ... fixed 
      p2[3] = par[13]; // cu ka1 sigma ... free 
      p2[4] = par[14]; // cu kb gain ... free
      p2[5] = par[15]; // cu kb mean ... fixed
      p2[6] = par[16]; // cu kb sigma ... free
      p2[7] = par[17]; // ni ka1 gain ... free
      p2[8] = par[18]; // ni ka1 mean ... fixed
      p2[9] = par[19]; // ni ka1 sigma ... free
      p2[10] = par[10]; // background slope ... common
      p2[11] = par[20]; // forbidden gauss gain ... free
      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   } 

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

void combinedFit(Int_t reBin) { 

  TFile *fIN = new TFile("energyHistograms_40keV.root");
  Int_t nPar = 21; // does not change in the whole code automatically -> strg + f 21
  Int_t lowerL = 7000;
  Int_t upperL = 9500;
  
  TH1F *hSB = (TH1F*)fIN->Get("withCurrentSum");
  TH1F *hB   = (TH1F*)fIN->Get("noCurrentSmallSum"); 
  
  // setting the initial parameters for the fit
  Double_t parInit[21] = { 510. * reBin/25 , 310000. * reBin/25 , 8047.78 , 75. , 79000. * reBin/25 , 8905.29 , 80. , 14900. * reBin/25 , 7478.15 , 72. , -0.02 ,// parameters 0 ... 10 for background and common param
                                             312000. * reBin/25 , 8047.78 , 80. , 78500. * reBin/25 , 8905.29 , 78. , 12000. * reBin/25 , 7478.15 , 71. , 0. }; // parameters 11 ... 20 including signal

  hSB->Rebin(reBin);
  hB->Rebin(reBin);
  

  TF1 *fB = new TF1("fB", funcBg, lowerL, upperL, 11 );
  TF1 *fSB = new TF1("fSB", funcSigBg, lowerL, upperL, 12 );
  
//  for(Int_t i = 0; i < 10; i++){
//      
//      fB->SetParameter(i,parInit[i]);
//      
//  }
//  
  //fB->Draw();

  // perform now global fit 


  ROOT::Math::WrappedMultiTF1 wfB(*fB,1); // no idea 
  ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);

  ROOT::Fit::DataOptions opt; // simple structure holding the options on how the data are filled
  ROOT::Fit::DataRange rangeB; // class desribing ranges of data; SetRange(unsigned int icoord, double xmin, double xmax) 
  // set the data range
  rangeB.SetRange(lowerL,upperL); // sets range for the first coordinate -> it is apparently the range of the fit function
  ROOT::Fit::BinData dataB(opt,rangeB); 
  ROOT::Fit::FillData(dataB, hB); // fill the data vector from a TH1

  ROOT::Fit::DataRange rangeSB; 
  rangeSB.SetRange(lowerL,upperL);
  ROOT::Fit::BinData dataSB(opt,rangeSB); 
  ROOT::Fit::FillData(dataSB, hSB);

  ROOT::Fit::Chi2Function chi2_B(dataB, wfB); // no idea
  ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);

  GlobalChi2 globalChi2(chi2_B, chi2_SB);

  ROOT::Fit::Fitter fitter;

  

  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(21,parInit); // i guess this sets the initial parameters
  // fix some parameters 
  fitter.Config().ParSettings(2).Fix();
  fitter.Config().ParSettings(5).Fix();
  fitter.Config().ParSettings(8).Fix();
  fitter.Config().ParSettings(12).Fix();
  fitter.Config().ParSettings(15).Fix();
  fitter.Config().ParSettings(18).Fix();
  
//  fitter.Config().ParSettings(3).Fix();
//  fitter.Config().ParSettings(6).Fix();
//  fitter.Config().ParSettings(9).Fix();
//  fitter.Config().ParSettings(13).Fix();
//  fitter.Config().ParSettings(16).Fix();
//  fitter.Config().ParSettings(19).Fix();
  
  
  fitter.Config().ParSettings(0).SetName("Common Background Constant");
  fitter.Config().ParSettings(1).SetName("Cu Ka1 Gain BG");
  fitter.Config().ParSettings(2).SetName("Cu Ka1 Mean");
  fitter.Config().ParSettings(3).SetName("Cu Ka1 Sigma BG");
  fitter.Config().ParSettings(4).SetName("Cu Kb Gain BG");
  fitter.Config().ParSettings(5).SetName("Cu Kb Mean");
  fitter.Config().ParSettings(6).SetName("Cu Kb Sigma BG");
  fitter.Config().ParSettings(7).SetName("Ni Ka1 Gain BG");
  fitter.Config().ParSettings(8).SetName("Ni Ka1 Mean");
  fitter.Config().ParSettings(9).SetName("Ni Ka1 Sigma BG");
  fitter.Config().ParSettings(10).SetName("Common Background Slope");
  fitter.Config().ParSettings(11).SetName("Cu Ka1 Gain SIG");
  fitter.Config().ParSettings(12).SetName("Cu Ka1 Mean");
  fitter.Config().ParSettings(13).SetName("Cu Ka1 Sigma SIG");
  fitter.Config().ParSettings(14).SetName("Cu Kb Gain SIG");
  fitter.Config().ParSettings(15).SetName("Cu Kb Mean");
  fitter.Config().ParSettings(16).SetName("Cu Kb Sigma SIG");
  fitter.Config().ParSettings(17).SetName("Ni Ka1 Gain SIG");
  fitter.Config().ParSettings(18).SetName("Ni Ka1 Mean");
  fitter.Config().ParSettings(19).SetName("Ni Ka1 Sigma SIG");
  fitter.Config().ParSettings(20).SetName("Forbidden Gauss Gain");

  
  
  
  // first fit results:
  

  
  // set limits 
  fitter.Config().ParSettings(0).SetLimits(400. * reBin/25,600. * reBin/25);
  fitter.Config().ParSettings(1).SetLimits(200000. * reBin/25,400000. * reBin/25);

  fitter.Config().ParSettings(3).SetLimits(72.,80.);
  fitter.Config().ParSettings(4).SetLimits(70000. * reBin/25,90000. * reBin/25);

  fitter.Config().ParSettings(6).SetLimits(75.,85.);
  fitter.Config().ParSettings(7).SetLimits(10000. * reBin/25,20000. * reBin/25);

  fitter.Config().ParSettings(9).SetLimits(60. ,85. );
  fitter.Config().ParSettings(10).SetLimits(-0.05,0.05);
  fitter.Config().ParSettings(11).SetLimits(200000. * reBin/25,400000. * reBin/25);
  
  fitter.Config().ParSettings(13).SetLimits(75.,85.);
  fitter.Config().ParSettings(14).SetLimits(70000. * reBin/25,90000. * reBin/25);

  fitter.Config().ParSettings(16).SetLimits(75.,85.);
  fitter.Config().ParSettings(17).SetLimits(10000. * reBin/25,15000. * reBin/25);

  fitter.Config().ParSettings(19).SetLimits(65.,75.);
  //fitter.Config().ParSettings(20).SetLimits(0.,5000.);

  fitter.FitFCN(nPar,globalChi2,parInit,dataB.Size()+dataSB.Size()); // bool ROOT::Fit::Fitter::FitFCN (unsigned int  	npar,Function &  fcn,const double * params = 0,unsigned indataSize = 0,bool chi2fit = false )
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  TCanvas * c1 = new TCanvas("Simultaneous fit");
  c1->Divide(1,2);
  c1->cd(1);
  gStyle->SetOptFit(1111);

  fB->SetParameter(0, result.Parameter(0) ); 
  fB->SetParameter(1, result.Parameter(1) ); 
  fB->SetParameter(2, result.Parameter(2) ); 
  fB->SetParameter(3, result.Parameter(3) ); 
  fB->SetParameter(4, result.Parameter(4) ); 
  fB->SetParameter(5, result.Parameter(5) ); 
  fB->SetParameter(6, result.Parameter(6) ); 
  fB->SetParameter(7, result.Parameter(7) ); 
  fB->SetParameter(8, result.Parameter(8) ); 
  fB->SetParameter(9, result.Parameter(9) ); 
  fB->SetParameter(10,result.Parameter(10));
  
  fB->SetParError(0, result.Error(0) );
  fB->SetParError(1, result.Error(1) );
  fB->SetParError(2, result.Error(2) );
  fB->SetParError(3, result.Error(3) );
  fB->SetParError(4, result.Error(4) );
  fB->SetParError(5, result.Error(5) );
  fB->SetParError(6, result.Error(6) );
  fB->SetParError(7, result.Error(7) );
  fB->SetParError(8, result.Error(8) );
  fB->SetParError(9, result.Error(9) );
  fB->SetParError(10, result.Error(10) );
  
  fB->SetChisquare(result.MinFcnValue() );
  fB->SetNDF(result.Ndf() );
  hB->GetXaxis()->SetRangeUser(lowerL,upperL);
  hB->GetListOfFunctions()->Add(fB);
  //gPad->SetLogy();
  //fB->SetRange(rangeB().first, rangeB().second);   
  hB->Draw(); 

  c1->cd(2);
  
  fSB->SetParameter(0, result.Parameter(0)  ); 
  fSB->SetParameter(1, result.Parameter(11) ); 
  fSB->SetParameter(2, result.Parameter(12) ); 
  fSB->SetParameter(3, result.Parameter(13) ); 
  fSB->SetParameter(4, result.Parameter(14) ); 
  fSB->SetParameter(5, result.Parameter(15) ); 
  fSB->SetParameter(6, result.Parameter(16) ); 
  fSB->SetParameter(7, result.Parameter(17) ); 
  fSB->SetParameter(8, result.Parameter(18) ); 
  fSB->SetParameter(9, result.Parameter(19) ); 
  fSB->SetParameter(10,result.Parameter(10) );
  fSB->SetParameter(11,result.Parameter(20) );
  
  fSB->SetParError(0, result.Error(0)  );
  fSB->SetParError(1, result.Error(11) );
  fSB->SetParError(2, result.Error(12) );
  fSB->SetParError(3, result.Error(13) );
  fSB->SetParError(4, result.Error(14) );
  fSB->SetParError(5, result.Error(15) );
  fSB->SetParError(6, result.Error(16) );
  fSB->SetParError(7, result.Error(17) );
  fSB->SetParError(8, result.Error(18) );
  fSB->SetParError(9, result.Error(19) );
  fSB->SetParError(10,result.Error(10) );
  fSB->SetParError(11,result.Error(20) );
  
  fSB->SetChisquare(result.MinFcnValue() );
  fSB->SetNDF(result.Ndf() );
  fSB->SetRange(rangeSB().first, rangeSB().second); 
  hSB->GetXaxis()->SetRangeUser(lowerL,upperL);
  hSB->GetListOfFunctions()->Add(fSB);
  //gPad->SetLogy();
  hSB->Draw(); 
  
  Double_t sigCounts = fSB->GetParameter(11);
  Double_t sigError = fSB->GetParError(11);
  cout << endl << "The number of counts from PEP violating transitions is: " << sigCounts / reBin << " +- " << sigError / reBin << endl;
  cout << "The 3 sigma upper limit on PEP violating tranitions is therefore: " << 3 * sigError / reBin << " counts" << endl;

}