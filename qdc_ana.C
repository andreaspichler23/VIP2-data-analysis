///////////////////////////
//    qdc_ana.C
//    03 Mar. 2014 
//        H. Shi
///////////////////////////

#include  <stdio.h>
#include  <stdlib.h>
#include  <iostream>

#include  "TMath.h"
#include  "TF1.h"
#include  "TH1.h"

#include  "common.h"
#include  "qdc_ana.h"
#include  "PhaseOneScintiConnectionTable.h"

/*

////// Fit qdc spectrum and return the valley value. 
Int_t  FitSpectrum( TH1F *h )
{

    const static Int_t    nPar  = 8;

    Double_t  par[nPar]   = { 0. };
    Double_t  parEr[nPar] = { 0. };

    Double_t l_bin = 100.,  u_bin = 4000.; 

    Double_t ped_m= 10.,    qdc_mpv = 1000.; 
    Int_t    cut_value = 0;   // cut value threshold for pedestal


    TF1 *fqdc = new TF1("fqdc", QdcFunction, l_bin, u_bin, nPar );

    fqdc->SetParameters( );
    fqdc->SetParLimits(0, );
    fqdc->SetParLimits(1, );
    fqdc->SetParLimits(2, );
    fqdc->SetParLimits(3, );
    fqdc->SetParLimits(4, );
    fqdc->SetParLimits(5, );
    fqdc->SetParLimits(6, );
    fqdc->SetParLimits(7, );

    h->Fit("fqdc", "RIE");

    fqdc->GetParameters( par );

    for( Int_t j=0; j<nPar; j++ ){
       parEr[j] = fqdc->GetParError(j);
    }


    ped_m   = par[0]; 
    qdc_mpv = par[6];

    cut_value = GetThreshold( h, ped_m, qdc_mpv );

    return cut_value;
}


Int_t GetThreshold( TH1F *h, Double_t ped_peak,  Double_t signal_peak )
{

    Int_t    min_bin = 100;    // The cut value for pedestal

    Double_t min     = 1e10;
    Double_t binCont = 0;

    Int_t lbin = 0.,  ubin = 0.;

    lbin = h->FindBin( ped_peak );
    ubin = h->FindBin( signal_peak );

    for( Int_t ibin = lbin; ibin <= ubin; ibin ++ ){

        binCont = h->GetBinContent( ibin );
        
        if ( binCont < min ) {
            min_bin = ibin;
            min = binCont;
        }
    }

    return min_bin;
}


//////// Gaussian pedestal function ////////
Double_t PedestalFunction( Double_t *x, Double_t *par ){
    Double_t  arg = 0.;
    Double_t  sig = par[1];  
    if ( sig != 0 ){
        arg = ( x[0] - par[0] ) /sig;
    }
    Double_t pedFunc = par[2] / sig * TMath::Exp( -0.5 * arg * arg );

    return pedFunc; 
}

//////// Vavilov cosmic ray hit function //////
// Vavilov function parameters : 
//         par[0] : kappa;     par[1] : beta;     par[3] : most probable value; 
Double_t CosmicHitFunction( Double_t *x, Double_t *par ) 
{ 

    Double_t cosmicHitFunc = TMath::Vavilov( ( x[0] - par[3] ) * par[4], par[0], par[1] ) * par[2] ;

    return cosmicHitFunc; 
}


/////// Fit function for Qdc spectrum //////
//  Parameters : 
//         par[0] : ped mean;  par[1] : ped sig;  par[2] : ped amp; 
//         par[3] : kappa;     par[4] : beta;     par[5] : vav amp;   par[6] : most probable value; 
Double_t QdcFunction( Double_t *x, Double_t *par ){
    Double_t qdcFunc = 0.;

    Double_t  arg = 0.;
    Double_t  sig = par[1];  
    if ( sig != 0 ){
        arg = ( x[0] - par[0] ) /sig;
    }
    Double_t pedFunc = par[2] / sig * TMath::Exp( -0.5 * arg * arg );

    Double_t cosmicHitFunc = TMath::Vavilov( ( x[0] - par[6] ) * par[7], par[3], par[4] ) * par[5] ;

    qdcFunc = pedFunc + cosmicHitFunc; 

    return qdcFunc;
}

*/

void GetQdcCounts(TString rootfile){

   TString rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

   Int_t counts;
   Int_t xmax = 4000;

   
for (Int_t i = 0; i < NSiPM; i++){

   Int_t qdc_ch = i;
   Int_t xmin = QdcCuts_80nsGate[qdc_ch];
 
   TH1F* hsq;
   TFile *f = new TFile(rootfilename, "READ");
   hsq = (TH1F*)f->Get(Form("hq[%d]",qdc_ch));

   TAxis *axis = hsq->GetXaxis();

   Int_t bmin = hsq->GetXaxis()->FindBin(xmin);
   Int_t bmax = hsq->GetXaxis()->FindBin(xmax);
  

   counts = hsq->Integral(bmin,bmax);

   counts -= hsq->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);

   cout << "In channel " << qdc_ch << ": " << counts << " counts" << endl;

   hsq->Delete();

}


}

