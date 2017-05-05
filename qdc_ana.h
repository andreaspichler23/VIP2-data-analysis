////////////////////////////////
//      qdc_ana.h    
//       04 Mar. 2014 
//       H. Shi
///////////////////////////////
#ifndef _QDC_ANA_H_
#define _QDC_ANA_H_

#include  "common.h"

/////////////////
// Functions 

Int_t     FitSpectrum( TH1F *h );        // Fit spectrum and return histo threshold bin for pedestal cut 
Int_t     GetThreshold( TH1F *h,  Double_t pedestal_peak,  Double_t signal_peak );     // Return the minimum value bin. 

Double_t  PedestalFunction( Double_t *x, Double_t *par );       // Gaussian for pedestal.
Double_t  CosmicHitFunction( Double_t *x, Double_t *par );      // Vavilov for cosmic ray hits. 
Double_t  QdcFunction( Double_t *x, Double_t *par );            // Total fit function

#endif
