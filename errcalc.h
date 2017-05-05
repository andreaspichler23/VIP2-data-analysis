//////////////////////////////////
//     errcalc.h
//     Apr 1 2011 
//     error propagation for add,
//       sub, mul, and div, with 
//       asymmetric error. 
//     Apr 12 2011 
//     added function for nolinear 
//       prop with one VE parameter 
//////////////////////////////////
#ifndef _ERRCALC_H_INCLUDED_
#define _ERRCALC_H_INCLUDED_

#include "ve.h"

Double_t GetError( Double_t, Double_t, Double_t,
		   Double_t, Double_t, Double_t, 
                   TString,  TString );

VE ErrProp( VE, VE, TString );
VE ErrProp( VE, Double_t, Double_t, TString);   // April 12 2011 

#endif


