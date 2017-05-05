//////////////////////////////////
//     errcalc.C
//     Mar 25 2011 
//     error propagation for add,
//       sub, mul, and div, with 
//       asymmetric error. 
//     Apr 12 2011 
//     added nolinear propagation
//       for sqrt and log, expo
//////////////////////////////////
#ifndef _ERRCALC_INCLUDED
#define _ERRCALC_INCLUDED

#include "TMath.h"
#include "TString.h"
#include "Riostream.h"
#include "ve.h"
//#include "errcalc.h"

Double_t GetError( Double_t x1, Double_t x1errp, Double_t x1errm,
		   Double_t x2, Double_t x2errp, Double_t x2errm, 
                   TString ope="ADD",  TString erType="plus" )
{
    Double_t r1p = x1errp / x1;
    Double_t r2p = x2errp / x2;
    Double_t r1m = x1errm / x1;
    Double_t r2m = x2errm / x2;
    Double_t errp = 0;
    Double_t errm = 0;
    Double_t err  = 0;
    if ( (ope == "ADD") || (ope == "SUB")) {
        errp = sqrt( x1errp*x1errp + x2errp*x2errp );
        errm = sqrt( x1errm*x1errm + x2errm*x2errm );
    } else if ( ope == "MUL" ) {
        errp = (x1*x2)*sqrt( r1p*r1p + r2p*r2p );
        errm = (x1*x2)*sqrt( r1m*r1m + r2m*r2m );
    } else if ( ope == "DIV" ) {
        errp = (x1/x2)*sqrt( r1p*r1p + r2p*r2p );
        errm = (x1/x2)*sqrt( r1m*r1m + r2m*r2m );
    } else {
        cout << "--- Error : 'ope' should be 'ADD','SUB','MUL','DIV' ---" << endl;
    }

    if ( erType == "plus" ){
        err = errp;
    } else if ( erType == "minus" ){
        err = errm;
    } else {
        cout << "--- Error : 'erType' should be 'plus' or 'minus' ---" << endl;
    }
    return err;
}


// ErrProp function for add, sub, mul, div
VE ErrProp( VE v1, VE v2, TString ope="ADD")
{
  VE val;
  if      (ope == "ADD") val.v = v1.v + v2.v;
  else if (ope == "SUB") val.v = v1.v - v2.v;
  else if (ope == "MUL") val.v = v1.v * v2.v;
  else if (ope == "DIV") val.v = v1.v / v2.v;
  else {
    cout << "--- Error : 'ope' should be 'ADD','SUB','MUL','DIV' ---" << endl;
  }
  val.ep = GetError( v1.v, v1.ep, v1.em, v2.v, v2.ep, v2.em, ope, "plus");
  val.em = GetError( v1.v, v1.ep, v1.em, v2.v, v2.ep, v2.em, ope, "minus");
  return val;
}

// ErrPro for nolinear function propagation
//    Function                  Variance
// f = a A^(+-b)          sig_f / f = b * sig_A / A
// f = a e^(+-bA)         sig_f / f = b * sig_A
// f = a log( +- bA)      sig_f     = a * sig_A / A
// f = a^(+-b A)          sig_f / f = b * log(a) * sig_A
//     coef1 : parameter a, default as 1.;
//     coef2 : parameter b, default as 1.;
VE ErrProp( VE v0, Double_t coef1=1., Double_t coef2=1., TString ope="POW" )
{
    VE val; 
    if    ( ope == "POW" ) {
        val.v  = coef1 * TMath::Power( v0.v,  coef2 );
        val.ep = val.v * sqrt( coef2 * coef2 ) * v0.ep / v0.v;
        val.em = val.v * sqrt( coef2 * coef2 ) * v0.em / v0.v;
    } else if ( ope == "EXP" ) {
        val.v  = coef1 * TMath::Exp( coef2 * v0.v );
        val.ep = val.v * coef2 * v0.ep;
        val.em = val.v * coef2 * v0.em;
    } else if ( ope == "LOG" ) {
        val.v  = coef1 * TMath::Log( coef2 * v0.v );
        val.ep = coef1 * v0.ep / v0.v;
        val.em = coef1 * v0.em / v0.v;
    }

    return val;
}

#endif


