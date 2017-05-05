////////////////////////////
//   ve.C 
//   Mar 25 2011 
////////////////////////////
#ifndef _VE_C_INCLUDED
#define _VE_C_INCLUDED

#include "ve.h"

VE SetVE( Double_t v, Double_t e )         
{
  VE SetVal;
  SetVal.v  = v;
  SetVal.ep = e;
  SetVal.em = e;
  return SetVal;
}

VE SetVE( Double_t v, Double_t ep, Double_t em )         
{
  VE SetVal;
  SetVal.v  = v;
  SetVal.ep = ep;
  SetVal.em = em;
  return SetVal;
}

VE SetVE(VE ve)         
{
  VE SetVal;
  SetVal.v = ve.v;
  SetVal.ep = ve.ep;
  SetVal.em = ve.em;
  return SetVal;
}

void VE2Double(VE ve[], Double_t v[], 
               Double_t ep[], Double_t em[], Int_t N)
{
  for (Int_t i=0; i<N; i++){
    v[i]  = ve[i].v;
    ep[i] = ve[i].ep;
    em[i] = ve[i].em;
  }
}

#endif
