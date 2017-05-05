/*******************************

  ReadFullDaqBinData.h 
  2013/12/16
  H. Shi
  04 Mar. 2014 
     Separate the modules and 
     detectors part to common.h

********************************/

#ifndef READFULLDAQDATA_H_INCLUDED_
#define READFULLDAQDATA_H_INCLUDED_

#include  "common.h"

////////////////////
// For data ana

typedef struct{ 
    Int_t  header; 

    Int_t  rate; 
    Short_t  year; 
    Short_t  month; 
    Short_t  day; 
    Short_t  hour; 
    Short_t  min;
    UInt_t    ut;    // unix time 

    UInt_t   evid;   // incremented every loop iteration
    Int_t  evid1;  // From qdc module
    Int_t  evid2;  // From padc module
    Int_t  evid3;  // From tdc module
    Int_t  mul;    // "Multiplicity" for number of SiPM hits per event

    Int_t  qdc[V792_QDC_CH];
    Int_t  tdc[V1190_TDC_CH];
    Int_t  tot[V1190_TDC_CH];   // Time over threshold for SiPM signal
    Int_t  padc[V785_ADC_CH];   // SDD peak ADC 

    Int_t  layer;  // position of scinti. bot to top   0, 1, 2, 3, 4, 5
    Int_t  col;    // position of scinti. side to side 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    Int_t  trgid;     // type of trigger
    Int_t  hr; 

} EventStruct;


/////////////////
// Functions 
Int_t ReadData( TString filelist, TString rootfile); 
void  MakeTree( TString datafile, TString rootfile);
void  MakePlots( TString rootfile);


#endif
