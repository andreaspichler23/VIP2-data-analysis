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
    Short_t  header; 

    Short_t  rate; 
    Short_t  year; 
    Short_t  month; 
    Short_t  day; 
    Short_t  hour; 
    Short_t  min;
    UInt_t    ut;    // unix time 

    UInt_t   evid;   // incremented every loop iteration
    Short_t  evid1;  // From qdc module
    Short_t  evid2;  // From padc module
    Short_t  evid3;  // From tdc module
    Short_t  mul;    // "Multiplicity" for number of SiPM hits per event

    Short_t  qdc[V792_QDC_CH];
    Short_t  tdc[V1190_TDC_CH];
    Short_t  tot[V1190_TDC_CH];   // Time over threshold for SiPM signal
    Short_t  padc[V785_ADC_CH];   // SDD peak ADC 

    Short_t  layer;  // position of scinti. bot to top   0, 1, 2, 3, 4, 5
    Short_t  col;    // position of scinti. side to side 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    Short_t  trgid;     // type of trigger
    Short_t  hr;


} EventStruct;


/////////////////
// Functions 
Int_t ReadData( TString filelist, TString rootfile); 
void  MakeTree( TString datafile, TString rootfile);
void  MakePlots( TString rootfile);


#endif
