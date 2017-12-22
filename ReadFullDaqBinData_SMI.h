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
    Short_t  sec;
    Short_t  time[6];
    UInt_t    ut;    // unix time 
    Int_t    clk;

    UInt_t   evid;   // incremented every loop iteration
    Short_t  evid1;  // From qdc module
    Short_t  evid2;  // From padc module
    Short_t  evid3;  // From tdc module
    UInt_t   evnr;   // event number in rootfile
    //Short_t  mul;    // "Multiplicity" for number of scintillator hits per event

    Short_t  qdc[V792_QDC_CH];
    Short_t  tdc[V1190_TDC_CH];
    Short_t  tot[V1190_TDC_CH];   // Time over threshold for SiPM signal
    Short_t  padc[V785_ADC_CH];  // SDD peak ADC , branch is named "adc"
    Short_t  tdc_dig[9]; // 4 x tdc channel data [to, ti, bi, bo]; 0...tdc < 0; 1: 4000  tdc < 7000; 2: tdc> 7000; 1 x outer layer 1/0; 1 x inner layer 1/0; 1 x OR of inner + outer; 1 x AND of inner + outer; 1 x tdc layer mul
    Short_t  tdc_dig_smi[13];
    Short_t  qdc_dig[37]; // 32 x qdc channel data [0...31]; 0: qdc < thr; 1: qdc > thr; 1 x outer layer 1/0; 1 x inner layer 1/0; 1 x OR of inner + outer; 1 x AND of inner + outer; 1 x scinti mul
    Short_t  adc_dig[8]; //  6 x adc > thr; 1 x OR of all of them; 1 x sdd mul
//Short_t  padc[32]; 

    Short_t  layer;  // position of scinti. bot to top   0, 1, 2, 3, 4, 5
    Short_t  col;    // position of scinti. side to side 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    Short_t  trgid;     // type of trigger: -1: tdc[12] OR tdc[13]; 1: SDD only; 2: scintillator only; 3 coincidence; 4: SDD + 1 tdc layer
    //Short_t  hr;
    
    Double_t slow[19];
    //TObjString binFile;

 //   Short_t qdc_flags[3];  // qdc_outer flag, qdc_ inner flag and qdc-coincidence flag (inner AND outer)
 //   Short_t 


} EventStruct;


/////////////////
// Functions 
Int_t ReadData( TString datafilename, TString root_file, TString place, int& root_file_event_counter ); 
void MakeTree( TString filelist, TString rootfilename, TString place );
//void MakePlots(TString rootfile, Int_t qdc_th[NScinti] );
Int_t GetSlowParameter(UInt_t time, Double_t *slow);


#endif
