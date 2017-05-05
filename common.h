/*******************************

  common.h 
  04 Mar. 2014
  H. Shi

********************************/

#ifndef _COMMON_H_
#define _COMMON_H_

////////////////////
// Modules
const static Int_t   V1190_TDC_CH  = 16;
const static Int_t   V792_QDC_CH   = 32;
//const static Int_t   V785N_ADC_CH  = 16;
const static Int_t   V785_ADC_CH   = 16;

const static Int_t   TDC_REF_CH    = 14;

const static Int_t   NQDC_CH       = 32;

//const static Float_t NSECTDCCH       = 0.8;  // with pair mode before 2014 02 25
const static Float_t NSECTDCCH       = 0.1;  // with pair mode before 2014 02 25


////////////////////
// Detectors 
//const static Int_t   NSiPM    = 64;   // For final setup
const static Int_t   NSiPM    = 32;   // For test setup of Phase one
const static Int_t   NScinti  = 32;   // For final setup
const static Int_t   NPADC_CH = 32;   
const static Int_t   ENABLED_NPADC_CH = 16;   
const static Int_t   NCHIP    = 1;    // Number of SDD chip
const static Int_t   NARRAY   = 2;    // Two arrays each with 3 SDDs
const static Int_t   NSDD     = 6;    // For Phase one setup

const static Double_t SiW     = 3.81;  //energy for production of electron hole pair of Silicon @ 77 Kelvin
//const static Double_t SiW     = 1.16;  // Band gap width for Silicon
const static Double_t SiEsc   = 1.74;  // Silicon Ka energy for escape peak

///////////////////
// For data decode
const static Int_t   MAXCHAR = 512;
const static Int_t   SEC_DAY = 86400;
const static Int_t   FIRSTJAN2014 = 1388530800;   // unix time for 2014 Jan 01 

const static Double_t SQRT2       = 1.414213562373095;
const static Double_t SIG2FWHM    = 2.35; 


/////////////////////
// File directories:
const static TString  WORK_PATH  = "/home/andreas/vip2";
const static TString  GRAPH_PATH = WORK_PATH + "/graphs";
const static TString  DATA_PATH  = WORK_PATH + "/data";
const static TString  LIST_PATH  = WORK_PATH + "/filelist";
const static TString  TXT_PATH   = DATA_PATH + "/txt";
const static TString  BIN_PATH_SMI   = DATA_PATH + "/bin/SMI";
const static TString  ROOT_PATH  = DATA_PATH + "/root";
const static TString  SLOW_TEXT_PATH  = DATA_PATH + "/slow/text";
const static TString  SLOW_ROOT_PATH  = DATA_PATH + "/slow/root";
const static TString  BIN_PATH_LNGS   = DATA_PATH + "/bin/LNGS";
const static TString  ROOT_PATH_LNGS  = DATA_PATH + "/root/LNGS";
const static TString  ROOT_PATH_SMI  = DATA_PATH + "/root/SMI";



//////////////////////////////////////
// Constants from common.h for sidt
const static Double_t  SidtSddSigAt6keV = 63.;      // typical sidt SDD resolution - sigma
const static Double_t  DearCcdSigAt6keV = 63.;      // typical DEAR CCD resolution - sigma
const static Double_t  KpXSiLiSigAt6keV = 170.;     // typical KpX  SDD resolution - sigma

const static Double_t  SDDTDCNSCH = 8.3;        // Channel to ns correspondence of SDD TDC
const static Double_t  KDTDCNSCH  = 0.035;      // Channel to ns correspondence of KD  TDC



////////////////////
// Macro
#define sq(x) ((x)*(x))


#endif
