/*********************************************
            SDDcalib.h 
              Aug. 2014
              H. Shi
**********************************************/
#ifndef SDDcalib_h_included_
#define SDDcalib_h_included_

// For Ti Mn calibration with 55Fe source and Ti foil, "source" can be void by default.
// 2014 09 08 
void  SingleSDDCalib( ID nSDD,  ID nPADC,  TString rootf, TString source, TString place );
void  PhaseOneCalib( TString rootf,  TString source, TString place );

#endif
