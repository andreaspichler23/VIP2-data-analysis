/********************************************

      PhaseOneSDDConnectionTable.h 
           16 May 2014
           H. Shi

*********************************************/ 
#ifndef _PHASEONESDDCONNECTION_TABLE_H_
#define _PHASEONESDDCONNECTION_TABLE_H_

#include  "common.h"

// six SDDs from SIDDHARTA chip1000 for Phase One setup 
const static Int_t PadcToSDD[8] = {
     1,  -1,  2,  3,  -1,  4,  5,  6  
};

const static Int_t SDDToPadc[7] = {
     -1, 0,  2,  3,  5,  6,  7  
};

const static Int_t TdcChSDD[NSDD] = {
     53,  52, 51, 50, 40, 49        // 40 for SDD 5 is dummy, actually not connected to TDC
};


const static Int_t PcbSDD[NSDD] = {    // Board number 
    2, 2, 2, 1, 1, 1
};

const static Int_t BcSDD[NSDD] = {     // BC voltage for SDD
    1, 2, 3, 1, 2, 3

};

const static Int_t GroupSDD[NSDD] = {   // 1 for A side, 2 for B side in SIDDHARTA configuration
                                        // for example, Vb1A voltage in SIDDHARTA slow log is for SDD 1 on board 2 in VIP2
    1, 1, 1, 2, 2, 2
};

const static Int_t AdcCuts[6] = { 320, 300, 320, 290, 280, 270};

#endif

