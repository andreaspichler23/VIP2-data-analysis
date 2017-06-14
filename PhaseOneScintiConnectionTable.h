/********************************************

      PhaseOneScintiConnectionTable.h 
           01 Mar 2014
           H. Shi

          altered:
        all the time
      Andreas Pichler

*********************************************/ 
#ifndef _PHASE_ONE_CONNECTION_TABLE_H_
#define _PHASE_ONE_CONNECTION_TABLE_H_

#include  "ReadFullDaqBinData.h"

// scintillator layer to tdc connections:

const static Int_t Tdc2Column[8] = { 0, 2, 3, 1, 2, 2, 4, 2};

const static Int_t Tdc2Layer[8] = { 2, 1, 2, 2, 4, 3, 2, 0};

const static Int_t Tdc2Column_lngs[11] = { 0, 2, 3, 1, 2, 2, -1, -1, 2, -1, 4};

const static Int_t Tdc2Layer_lngs[11] = { 2, 1, 2, 2, 4, 3, -1, -1, 0, -1, 2};

const static Int_t TdcChannels[8] = {0, 1, 2, 3, 4, 5, 8, 10};

const static Int_t TdcConnections_smi[8] = {4, 5, 1, 7, 0, 3, 2, 6}; // until end of january 2016 at lngs; in the order: top outer, top inner, bottom inner, bottom outer, left outer, left inner, right inenr, right outer

const static Int_t TdcConnections_lngs[8] = {4, 5, 1, 8, 0, 3, 2, 10}; // starting end of january 2016 at lngs; in the order: top outer, top inner, bottom inner, bottom outer, left outer, left inner, right inenr, right outer

//////////////

//                                    bi 1           bi 6    bo 1           bo 6
const static Int_t Qdc2Outer[NSiPM] = { 0, 0, 0, 0, 0, 0,     1, 1, 1, 1, 1, 1,   
//ti 1         ti 6   to 1          to 6
  0, 0, 0, 0, 0, 0,    1, 1, 1, 1, 1, 1,     0, 0,    1, 1,    0, 0,   1, 1      
};

//////////////


const static Int_t Tdc2Outer_smi[8] = { 1, 0, 0, 0, 1, 0, 1, 1 };


/////////////

const static Int_t QdcCuts[NSiPM] = {682,889,1000,709,841,779,516,517,607,551,492,432,650,386,
900,743,720,378,626,567,608,746,722,700,800,760,563,600,798,796,880,800};

const static Int_t QdcCuts_80nsGate[NSiPM] = {270,325,345,370,315,310,330,390,385,370,410,340,360,320,
420,335,400,325,400,345,360,390,430,350,290,370,380,315,360,380,270,370};
//12:430->370                               QDC 0   1   2   3   4   5   6   7   8   9   10  11  12 13  14  15 (7: 435->375, 8: 435->360, 10: 465->400, 11: 410_>350)
const static Int_t QdcCuts_100nsGate[NSiPM] = {295,370,385,430,350,345,420,375,360,440,400,350,430,400,470,385,
450,360,400,400,405,425,490,385,340,420,425,350,400,420,330,420};
//16 17 18  19  20   21  22 23  24   25  26 27  28  29   30  31
// 16:480->450, 18:450->400, 19:380->400
// NSiPM 20 on 01 Mar 2014 for Phase One setup
const static Int_t QdcChSiPM[NSiPM] = {
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
     16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
};


// general serial ID for the veto scintillators
// For final configuration 
// 29 04 2015
const static Int_t Qdc2Scinti[NScinti] = {
       // bottom
       6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5, 
       // top
       25,24, 23, 22, 21, 20, 31, 30, 29, 28, 27, 26,
       // left
       13, 15, 12, 14,
       // right
       16, 18, 17, 19
//    17, 17, 18, 18, 19, 19, 20, 20, 21, 21,    22, 22, 24, 24, 25, 25, 26, 26, 27, 27
};


////////
// SiPM QDC channel corresponds to TDC 
const static Int_t Qdc2Tdc_smi[32] = {
   // 0               5    6           11   12          17
    1, 1, 1, 1, 1, 1,   7, 7, 7, 7, 7, 7,   5, 5, 5, 5, 5, 5, 
   // 18           23   24                           31
    4, 4, 4, 4, 4, 4,   3, 3,    0, 0,    2, 2,    6, 6
};

const static Int_t Qdc2Tdc_lngs[32] = {
   // 0               5    6           11   12          17
    1, 1, 1, 1, 1, 1,   8, 8, 8, 8, 8, 8,   5, 5, 5, 5, 5, 5, 
   // 18           23   24                           31
    4, 4, 4, 4, 4, 4,   3, 3,    0, 0,    2, 2,    10, 10
};

const static Int_t Tdc2Scinti[32] = {
   // 0               4    5                   10   11              15
    27, 26, 25, 24, 22,   -1, -1, -1, -1, -1, -1,  17, 18, 19, 20, 21, 
   // 16             20    21                   26   27             31 
    27, 26, 25, 24, 22,   -1, -1, -1, -1, -1, -1,  17, 18, 19, 20, 21
};




/////// Bottom layer is layer 0, then 1, 2, 3, 4, 5,  including the final all scintillators
// For final configuration 
// Scintillator to position : 
// 29 04 2015
const static Int_t ScintiLayer[NScinti] = {
    // 0,                       11
    0, 0, 0, 0, 0, 0,   1, 1, 1, 1, 1, 1,
    // 12,                      19 
    2, 2, 3, 3,    2, 2, 3, 3,
    // 20,                      31
    4, 4, 4, 4, 4, 4,   5, 5, 5, 5, 5, 5
};

/////// Number of columns includes the final all scintillators
// parameter is the scinti ID, table for PhaseOne only 
const static Int_t ScintiColumn[NScinti] = {
    // 0                                         11
     1, 2, 3, 4, 5, 6,    1, 2, 3, 4, 5, 6,
    // 12,                                       19 
     1, 2, 1, 2,          5, 6, 5, 6,
    // 20,                      31
     1, 2, 3, 4, 5, 6,    1, 2, 3, 4, 5, 6
};





// parameter is the qdc channel, table for PhaseOne only 
// For final configuration 
// 29 04 2015
// QDC channel to position :
const static Int_t Qdc2Layer[NQDC_CH] = {
   // 0                11 
     1, 1, 1, 1, 1, 1,   0, 0, 0, 0, 0, 0,
   // 12               23
     4, 4, 4, 4, 4, 4,   5, 5, 5, 5, 5, 5,
   // 24               31 
     3, 2, 2, 3,         3, 2, 3, 2 
};

const static Int_t Qdc2Column[NQDC_CH] = {
   // 0                11
     1, 2, 3, 4, 5, 6,    1, 2, 3, 4, 5, 6,
   // 12               23 
     6, 5, 4, 3, 2, 1,    6, 5, 4, 3, 2, 1,
   // 24               31
     1, 1, 0, 0,          6, 6, 7, 7

};

#endif

