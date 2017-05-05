/***************************************** 
                SDDclass.h

                by H. Shi

      Originally from sidt/final/
      for SIDDHARTA setup. 
                
      Modifications:
        Feb. 23 2010: Added Member funciton
                      and variable for average 
                      rate and x tube current 
                      information. 
        June 22 2010: Added Ti Mn fit methods.

        May 2014    : modified for VIP-2 use.

        July 2016   : modified again for VIP-2 by Andreas Pichler

*****************************************/ 
#ifndef _SDDCLASS_H_INCLUDED_
#define _SDDCLASS_H_INCLUDED_

#include  "TLine.h"
#include  "TF1.h"
#include  "TFile.h"
#include  "TString.h"
#include  "TCanvas.h"

#include  "common.h"

typedef Double_t AnaPar;           // Response function parameter
typedef Double_t SlowPar;          // Slow log parameter 
typedef Int_t    ID;               // Geometry 
typedef UInt_t   Time;             // Date time

static const int maxPar = 100;          // upper limit for parameter number 2011/04/19
//static const int nPar = 25;
static const int nPar = 27; // current number of parameters for fit function

class SDDclass  {
    public:
        SDDclass( ID id,   ID padc_ch);            // Constructor
        ~SDDclass();                                         // Destructor

        //bool SetReadList( TString input_list, Int_t nf );
        void SetRootfile( TString ith_run, TString place );

        // Constant member functions return parameters
        ID GetBus()  const ;
        ID GetID()   const ;
        ID GetChip() const ;
        ID GetGrp()  const ;
        ID GetPos()  const ;
        ID GetRow()  const ;
        ID GetCol()  const ;
        ID GetRnk()  const ; 

        Time   GetTlabel()   const ;
        Time   GetUnixT()    const ;
        Time   GetRunTime() ;

        Int_t  GetNdf()       const ;
        AnaPar GetChi2()      const ;
        AnaPar GetFano()      const ;
        AnaPar GetFanoEr()    const ;
        AnaPar GetCstn()      const ;
        AnaPar GetCstnEr()    const ;
        AnaPar GetTika1()     const ;
        AnaPar GetTika1Er()   const ;
        AnaPar GetCuka1()     const ;
        AnaPar GetCuka1Er()   const ;
        AnaPar GetCuKaG()     const ;
        AnaPar GetMnka1()     const ;
        AnaPar GetMnka1Er()   const ;
        AnaPar GetMnKaG()     const ;
        AnaPar GetMnD()       const ;
        AnaPar GetMnDEr()     const ;
        AnaPar GetTiKaG()     const ;
        AnaPar GetFeG()       const ;
        AnaPar GetSigTika()   const ;
        AnaPar GetSigTikaEr() const ;
        AnaPar GetSigMnka()   const ;
        AnaPar GetSigMnkaEr() const ;
        AnaPar GetFwhm()      const ;
        AnaPar GetFwhmEr()    const ; 
        AnaPar GetEv2ch()     const ;
        AnaPar GetEv2chEr()   const ;
        AnaPar GetE2ChOffset()    const ;
        AnaPar GetE2ChOffsetEr()  const ;
        AnaPar GetTiTailG()   const ; // Gain of tail component
        //AnaPar GetMnTailG()   const ;
        AnaPar GetCuTailG()   const ;
        AnaPar GetPileG()     const ;
        AnaPar GetPileGEr()   const ;
        AnaPar GetTTiBeta()   const;
        AnaPar GetTMnBeta()   const;
        AnaPar GetTShift()    const;
        AnaPar GetPSigR()     const;
        AnaPar GetPSigREr()   const;

        /// Information for tail component
        AnaPar GetNtiK()      const;
        AnaPar GetNtiKEr()    const;
        AnaPar GetNcuK()      const;
        AnaPar GetNcuKEr()    const;
        //AnaPar GetNmnK()      const;
        //AnaPar GetNmnKEr()    const;
        AnaPar GetNtiT()      const;
        AnaPar GetNtiTEr()    const;
        //AnaPar GetNmnT()      const;
        //AnaPar GetNmnTEr()    const;
        AnaPar GetNcuT()      const;
        AnaPar GetNcuTEr()    const;
        //AnaPar GetRtiT()      const;
        //AnaPar GetRmnT()      const;
        //AnaPar GetRmnTEr()    const;
        //AnaPar GetRcuT()      const;
        //AnaPar GetRtiTEr()    const;
        //AnaPar GetRcuTEr()    const;
        //AnaPar GetAverageRate();


        /*
        SlowPar  GetVrefPOS()    const ;
        SlowPar  GetVrefPULSE()  const ;
        SlowPar  GetDacOffset()  const ;
        SlowPar  GetVsa()        const ;
        SlowPar  GetVsss()       const ;
        SlowPar  GetRxA()        const ;
        SlowPar  GetRxB()        const ;
        SlowPar  GetBfA()        const ;
        SlowPar  GetBfB()        const ;
        SlowPar  GetVb1A()       const ;
        SlowPar  GetVb1B()       const ;
        SlowPar  GetVb2A()       const ;
        SlowPar  GetVb2B()       const ;
        SlowPar  GetVb3A()       const ;
        SlowPar  GetVb3B()       const ;
        SlowPar  GetR1A()        const ;
        SlowPar  GetR1B()        const ;
        SlowPar  GetP12VA()      const ;
        SlowPar  GetVdd()        const ;
        SlowPar  GetM12VA()      const ;
        SlowPar  GetVss()        const ;
        SlowPar  GetP12VB()      const ;
        SlowPar  GetHv()         const ;
        SlowPar  GetM12VB()      const ;
        SlowPar  GetM5V()        const ;
        SlowPar  GetTsenseA()    const ;
        SlowPar  GetTsenseB()    const ;
        SlowPar  GetPt100()      const ;
        SlowPar  GetTADC()       const ;
        SlowPar  GetXtubeAm()    const ;

        void      SetSlowPar();
        */

        // Histogram analysis functions
        bool      SearchPeak();                       // return true for 2 peaks,also fit bg
        Int_t     FitTiCu(Bool_t saveflag );          // Add the part for tail analysis
                                                      //  on 11/13/2009
        Int_t     FitTiMnCu(Bool_t saveflag );        // Add the part for Mn 
        Int_t     FitTiMn(Bool_t saveflag );          // For Ti Mn fit 
        Int_t     FitTiCuZr(Bool_t saveflag );        // For X ray tube data
        
        void      GetHistoRange();
        void      GetFitParameters( TString source );
        void      CalcFwhmMn( TString source );
        void      PlotResidue( TString source );
        void      OpenCanvas();
        void      CloseCanvas(Bool_t saveplot);
        void      CloseRootFile();
        void      FitE2ChLine( TString source );  // here I added the parameter nSDD
        //void      GetParameter();


    private:
        // Hisotgram analysis variables
        TFile *froot;
        TH1F  *hcal;
        TH1F  *hraw;
        TH1F  *hcut;
        TH1F  *hkh;
        TH1F  *hev;

        TF1   *fcnFit;
        TF1   *fcnResult;

        TF1   *fcnCaka1;
        TF1   *fcnTika1;
        TF1   *fcnMnka1;
        TF1   *fcnCuka1;
        TF1   *fcnZrka1;

        TF1   *fcnTCaK1;
        TF1   *fcnTTiK1;
        TF1   *fcnTMnK1;
        TF1   *fcnTCuK1;
        TF1   *fcnTZrK1;

        TF1   *fcnTCaK2;
        TF1   *fcnTTiK2;
        TF1   *fcnTMnK2;
        TF1   *fcnTCuK2;
        TF1   *fcnTZrK2;

        TF1   *fcnPCaK1;
        TF1   *fcnPTiK1;
        TF1   *fcnPMnK1;
        TF1   *fcnPCuK1;

        TF1   *fcnEscCuka1;
        TF1   *fcnEscMnka1;
        TF1   *fcnEscTika1;
        TF1   *fcnEscZrka1;

        TF1   *fcnBck;
        TF1   *fe2c;

        TF1   *fcnShelfMn1;
        TF1   *fcnShelfTi1;
	TF1   *fcnShelfZr1;

	TF1   *fcnBrems;

      //  TLine *lCaka1,  *lCakb;
        TLine *lTika1,  *lTikb;
        TLine *lMnka1,  *lMnkb;
        TLine *lCuka1,  *lCukb;
        TLine *lZrka1,  *lZrkb;
	TLine *lEscZrka1, *lEscMnka1, *lEscCuka1, *lEscTika1;
        TLine *lFe,     *lMn;

        TCanvas *c1;
        TCanvas *ctmp;
        TCanvas *cev; 

        ////// Variables for fit //////
        Int_t  status;      // Status of Minimization in fit.
        Int_t  peakN;     // numbers of peaks found by TSpectrum in searchpeak()
        Int_t  rebin;    // histogram for fitting the spectrum is rebinned from channel 0 - 4096 (one channel per bin) -> to "rebin" less bins
        Int_t  ndf;
        Short_t  adc[NPADC_CH];

        AnaPar chi2;
        bool   k2Peaks; // saying if 2 peaks were found in searchpeak()
        AnaPar ll;  // lower channel limit hcal + searchpeak histograms, set in gethistorange() -> mainly limits for the hcal histogram, which is the one with the fit!!, 350 atm
        AnaPar ul;  // ... 1500 atm
        AnaPar peakX[2]; // starting value for the peak fitting algorithm found by tspectrum in searchpeak(); peakX[0,1] ... 2 biggest peaks
        AnaPar peakY[2];
        AnaPar bpar[3];
        Double_t pX[6], pY[6], pXEr[6] = {0}, pYEr[6], pXmn[6], pYmn[6], pXmnEr[6] = {0}, pYmnEr[6], pX_line[2], pY_line[2]; // fields for drawing the ev2channel line fit in
        AnaPar pXcaka1,    pXtika1,   pXcuka1,   pXmnka1, pXzrka1;
        AnaPar pXcaka1Er,  pXtika1Er, pXcuka1Er, pXmnka1Er, pXzrka1Er;
        AnaPar pXcakb13,   pXtikb1, pXcukb1, pXmnkb1, pXzrkb1;
        AnaPar pXcakb13Er, pXtikb1Er, pXcukb1Er, pXmnkb1Er, pXzrkb1Er;
        AnaPar pXfe,      mnD,       mnDEr;
        AnaPar cukaG,   mnkaG,   tikaG,  feG;
        AnaPar cukaGER, mnkaGEr, tikaGEr;
        AnaPar bck1;                                 // Constant part of back func
        AnaPar bck2;                                 // y-axis shift of back func
        AnaPar bck3;                                 // x coefficient of back func
        AnaPar bck1Er;
        AnaPar bck2Er;
        AnaPar bck3Er;
        AnaPar fano;
        AnaPar fanoEr;
        AnaPar cstn;
        AnaPar cstnEr;
        AnaPar ttiG;                                 // Gain of Ti Tail 
        AnaPar ttiGEr;
        AnaPar tmnG;                                 // Gain of Mn Tail
        AnaPar tmnGEr;                        
        AnaPar tcuG;                                 // Gain of Cu Tail
        AnaPar tcuGEr;
        AnaPar tTiBeta;                                // Beta factor of tail 
        AnaPar tTiBetaEr;
        AnaPar tMnBeta;                                // Beta factor of tail 
        AnaPar tMnBetaEr;
        AnaPar pileG;                                // Gain of Pileup 
        AnaPar pileGEr;
        AnaPar pSigR;                                // Broadening factor for Pile-up sigma
        AnaPar pSigREr;       
        AnaPar tshift;                               // Tail mean shift from Ka1
        AnaPar tshiftEr;
        AnaPar rate;                                 // Avarage rate in one run
        AnaPar shelfG;                               // gain of the shelf function


        /// Parameters for tail analysis ///
        AnaPar rtiT,     rmnT;
        AnaPar rtiTEr,   rmnTEr;
        AnaPar rcuT;
        AnaPar rcuTEr;
        AnaPar ntiT,     nmnT;
        AnaPar ntiTEr,   nmnTEr;
        AnaPar ntiK,     nmnK;
        AnaPar ntiKEr,   nmnKEr;
        AnaPar ncuT;
        AnaPar ncuTEr;
        AnaPar ncuK;
        AnaPar ncuKEr;

        AnaPar lbtl[2];
        AnaPar ubtl[2];
        AnaPar lbpk[2];
        AnaPar ubpk[2];

        ////// Variables for files //////
        TString  file_list;
        TString  rootfile;
        TString  logfolder;
        TString  slowlogfile;
        
        FILE*    flist; 
        FILE*    fslow;
        Int_t    nfile;

        ID bus;
        ID chip;
        ID id;
        ID grp;
        ID row;
        ID col;
        ID pos;                                         // pos : layer of SDD
        ID rnk;                                         // rnk : 0 ~ 5, 0 not working, 5 best
        ID padc_ch;                                     // channel number at PADC, 0-15

        Time tlabel;  
        Time year;
        Time month;
        Time day;
        Time hra;
        Time hrb;
        Time mina;
        Time minb;
        Time ut;
        Time tstart;
        Time tstop;
        Time tdiff;

        AnaPar pfit[maxPar]; // array of the fit parameters
        AnaPar sig_tika;
        AnaPar sig_tikaEr;
        AnaPar sig_mnka;
        AnaPar sig_mnkaEr;
        AnaPar fwhm;                                  // calculated at 6 keV
        AnaPar fwhmEr;                                // Error of calculated fwhm
        AnaPar ev2ch;          // ev2ch ratio from the fit of the ev2ch line in the fite2chLine() function; [ch/eV]
        AnaPar ev2chEr; 
        AnaPar offset;
        AnaPar offsetEr;

        /*
        SlowPar  vrefPOS;
        SlowPar  vrefPULSE;
        SlowPar  dacOffset;
        SlowPar  vsa;
        SlowPar  vsss;
        SlowPar  rxA;
        SlowPar  rxB;
        SlowPar  bfA;
        SlowPar  bfB;
        SlowPar  vb1A;
        SlowPar  vb1B;
        SlowPar  vb2A;
        SlowPar  vb2B;
        SlowPar  vb3A;
        SlowPar  vb3B;
        SlowPar  r1A;
        SlowPar  r1B;
        SlowPar  p12VA;
        SlowPar  vdd;
        SlowPar  m12VA;
        SlowPar  vss;
        SlowPar  p12VB;
        SlowPar  hv;
        SlowPar  m12VB;
        SlowPar  m5V;
        SlowPar  tsenseA;
        SlowPar  tsenseB;
        SlowPar  pt100;
        SlowPar  tADC;
        SlowPar  xtubeAm;

        SlowPar  current;                       // The constant current set for the target
        */

};
#endif
