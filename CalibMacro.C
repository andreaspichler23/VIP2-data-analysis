
#include  <stdlib.h>
#include  <stdio.h>
#include  <iostream>
#include  <errno.h>

#include  "SDDclass.h"
#include  "SDDclass.C"
#include  "SDDcalib.h"
#include  "common.h"


using namespace std;


void CalibMacro( ID nSDD,  ID nPADC,  TString rootf, TString source, TString place )
{

    Bool_t   saveflag = kTRUE,  saveplot = kFALSE;

    //ofstream logfile;
    ofstream datafile;
    //logfile.open(WORK_PATH + "/reports/Analysis/logifle.txt",ios::app);
    //datafile.open(WORK_PATH + "/reports/Analysis/BestTimeForRootfile/01day-datafile.txt",ios::app);
    datafile.open(WORK_PATH + "/calibrationlist/1-618Files-LNGS.txt",ios::app);
    TString fileNamePart = rootf(20,3); // change here!! 
    Int_t fileNumber = fileNamePart.Atoi();

    Int_t    fstatus,  ndf;  
    /*
    Int_t    bus, id,  fstatus,  ndf,  rnk;
    Int_t    year,  month, day;
    */
    Float_t  fano,  fanoEr,  cstn,  cstnEr,  sig_tika,  sig_tikaEr,  fwhm, fwhmEr; 
    Float_t  backG, cukaG;
    Float_t  tika1, tika1Er, mnka1, mnka1Er, cuka1, cuka1Er;
    Float_t  tikay, mny, ttiG, pG, pGEr;
    Float_t  ev2ch,  ev2chEr, offset,  offsetEr, chi2, tBeta, tshift, slope, slopeEr;

    SDDclass::SDDclass* sdd = new SDDclass::SDDclass( nSDD, nPADC );

    sdd->SetRootfile( rootf, place );
    sdd->OpenCanvas();

    if( sdd->SearchPeak() ){
        if(source == "TiMnCu") {fstatus = sdd->FitTiMnCu( saveflag );}
        if(source == "TiMn")   {fstatus = sdd->FitTiMn( saveflag );}
	if(source == "TiCuZr") {fstatus = sdd->FitTiCuZr( saveflag );}
        sdd->GetFitParameters( source );
        sdd->PlotResidue(source);
        sdd->CalcFwhmMn(source);
        sdd->FitE2ChLine( source );
        
        backG   = (Float_t)sdd->GetBackG();
        cukaG   = (Float_t)sdd->GetCuKaG();
        fano    = (Float_t)sdd->GetFano();
        fanoEr  = (Float_t)sdd->GetFanoEr();
        cstn    = (Float_t)sdd->GetCstn();
        cstnEr  = (Float_t)sdd->GetCstnEr();
        sig_tika   = (Float_t)sdd->GetSigTika();
        sig_tikaEr = (Float_t)sdd->GetSigTikaEr();
        fwhm    = (Float_t)sdd->GetFwhm();
        fwhmEr    = (Float_t)sdd->GetFwhmEr();
        tika1   = (Float_t)sdd->GetTika1();
        tika1Er = (Float_t)sdd->GetTika1Er();
        mnka1   = (Float_t)sdd->GetMnka1();
        mnka1Er = (Float_t)sdd->GetMnka1Er();
        cuka1   = (Float_t)sdd->GetCuka1();
        cuka1Er = (Float_t)sdd->GetCuka1Er();
        //mnDt    = (Float_t)sdd->GetMnD();
        //mnDtEr  = (Float_t)sdd->GetMnDEr();
        tikay   = (Float_t)sdd->GetTiKaG();
        mny     = (Float_t)sdd->GetMnKaG();
        ttiG    = (Float_t)sdd->GetTiTailG();
        //pG      = (Float_t)sdd->GetPileG();
        //pGEr    = (Float_t)sdd->GetPileGEr();
        ev2ch   = (Float_t)sdd->GetEv2ch();
        ev2chEr = (Float_t)sdd->GetEv2chEr();
        offset  = (Float_t)sdd->GetOffset();
        offsetEr= (Float_t)sdd->GetOffsetEr();
        slope =   (Float_t)sdd->GetSlope();
        slopeEr = (Float_t)sdd->GetSlopeEr();
        chi2    = (Float_t)sdd->GetChi2();
        tshift  = (Float_t)sdd->GetTShift();
        tBeta   = (Float_t)sdd->GetTTiBeta();
        ndf     = sdd->GetNdf();
 
//        if( fstatus != 0 ){
//       
//            //logfile << "Fit error in file " << rootf << " with SDD " << nSDD  << " status: " << fstatus << endl;
//            cout << endl << "Fit error in file " << rootf << " with SDD " << nSDD  << " status: " << fstatus << endl;
//            
//        } 
        
        // this is for BestTimeForRootfile: datafile << nSDD << " " << fwhm << " " << fwhmEr << " " << tika1 << " " << tika1Er << " " << mnka1 << " " << mnka1Er << " " << cuka1 << " " << cuka1Er << " " << fstatus << endl;
        //cout << endl << " THE FIT STATUS IS: " << fstatus << endl;
        datafile << fileNumber << " " << nSDD << " " << slope << " " << slopeEr << " " << offset << " " << offsetEr << " " << fwhm << " " << fwhmEr << " " << fano << " " << fanoEr << " " 
                << cstn << " " << cstnEr << " " << mnka1 << " " << mnka1Er << " " << mny << " " << tika1 << " " << tika1Er << " " << tikay << " " << cuka1 << " " << cuka1Er << " " << backG << " "
                << cukaG << " " << chi2 << " " << fstatus << " " << ndf << endl;
        //logfile.close();
        datafile.close();
 
    }else{
        fstatus = 9;
        /*
        tikay = 0.,   ttiG  = 0.;  
        pG      = 0.,  pGEr  = 0.,   fano = 0.,  fanoEr  = 0.,  cstn  = 0.,  cstnEr = 0.; 
        mny     = 0.,  sig_tika   = 0.,          sig_tikaEr = 0.,             fwhm  = 0.; 
        tika1   = 0.,  tika1Er = 0.; 
        mnka1   = 0.,  mnka1Er = 0.,             ev2ch   = 0.,  ev2chEr = 0.; 
        mnDt    = 0.,  mnDtEr  = 0.;
        offset  = 0.,  offsetEr= 0., chi2 = 0.,  tshift  = 0.,  tBeta = 0.,  ndf    = 0;
        */
    }

    sdd->CloseCanvas( saveplot );
    if( fstatus != 9 ){              // When histogram fit was done.
        sdd->CloseRootFile();   
    }
    else{
        sdd->CloseCanvas( false );
    }
    delete sdd;

    cout << "# END of one run for one SDD" << endl;

    return; 
}