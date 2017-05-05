/*****************************************
               SDDStability.C
     Make root file for SDD parameters 
     including calibration data with 
     slow log .
             July 2009
             H. Shi
******************************************/
#include  <stdlib.h>
#include  <stdio.h>
#include  <iostream>
#include  <errno.h>
#include  "SDDtable.h"
#include  "SDDclass.C"
#include  "TTree.h"
#include  "TBranch.h"
#include  "SDDStability.h"
using namespace std;


void SDDProfile( TString filelist, 
                 Bool_t saveflag, Bool_t saveplot ){
    Int_t    nf = 0;
    Int_t    count = 0;
    TString  ith_file;
    TString fname;
    nf = GetNumberOfRuns( filelist );
    fname = filelist;
    char    f_name[9];
    sscanf( fname, "%*19s%9s\n", f_name );
    TString  parfile = Form( CALPAR_PATH + "/txt/%s_log.txt", f_name );
    for( Int_t ibus = 1; ibus <= NBUS; ibus ++ ){
        for( Int_t isdd = 0; isdd < NSDD; isdd ++ ){
            if( RnkSDD[ibus][isdd] >= 1 ){ //&& RowSDD[ibus][isdd] <=2 ){ //  && ibus == 3 && isdd == 11 ){
                for( Int_t ifile = 1; ifile <= nf; ifile ++ ){
                    ith_file = GetRunNames( filelist, ifile );
                    SingleRunOneSDDCalib(ibus, isdd, filelist, nf, ith_file, 
                            saveflag, saveplot, parfile ); 
                    count ++;
                    cout << "# Total number of calibration : " << count << endl;
                }
            }
        }
    }
    MakeParRootFile( parfile,  filelist,  saveflag );
    return;
}


void SingleRunOneSDDCalib( ID nBus,  ID nSDD,   TString input_list,  
                           Int_t nf, TString    ith_run,  
                           Bool_t  saveflag,    Bool_t saveplot,
                           TString parfile){
    Int_t     bus, id,  chip,   grp,  row,  col,  pos,  fstatus,  ndf,  rnk;
    UInt_t    tlabel,   ut;
    Float_t  fano,  fanoEr,  cstn,  cstnEr,  sig_tika,  sig_tikaEr,  fwhm; 
    Float_t  tika1, tika1Er, cuka1, cuka1Er, mnka1, mnka1Er, mnDt, mnDtEr;
    Float_t  cukay, tikay, fey, mny, rcuti, ttiG, tcuG, pG, pGEr;
    Float_t  ev2ch,  ev2chEr, offset,  offsetEr, chi2, tBeta, tshift;
    Float_t  vrefPOS,  vrefPULSE,  dacOffset,  vsa,  vsss,  rxA,  rxB;
    Float_t  bfA,  bfB,  vb1A,  vb1B,  vb2A,  vb2B,  vb3A,  vb3B;
    Float_t  r1A,  r1B,  p12VA,  vdd,  m12VA,  vss, p12VB,  hv,  m12VB;
    Float_t  m5V,  tsenseA,  tsenseB,  pt100,  tADC, rate,  xAm;  

    SDDclass::SDDclass* sdd = new SDDclass::SDDclass( nBus, nSDD );
    if( sdd->SetReadList( input_list,  nf ) ){
        sdd->SetRootfile( ith_run );
        sdd->SetSlowPar();
        //sdd->GetRunTime();
        sdd->OpenCanvas();
        // Get SDD labels
        bus     = sdd->GetBus();
        id      = sdd->GetID();
        chip    = sdd->GetChip();
        grp     = sdd->GetGrp();
        pos     = sdd->GetPos();
        row     = sdd->GetRow();
        col     = sdd->GetCol();
        ut      = sdd->GetUnixT();
        tlabel  = sdd->GetTlabel();
        rnk     = sdd->GetRnk();
        if( sdd->SearchPeak() && rnk > 1 ){
            //fstatus = sdd->FitTiCu( saveflag );
            fstatus = sdd->FitTiMn( saveflag );
            //fstatus = sdd->FitTiMnCu( saveflag );
            sdd->GetFitParameters();
            rcuti   = sdd->GetCuTi();
            sdd->PlotResidue();
            sdd->CalcSigTi();
            sdd->FitE2ChLine();

            fano    = (Float_t)sdd->GetFano();
            fanoEr  = (Float_t)sdd->GetFanoEr();
            cstn    = (Float_t)sdd->GetCstn();
            cstnEr  = (Float_t)sdd->GetCstnEr();
            sig_tika   = (Float_t)sdd->GetSigTika();
            sig_tikaEr = (Float_t)sdd->GetSigTikaEr();
            fwhm    = (Float_t)sdd->GetFwhm();
            tika1   = (Float_t)sdd->GetTika1();
            tika1Er = (Float_t)sdd->GetTika1Er();
            cuka1   = (Float_t)sdd->GetCuka1();
            cuka1Er = (Float_t)sdd->GetCuka1Er();
            mnka1   = (Float_t)sdd->GetMnka1();
            mnka1Er = (Float_t)sdd->GetMnka1Er();
            mnDt    = (Float_t)sdd->GetMnD();
            mnDtEr  = (Float_t)sdd->GetMnDEr();
            tikay   = (Float_t)sdd->GetTiKaG();
            fey     = (Float_t)sdd->GetFeG();
            mny     = (Float_t)sdd->GetMnKaG();
            cukay   = (Float_t)sdd->GetCuKaG();
            ttiG    = (Float_t)sdd->GetTiTailG();
            tcuG    = (Float_t)sdd->GetCuTailG();
            pG      = (Float_t)sdd->GetPileG();
            pGEr    = (Float_t)sdd->GetPileGEr();
            ev2ch   = (Float_t)sdd->GetEv2ch();
            ev2chEr = (Float_t)sdd->GetEv2chEr();
            offset  = (Float_t)sdd->GetE2ChOffset();
            offsetEr= (Float_t)sdd->GetE2ChOffsetEr();
            chi2    = (Float_t)sdd->GetChi2();
            tshift  = (Float_t)sdd->GetTShift();
            tBeta   = (Float_t)sdd->GetTBeta();
            ndf     = sdd->GetNdf();
        }else{
            fstatus = 9;
            cukay   = 0.,  tikay = 0.,   fey  = 0.,  rcuti   = 0.,  ttiG  = 0.,  tcuG   = 0.;
            pG      = 0.,  pGEr  = 0.,   fano = 0.,  fanoEr  = 0.,  cstn  = 0.,  cstnEr = 0.; 
            mny     = 0.,  sig_tika   = 0.,          sig_tikaEr = 0.,             fwhm  = 0.; 
            tika1   = 0.,  tika1Er = 0.,             cuka1   = 0.,  cuka1Er = 0.; 
            mnka1   = 0.,  mnka1Er = 0.,             ev2ch   = 0.,  ev2chEr = 0.; 
            mnDt    = 0.,  mnDtEr  = 0.;
            offset  = 0.,  offsetEr= 0., chi2 = 0.,  tshift  = 0.,  tBeta = 0.,  ndf    = 0;
        }
        // Get Slow parameters
        vrefPOS   = (Float_t)sdd->GetVrefPOS();
        vrefPULSE = (Float_t)sdd->GetVrefPULSE(); 
        dacOffset = (Float_t)sdd->GetDacOffset();
        vsa       = (Float_t)sdd->GetVsa();
        vsss      = (Float_t)sdd->GetVsss();
        rxA       = (Float_t)sdd->GetRxA();
        rxB       = (Float_t)sdd->GetRxB();
        bfA       = (Float_t)sdd->GetBfA();
        bfB       = (Float_t)sdd->GetBfB();
        vb1A      = (Float_t)sdd->GetVb1A();
        vb1B      = (Float_t)sdd->GetVb1B();
        vb2A      = (Float_t)sdd->GetVb2A();
        vb2B      = (Float_t)sdd->GetVb2B();
        vb3A      = (Float_t)sdd->GetVb3A();
        vb3B      = (Float_t)sdd->GetVb3B();
        r1A       = (Float_t)sdd->GetR1A();
        r1B       = (Float_t)sdd->GetR1B();
        p12VA     = (Float_t)sdd->GetP12VA();
        vdd       = (Float_t)sdd->GetVdd();
        m12VA     = (Float_t)sdd->GetM12VA();
        vss       = (Float_t)sdd->GetVss();
        p12VB     = (Float_t)sdd->GetP12VB();
        hv        = (Float_t)sdd->GetHv();
        m12VB     = (Float_t)sdd->GetM12VB();
        m5V       = (Float_t)sdd->GetM5V();
        tsenseA   = (Float_t)sdd->GetTsenseA();
        tsenseB   = (Float_t)sdd->GetTsenseB();
        pt100     = (Float_t)sdd->GetPt100();
        tADC      = (Float_t)sdd->GetTADC();
        rate      = (Float_t)sdd->GetAverageRate();
        xAm       = (Float_t)sdd->GetXtubeAm();

        // Write parameter data to par file
        FILE *fpar = fopen( parfile, "a");
        Int_t  rtn_close,  rtn_fprintf;
        //cout << "# fpar in SingleRunOneSDDCalib() : " << fpar << endl;
        if( fpar != NULL ){
            fflush( fpar );
            if( saveflag ){
            rtn_fprintf = 
                fprintf( fpar, "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d,"
                               "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
                               "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
                               "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
                               "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
                               "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
                        bus, id, chip, grp, pos, row, col, rnk,  ut, tlabel, fstatus, ndf, 
                        fano, fanoEr, cstn, cstnEr, sig_tika, sig_tikaEr, fwhm, 
                        tika1, tika1Er, cuka1, cuka1Er, mnka1, mnka1Er, mnDt, mnDtEr,
                        cukay, tikay, fey, mny, rcuti, ttiG, tcuG, pG, pGEr,
                        ev2ch, ev2chEr, offset, offsetEr, chi2, tBeta, tshift,
                        vrefPOS, vrefPULSE, dacOffset, vsa, vsss, rxA, rxB, bfA, bfB, 
                        vb1A, vb1B, vb2A, vb2B, vb3A, vb3B, r1A, r1B, 
                        p12VA, vdd, m12VA,vss,  p12VB,hv, m12VB, m5V, 
                        tsenseA, tsenseB, pt100, tADC, rate, xAm );
            }
            rtn_close = fclose(fpar);
            //cout << "# Return value of fprintf in SingleRunOneSDDCalib() : " << rtn_fprintf << endl;
            //cout << "# Return value of fclose(fpar) in SingleRunOneSDDCalib() : " << rtn_close << endl;
        }else{
            perror("# The following error occurred");
            cerr << "# Could not open fpar in SingleRunOneSDDCalib() ... " << endl;
        }
        sdd->CloseCanvas( saveplot );
        if( fstatus != 9 ){              // When histogram fit was done.
            sdd->CloseRootFile();
        }
        //else{
        //sdd->CloseCanvas( false );
        //}
    }
    delete sdd;
    cout << "# END of one run for one SDD" << endl;
    return; 
}


void MakeParRootFile( TString  parfile, TString filelist , Bool_t saveflag){
    Int_t     bus,  chip,  id,  grp,  row,  col,  pos, fstatus, ndf,  rnk;
    UInt_t    tlabel,  ut;
    Float_t   fano,  fanoEr,  cstn,  sig_tika,  sig_tikaEr,  fwhm,  cstnEr; 
    Float_t   tika1, tika1Er, cuka1, cuka1Er, mnka1, mnka1Er, mnDt, mnDtEr;;
    Float_t   cukay, tikay, fey, mny, rcuti, ttiG, tcuG, pG, pGEr;
    Float_t  ev2ch,  ev2chEr, offset,  offsetEr, chi2, tBeta, tshift;
    Float_t  vrefPOS,  vrefPULSE,  dacOffset,  vsa,  vsss,  rxA,  rxB;
    Float_t  bfA,  bfB,  vb1A,  vb1B,  vb2A,  vb2B,  vb3A,  vb3B;
    Float_t  r1A,  r1B,  p12VA,  vdd,  m12VA,  vss, p12VB,  hv,  m12VB;
    Float_t  m5V,  tsenseA,  tsenseB,  pt100,  tADC, rate,  xAm;  

    TString fname;
    fname = filelist;
    char    f_name[9];
    sscanf( fname, "%*19s%9s\n", f_name );
    TString rootfile = Form( CALPAR_PATH + "/root_file/%s_sdd_profile.root", f_name );
    TFile   *f  = new TFile( rootfile, "RECREATE");
    f->cd();
    TTree *tr   = new TTree("tr", "tree of parameters");
    tr->Branch("bus",    &bus, "bus/I");
    tr->Branch("id",     &id,  "id/I");
    tr->Branch("chip",   &chip,"chip/I");
    tr->Branch("row",    &row, "row/I");
    tr->Branch("col",    &col, "col/I");
    tr->Branch("grp",    &grp, "grp/I");
    tr->Branch("pos",    &pos, "pos/I");
    tr->Branch("rnk",    &rnk, "rnk/I");
    tr->Branch("ut",     &ut,  "unixT/I");
    tr->Branch("tlabel", &tlabel, "time/I");
    tr->Branch("fstatus",&fstatus, "fStatus/I");
    tr->Branch("chi2",   &chi2,    "chi2/F");
    tr->Branch("ndf",    &ndf,     "ndf/I");
    tr->Branch("fano",   &fano, "fano/F");
    tr->Branch("fanoEr", &fanoEr,"fanoEr/F");
    tr->Branch("cstn",   &cstn, "cstn/F");
    tr->Branch("cstnEr", &cstnEr,"cstnEr/F");
    tr->Branch("sig_tika",&sig_tika,"sig_tika/F");
    tr->Branch("fwhm",   &fwhm,  "fwhmTiKa/F");
    tr->Branch("tika1",  &tika1, "tika1/F");
    tr->Branch("tika1Er",&tika1Er, "tika1Er/F");
    tr->Branch("cuka1",  &cuka1, "cuka1/F");
    tr->Branch("cuka1Er",&cuka1Er,"cuka1Er/F");
    tr->Branch("mnka1",  &mnka1,  "mnka1/F");
    tr->Branch("mnka1Er",&mnka1Er,"mnka1Er/F");
    tr->Branch("mnDt",   &mnDt,   "mnDt/F");
    tr->Branch("mnDtEr", &mnDtEr, "mnDtEr/F");
    tr->Branch("cukay",  &cukay,  "cukay/F");
    tr->Branch("tikay",  &tikay,  "tikay/F");
    tr->Branch("fey",    &fey,    "fey/F" );
    tr->Branch("mny",    &mny,    "mny/F" );
    tr->Branch("rcuti",  &rcuti,  "CuTiRatio/F");
    tr->Branch("ttiG",   &ttiG,   "TiTailGain/F");
    tr->Branch("tcuG",   &tcuG,   "CuTailGain/F");
    tr->Branch("pG",     &pG,     "PileUpGain/F");
    tr->Branch("pGEr",   &pGEr,   "PileGainEr/F");
    tr->Branch("tBeta",  &tBeta,  "tailBeta/F");
    tr->Branch("tshift", &tshift, "tailShift/F");
    tr->Branch("ev2ch",  &ev2ch,  "ev2ch/F");
    tr->Branch("ev2chEr",&ev2chEr,"ev2chEr/F");
    tr->Branch("offset", &offset, "offset/F");
    tr->Branch("offsetEr",&offsetEr,"offsetEr/F");
    tr->Branch("vrefPOS", &vrefPOS, "vrefPULSE/F");
    tr->Branch("vrefPULSE", &vrefPULSE, "vrefPULSE/F");
    tr->Branch("dacOffset", &dacOffset, "dacOffset/F");
    tr->Branch("vsa",     &vsa,    "vsa/F");
    tr->Branch("vsss",    &vsss,   "vsss/F");
    tr->Branch("rxA",     &rxA,    "rxA/F");
    tr->Branch("rxB",     &rxB,    "rxB/F");
    tr->Branch("bfA",     &bfA,    "bfA/F");
    tr->Branch("bfB",     &bfB,    "bfB/F");
    tr->Branch("vb1A",    &vb1A,   "vb1A/F");
    tr->Branch("vb1B",    &vb1B,   "vb1B/F");
    tr->Branch("vb2A",    &vb2A,   "vb2A/F");
    tr->Branch("vb2B",    &vb2B,   "vb2B/F");
    tr->Branch("vb3A",    &vb3A,   "vb3A/F");
    tr->Branch("vb3B",    &vb3B,   "vb3B/F");
    tr->Branch("r1A",     &r1A,    "r1A/F");
    tr->Branch("r1B",     &r1B,    "r1B/F");
    tr->Branch("p12VA",   &p12VA,  "p12VA/F");
    tr->Branch("vdd",     &vdd,    "vdd/F");
    tr->Branch("m12VA",   &m12VA,  "m12VA/F");
    tr->Branch("vss",     &vss,    "vss/F");
    tr->Branch("p12VB",   &p12VB,  "p12VB/F");
    tr->Branch("hv",      &hv,     "hv/F");
    tr->Branch("m12VB",   &m12VB,  "m12VB/F");
    tr->Branch("m5V",     &m5V,    "m5V/F");
    tr->Branch("tsenseA", &tsenseA,"tsenseA/F");
    tr->Branch("tsenseB", &tsenseB,"tsenseB/F");
    tr->Branch("pt100",   &pt100,  "pt100/F");
    tr->Branch("tADC",    &tADC,   "tADC/F");
    tr->Branch("rate",    &rate,   "rate/F");
    tr->Branch("xAm",     &xAm,    "Xtube/F");

    FILE *fpar = fopen( parfile, "r");
    //cout << "# fpar in MakeRootFile() : " << fpar << endl;
    if( fpar != NULL ){
        char   parline[MAXCHAR];
        while( fgets( parline, MAXCHAR, fpar ) != NULL ){
            //cout << "# Parline to read : " << endl << "     " << parline << endl;
            sscanf( parline, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
                             "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
                    &bus, &id, &chip, &grp, &pos, &row, &col, &rnk, &ut, &tlabel, &fstatus, &ndf, 
                    &fano, &fanoEr, &cstn, &cstnEr, &sig_tika, &sig_tikaEr, &fwhm, 
                    &tika1, &tika1Er, &cuka1, &cuka1Er, &mnka1, &mnka1Er, &mnDt, &mnDtEr,
                    &cukay, &tikay, &fey, &mny, &rcuti, &ttiG, &tcuG, &pG, &pGEr, 
                    &ev2ch, &ev2chEr, &offset, &offsetEr, &chi2, &tBeta, &tshift, 
                    &vrefPOS, &vrefPULSE, &dacOffset, &vsa, &vsss, &rxA, &rxB, &bfA, &bfB, 
                    &vb1A, &vb1B, &vb2A, &vb2B, &vb3A, &vb3B, &r1A, &r1B, 
                    &p12VA, &vdd, &m12VA, &vss, &p12VB, &hv, &m12VB, &m5V, 
                    &tsenseA, &tsenseB, &pt100, &tADC, &rate, &xAm );
            cout << "# fano, cstn, sig_tika, tika1: " 
                <<  fano << "   " << cstn << "    " << sig_tika << "   " 
                <<  tika1 << endl;
            tr->Fill();
        }
        fclose(fpar);
    }else{
        perror("# The following error occurred");
        cerr << "# Could not open par file in MakeRootFile() ... " << endl;
    }
    if( saveflag ){
        f->Write();
    }
    f->Close();
    cout << "#Made root file for SDD parameter data at " << rootfile << endl;
    return;
}


Int_t GetNumberOfRuns( TString filelist ){
    Int_t nfile = 0;
    Int_t rtn_close;
    char  fline[MAXCHAR];
    FILE* flist = fopen( filelist, "r" );
    cout << "# flist in GetNumberOfRuns() : " << flist << endl;
    if( flist != NULL ){
        while( fgets( fline, MAXCHAR, flist ) != NULL ){
            nfile ++ ;
        }
        rtn_close = fclose(flist);
        //cout << "# Return value of fclose(flist) in GetNumberOfRuns() : " << rtn_close << endl;
    }else{
        perror("# The following error occurred");
        cerr << "# Could not open flist in GetNumberOfRuns() ... " << endl;
    }
    return nfile;
}

TString GetRunNames( TString filelist,  Int_t ifile ){
    TString  ith_run_name;
    Int_t    iline = 1;
    Int_t    rtn_close;
    char  fline[MAXCHAR];
    FILE* fl = fopen( filelist,  "r" );
    cout << "# fl in GetRunNames() : " << fl << endl;
    if( fl != NULL ){
        while( fgets( fline, MAXCHAR, fl ) != NULL ){
            if( iline == ifile ){
                ith_run_name = fline; 
                break;
            }
            iline ++;
        }
        rtn_close = fclose(fl);
        //cout << "# Return value of fclose(fl) in GetRunNames() : " << rtn_close << endl;
    }else{
        perror("# The following error occurred");
        cout << "# In SDDStability.C, function GetRunNames(), list file not properly opened. " 
            << endl;
    }
    return ith_run_name; 
}
