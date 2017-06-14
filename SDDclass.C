/**************************************************************
              SDDclass.C

                  Jun 23 2009 
                  H. Shi
                  For SIDDHARTA SDDs

                  2009/12/30
                  Added fwhmEr with only propagation from slopeEr;

                  2013/09/25 
                  First Modification for VIP2 usage; 

                  2014/05/16
                  Second modification for VIP2 test setup;

                  2014/08/26
                  For VIP2 SDD calibration. 

                 July 2016: modified again for VIP-2 use by Andreas Pichler
**************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <iostream>

#include  "TAxis.h"
#include  "TCanvas.h"
#include  "TBranch.h"
#include  "TDatime.h"
#include  "TFile.h"
#include  "TGraph.h"
#include  "TGraphErrors.h"
#include  "TH1.h"
#include  "TSpectrum.h"
#include  "TLatex.h"
#include  "TStyle.h"
#include  "TString.h"
#include  "TTree.h"

#include  "common.h"
#include  "PhaseOneSDDConnectionTable.h"
#include  "XrayLines.h"
#include  "SDDclass.h"
#include  "CalibFunction.h"
#include  "CalibFunction.C"
#include  "errcalc.C"

using namespace std;

SDDclass::SDDclass( ID numSDD, ID adcCH )
{
    // cout < "# Constructor of SDDclass. " << endl;

    ///// ID definition ////////
   // bus     = numBus;
    id      = numSDD;
    padc_ch = adcCH;


    cout << Form("# Implement SDD number %d at ADC channel  %d.", id, padc_ch ) << endl;
}


SDDclass::~SDDclass(){
    //cout << "# Destructor of SDDclass. " << endl;
}


//**** Member functions ****// 



void SDDclass::SetRootfile( TString run_name, TString place ) // name should be in format year-month-day
{   
    cout << "# In SDDclass::SetRootfile() ... " << endl;

    Int_t size = place.Sizeof();

    if(size==4){rootfile =  ROOT_PATH_SMI + "/" + run_name;}
    if(size==5){rootfile =  ROOT_PATH_LNGS + "/" + run_name;}

  //  year = 2014,  month = 8; day = 21,  hra = 0,  mina = 0;

    hra = 0;
    mina = 0;
    sscanf( run_name, "%04d%02d%02d\n", 
                      &year, &month, &day );

    tlabel = day * 1000000 + month * 10000 + year;
    cout << "# Calibration file : " << rootfile << endl;
    TDatime daT(year, month, day, hra, mina, 0);
    ut = daT.Convert();
//    Int_t date = daT.GetDate();
    cout << "# tlabel,  ut  : " << tlabel << "   " << ut << endl << endl;

    return ;
}


bool SDDclass::SearchPeak()
{
    cout << "# In SDDclass::SearchPeak() ... " << endl;
    cout << "# Rootfile : " << rootfile << endl;
    k2Peaks = false;
    Double_t  thr = 0.01;
    Int_t     lb  = 400; 
    peakX[0] = 0., peakX[1] = 0., peakY[0] = 0., peakY[1] = 0.;

    for(Int_t i=0; i<NPADC_CH; i++){
        adc[i] = 0;
    }


    this->GetHistoRange();
           
    TFile *ftmp = new TFile( rootfile, "READ");
    ftmp->cd();
    TTree *tr_t = (TTree*)ftmp->Get("tr");
    Int_t ent = (Int_t)tr_t->GetEntries();
    cout << "# entries in tree: " << ent << endl;



    TH1F *ht = new TH1F("ht", "ht", 1024, -0.5, 4095.5); // DONT CHANGE BINNING HERE

   
  
    TBranch *adcEvent   = tr_t->GetBranch("adc"); 
    adcEvent->SetAddress(adc);

    

    for( Int_t ient=0; ient<ent; ient ++){
        adcEvent->GetEntry(ient);
        if( adc[padc_ch] > lb){
            ht->Fill((Double_t)(adc[padc_ch]));
        }
    }

    ht->Draw();
   

    TH1F *htmp = (TH1F*)ht->Clone("htmp");
    //TH1F *hb   = (TH1F*)ht->Clone("hb");
    htmp->SetStats(0);


   // ctmp->cd(1);
    //htmp->Rebin(8);


    htmp->GetXaxis()->SetRangeUser( ll, ul );
    //htmp->SetAxisRange( ll, ul, "X");

  

    TSpectrum *s = new TSpectrum();
    s->Search(htmp, 1, "", thr );                    // 1 sigma, 0.01 threshold
 

    peakN = s->GetNPeaks();
    Double_t *pX = s->GetPositionX(); // pX is the channel number of the peaks starting from the highest peak at pX[0] - for Fe source LNGS data this is always Ti
// for X ray tube data - this can be cu or ti
    Double_t *pY = s->GetPositionY();



    cout << "p1 : " << pX[0] << " p2 : " << pX[1] << " p3 : " << pX[2] << " p4 : " << pX[3] << " p5 : " << pX[4] << endl;


    if( peakN == 4 || peakN == 5){ // peakX[0] should be the Ti channel and peakX[1] the Mn/Cu channel
        if( pX[0] <= ul && pX[0] > ll ){
	  if( pX[0] < pX[1] ){ // if Ti peak bigger than Mn/ Cu(in the case of X ray tube)
            peakX[0] = pX[0];
            peakY[0] = pY[0];
            peakX[1] = pX[1];
            peakY[1] = pY[1];
	  }
	  else{ // if Mn / Cu peak is bigger than Ti
            peakX[0] = pX[1];
            peakY[0] = pY[1];
            peakX[1] = pX[0];
            peakY[1] = pY[0];

	  }
        }
    }

    //cout << "peakY[0] = " << peakY[0] << "peakY[1] = " << peakY[1] << endl;

    cout << "found peaks: " << peakN << "; " << endl;


    if( peakX[0] > 0 && peakX[1] > 0 ){
        k2Peaks = true;

    }else{ 
        k2Peaks = false;
        cout << "WE GOT NO PEAKS" << endl;
    }

    // Initial value for the Iron peak:
    /*
    // Ti Cu case
    Double_t sl = ( peakX[1] - peakX[0] ) / ( CuKa1 - TiKa1 );
    pXfe  = (peakX[0]+peakX[1])/2 + sl * (FeKa1-(CuKa1+TiKa1)/2);
    cout << "# Initial values for TiKa1 and CuKa1: " << peakX[0] 
            << ";  " << peakX[1] << endl;
    cout << "# Initial values for Fe mean :" << pXfe << endl;

    // Ti Mn case
    Double_t sl = ( peakX[1] - peakX[0] ) / ( MnKa1 - TiKa1 );
    pXfe  = (peakX[0]+peakX[1])/2 + sl * (FeKa1-(MnKa1+TiKa1)/2);
    cout << "# Initial values for TiKa1 and MnKa1: " << peakX[0] 
            << ";  " << peakX[1] << endl;
    cout << "# Initial values for Fe mean :" << pXfe << endl;

    //////////// Fit background excluding peak area/////////
    cout << "# Fitting initial values for background" << endl;
    TH1F *hbg = new TH1F("hbg","hbg", 3200, 0.5, 3200.5);
    Double_t  cont = 0.;
    for( Int_t i = (Int_t)ll; i < (Int_t)ul; i ++ ){
        if( (i<(Int_t)peakX[0]-100 && i>(Int_t)peakX[0]-200)  
                || (i > (Int_t)peakX[0] + 180 && i < (Int_t)pXfe - 250)
                || i > (Int_t)peakX[1] + 350.){
            cont = hb->GetBinContent(i);
            hbg->SetBinContent(i, cont);
        }
    }
    ctmp->cd(2);
    TF1    *fb  = new TF1("fb", backFunc, ll, ul , 3);
    hbg->Rebin(10);
    hbg->Fit("fb","R");
    hbg->GetXaxis()->SetRangeUser(700,3000);
    bpar[0] = fb->GetParameter(0);
    bpar[1] = fb->GetParameter(1);
    bpar[2] = fb->GetParameter(2);
    ////////////////////////////////////////////////////////
    */
    bpar[0] = 0;
    bpar[1] = 0;
    bpar[2] = 0;

    s->Delete();
 //   ctmp->Update();
    //ftmp->Close();

   
    if( k2Peaks ){
        cout << "# Two initial peaks for calib found. " << endl;
    }
    return k2Peaks;
}


Int_t SDDclass::FitTiMn( Bool_t saveflag )
{
    Int_t npar_timn = 22; // set here number of parameters for fit
    cout << "# In SDDclass::FitTiMn() ... " << endl;
    status = 0;          // Fit status:Return value of Fit. 0: converged; 4: abnormal termination
    chi2   = 0.;
    ndf    = 0;
    rebin  = 4;

    /// for tail ana /////
    ntiT   = 0.;
    ntiTEr = 0.;
    ntiK   = 0.;
    ntiKEr = 0.;
    nmnT   = 0.;
    nmnK   = 0.;
    nmnTEr = 0.;
    nmnKEr = 0.;
    rtiT   = 0.;
    rtiTEr = 0.;
    rmnT   = 0.;
    rmnTEr = 0.;
    for( Int_t i = 0; i <= 1; i++ ){
        lbtl[i] =  0.;
        ubtl[i] =  0.;
        lbtl[i] =  0.;
        ubtl[i] =  0.;
        lbpk[i] =  0.;
        ubpk[i] =  0.;
        lbpk[i] =  0.;
        ubpk[i] =  0.;
    }
    ////////////////////

    if( k2Peaks == true ){
        Double_t  l_bin = 350; // lower end of histogram
        Double_t  u_bin = peakX[1] + 300.; // 300 channels above smaller of the 2 peaks - sets range for fitting of the main fit function
        Double_t  ti_gain_setvalue = peakY[0] * rebin; // setting initial gains - but how/why so? - NUMERICAL ESTIMATION
        Double_t  mn_gain_setvalue = peakY[1] * rebin;
 //       Double_t  ca_gain_setvalue = ti_gain_setvalue / 10; // ca gain is 1/10 of the ti gain... ok?!
        Double_t  init_mean_ti     = peakX[0] ;
        Double_t  init_mean_mn     = peakX[1] ;
        Double_t  slope            = ( init_mean_mn - init_mean_ti )/( MnKa1 - TiKa1 ); // inital slope in channels/eV 
 //       Double_t  init_mean_ca     = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-CaKa1)*slope; // initial channel of Ca Kalpha
 //       Double_t  init_mean_cakb   = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-CaKb1_3)*slope;

        Double_t  init_fano        = 0.16; // initial fano factor (=variance of produced charges squared / mean number of produced charges)
        Double_t  init_cnst        = 10.;  // constant noise


        Int_t     range = (Int_t)u_bin - (Int_t)l_bin + 1; 
        //range     = ( Int_t )( range / rebin );
        Double_t  init_mean_tikb   = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-TiKb1)*slope;
        Double_t  tikb_lbound      = init_mean_tikb - 20.;
        Double_t  tikb_ubound      = init_mean_tikb + 20.;
        Double_t  init_mean_mnkb   = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-MnKb1)*slope;
        Double_t  mnkb_lbound      = init_mean_mnkb - 20.;
        Double_t  mnkb_ubound      = init_mean_mnkb + 20.;

        Double_t  mnka_lbound      = init_mean_mn - 20.;
        Double_t  mnka_ubound      = init_mean_mn + 20.;
        Double_t  tika_lbound      = init_mean_ti - 20.;
        Double_t  tika_ubound      = init_mean_ti + 20.;
/*
        Double_t  cakb_lbound      = init_mean_cakb - 30.;
        Double_t  cakb_ubound      = init_mean_cakb + 30.;
        Double_t  caka_lbound      = init_mean_ca - 30.;
        Double_t  caka_ubound      = init_mean_ca + 30.;
*/
        //cout << "## Initial values of CaKa1 and CaKb : " << init_mean_ca << ",  "  << init_mean_cakb << endl;

        froot = new TFile( rootfile, "READ" );

        hcal = new TH1F("hcal", "hcal", 4096, -0.5, 4095.5);

        TTree *tr_t = (TTree*)froot->Get("tr");
        Int_t ent = (Int_t)tr_t->GetEntries();
        cout << "# entries in tr: " << ent << endl;

        TBranch *adcEvent   = tr_t->GetBranch("adc");
        adcEvent->SetAddress(adc);

        for( Int_t ient=0; ient<ent; ient++ ){
            adcEvent->GetEntry(ient);
            if( adc[padc_ch] > 400 ){
                hcal->Fill((Double_t)(adc[padc_ch]));
            }
        }

        hcal->Rebin( rebin );
        cout << "# Fit histogram in :" << rootfile << endl;

        fcnFit = new TF1("fcnFit", TiMnFullFitFunc, l_bin, u_bin, npar_timn );// !!

        fcnFit->SetParNames("p0", "p1", "p2", "FANO", "CONST_N",
                            "TiKa1Ch", "TiKa1Gain", "TiKbCh", "RTiKbKa",
                            "MnKa1Ch");
        fcnFit->SetParName( 10, "MnKa1Gain ");
        fcnFit->SetParName( 11, "MnKbCh    ");
        fcnFit->SetParName( 12, "RMnKbKa   ");
        fcnFit->SetParName( 13, "EscapeGain");         // Escape/MnKa Jun 23 2010
  //      fcnFit->SetParName( 14, "FeGain    ");
  //      fcnFit->SetParName( 15, "pShift    ");
  //      fcnFit->SetParName( 16, "pGainRa   ");
        fcnFit->SetParName( 14, "tShift    ");
        fcnFit->SetParName( 15, "tGainTiKa   ");// ti tail gain formerly 15
        fcnFit->SetParName( 16, "tGainTiKb   ");
        fcnFit->SetParName( 17, "tGainMnKa   ");// mn tail gain formerly 16
        fcnFit->SetParName( 18, "tGainMnKb   ");//
        fcnFit->SetParName( 19, "tTiBeta     ");//17->19
        fcnFit->SetParName( 20, "tMnBeta     ");//17->19
    //    fcnFit->SetParName( 20, "rSigR     ");
   //     fcnFit->SetParName( 22, "CaKa1Ch"   );
   //     fcnFit->SetParName( 23, "CaKa1Gain ");
   //     fcnFit->SetParName( 24, "CaKbCh    ");
   //     fcnFit->SetParName( 25, "RCaKbKa   ");
  //      fcnFit->SetParName( 26, "tGainCa   ");
        fcnFit->SetParName( 21, "shelfGain     ");

        fcnFit->SetParameter(0, 1.);
        fcnFit->SetParLimits(0, 0.,  100.);
        fcnFit->FixParameter(1, 0.);
        fcnFit->FixParameter(2, 0.);
        fcnFit->SetParameter(3, init_fano );           // FANO
        fcnFit->SetParLimits(3, 0.09, 0.4);
        fcnFit->SetParameter(4, init_cnst );           // Constant noise 
        fcnFit->SetParLimits(4, 0, 100.);           
        fcnFit->SetParameter(5, init_mean_ti);         // TiKa1 mean in ch
        fcnFit->SetParLimits(5, tika_lbound, tika_ubound);
        fcnFit->SetParameter(6, ti_gain_setvalue);     // TiKa1 gain
        fcnFit->SetParLimits(6, 100, 1.e7);
        fcnFit->SetParameter(7, init_mean_tikb);       // TiKb mean in ch
        fcnFit->SetParLimits(7, tikb_lbound, tikb_ubound);
        fcnFit->SetParameter(8, 0.26);                  // Ti Kb to Ka gain ratio
        fcnFit->SetParLimits(8, 0.05, 0.8);
        fcnFit->SetParameter(9, init_mean_mn);         // MnKa1 mean in ch
        fcnFit->SetParLimits(9, mnka_lbound, mnka_ubound);
        fcnFit->SetParameter(10, mn_gain_setvalue);    // MnKa1 gain
        fcnFit->SetParLimits(10, 10., 1.e7);
        fcnFit->SetParameter(11, init_mean_mnkb);      // MnKb mean in ch
        fcnFit->SetParLimits(11, mnkb_lbound, mnkb_ubound);
        fcnFit->SetParameter(12, 0.7);                 // Mn Kb to Ka gain ratio -> high for LNGS data
        fcnFit->SetParLimits(12, 0.05, 0.8);
        //fcnFit->SetParameter(15, 60);                  // Pileup mean shift eV
        //fcnFit->SetParLimits(15, 10, 200);              
        //fcnFit->FixParameter(15, 150.);
        //fcnFit->FixParameter(16, 0.);
        //fcnFit->SetParameter(16, 0.01);                // Pileup gain ratio 
        //fcnFit->SetParLimits(16, 0.001, 0.10);         
        fcnFit->SetParameter(14, 50.);                 // Tail mean shift eV 
        fcnFit->SetParLimits(14, 0., 200.);           // lower limit 100 to 30. 22/02/2010
        fcnFit->SetParameter(13, 0.01);                 // Escape peak gain
        fcnFit->SetParLimits(13, 0.0001, 0.05);
        //fcnFit->FixParameter(18,  0.);

        fcnFit->SetParameter(15, 0.2);                 // Tail gain ratio ti Ka ---- set gains to 0.2 from 0.02
        fcnFit->SetParLimits(15, 0.001, 0.40);           // upper 0.2 to 0.60 
        //fcnFit->FixParameter(19,  0.);
        fcnFit->SetParameter(16, 0.2);                 // Tail gain ratio ti Kb
        fcnFit->SetParLimits(16, 0.001, 0.55);          // upper 0.2 to 0.40 22/02/2010
        fcnFit->SetParameter(17, 0.2);                 // Tail gain ratio mn Ka
        fcnFit->SetParLimits(17, 0.001, 0.40);           // upper 0.2 to 0.60 
        //fcnFit->FixParameter(19,  0.);
        fcnFit->SetParameter(18, 0.2);                 // Tail gain ratio mn Kb 
        fcnFit->SetParLimits(18, 0.001, 0.55);          // upper 0.2 to 0.40 22/02/2010

        fcnFit->SetParameter(19, 10.);                  // Ti Tail beta slope, from Sato report
        fcnFit->SetParLimits(19, 1., 100.);
        fcnFit->SetParameter(20, 10.);                  // Mn Tail beta slope, from Sato report
        fcnFit->SetParLimits(20, 1., 100.);
     //   fcnFit->SetParameter(18, 2.);                   // Sigma broadening of pileup
     //   fcnFit->SetParLimits(18, 1.5, 10);
/*
        fcnFit->SetParameter(22, init_mean_ca);         // CaKa1 mean in ch
        fcnFit->SetParLimits(22, caka_lbound, caka_ubound );         
        fcnFit->SetParameter(23, ca_gain_setvalue);     // CaKa1 gain
        fcnFit->SetParameter(24, init_mean_cakb);       // CaKb mean in ch
        fcnFit->SetParLimits(24, cakb_lbound, cakb_ubound);     
        fcnFit->SetParameter(25, 0.3);                  // Ca Kb to Ka gain ratio
        fcnFit->SetParLimits(25, 0.05, 0.6);
        fcnFit->SetParameter(26, 0.02);                 // Tail gain ratio ca 
        fcnFit->SetParLimits(26, 0.001, 0.20);          // upper 0.2 to 0.40 
*/
        fcnFit->SetParameter(21, 0.001 );               // shelf gain ratio
        fcnFit->SetParLimits(21, 0.0004, 0.05);


        //fcnFit->FixParameter(14, 0);
        //fcnFit->SetParameter(14, 0.001);               // Fe to MnKa1 gain ratio
        //fcnFit->SetParLimits(14, -0.000001, 0.5);

        fcnFit->SetNpx(range);
        cout << "# Npx :" << range << endl;
        cout << "# Set parameters for fit function." << endl;

        c1->cd(1);
        gStyle->SetLabelSize(0.06, "x");
        gStyle->SetLabelSize(0.06, "y");
        gPad->SetTicks();
        gPad->SetGridy();
        gPad->SetGridx();
        gPad->SetLogy();
        hcal->GetXaxis()->SetRangeUser( ll, ul);
        c1->Update();

        status = hcal->Fit("fcnFit", "RIE"); // here the magic (fit) happens -------- 
    //    cout << "Fit done!!!!" << endl << endl;
        fcnFit->GetParameters(pfit);
        chi2 = fcnFit->GetChisquare();
        ndf  = fcnFit->GetNDF();

        fcnResult = new TF1("fcnResult", TiMnFullFitFunc, l_bin, u_bin, npar_timn );
        fcnResult->SetLineColor(2);
        fcnResult->SetLineWidth(2);
        fcnResult->SetParameters(pfit);
        fcnResult->SetNpx(range);
        fcnResult->Draw("same");

        fcnBck   = new TF1("fcnBck",   backFunc,  l_bin,  u_bin, 3);
        fcnShelfMn1   = new TF1("fcnShelfMn1",   shelfFunc,  l_bin,  u_bin, 7);
        fcnShelfTi1   = new TF1("fcnShelfTi1",   shelfFunc,  l_bin,  u_bin, 7);
        //fcnCaka1 = new TF1("fcnCaka1", caka1Func, l_bin,  u_bin, 5);
        fcnTika1 = new TF1("fcnTika1", tika1Func, l_bin,  u_bin, 5);
        fcnMnka1 = new TF1("fcnMnka1", mnka1Func, l_bin,  u_bin, 5);

        fcnTTiK1 = new TF1("fcnTTiK1", titailFunc,  l_bin,  u_bin, 8);
        fcnTMnK1 = new TF1("fcnTMnK1", mntailFunc,  l_bin,  u_bin, 8);

        fcnTTiK2 = new TF1("fcnTTiK2", tibetatailFunc,  l_bin,  u_bin, 9);
        fcnTMnK2 = new TF1("fcnTMnK2", mnbetatailFunc,  l_bin,  u_bin, 9);

        //fcnPTiK1 = new TF1("fcnPTiK1", tipileFunc,  l_bin,  u_bin, 8);
        //fcnPMnK1 = new TF1("fcnPMnK1", mnpileFunc,  l_bin,  u_bin, 8);

        fcnEscTika1 = new TF1("fcnEscTika1", tiescFunc, l_bin, u_bin, 6);
        fcnEscMnka1 = new TF1("fcnEscMnka1", mnescFunc, l_bin, u_bin, 6);

        // Draw the annotations for the fitted lines
        Double_t  tmn = pfit[9] + 40; 
        Double_t  tti = pfit[5] + 40; 
        //Double_t  tca = pfit[22] + 40; 
        Double_t  ty_mn = pfit[10] / 10.;
        Double_t  ty_ti = pfit[6] / 10.;
 
        TLatex *Tmn = new TLatex( tmn, ty_mn, "MnKa1");
        TLatex *Tti = new TLatex( tti, ty_ti, "TiKa1");
     //   TLatex *Tca = new TLatex( tca, ty, "CaKa1");
        Tmn->SetTextSize(0.04); Tmn->SetTextAngle(90);
        Tti->SetTextSize(0.04); Tti->SetTextAngle(90);
     //   Tca->SetTextSize(0.04); Tca->SetTextAngle(90);
        Tmn->Draw("same");
        Tti->Draw("same");
     //   Tca->Draw("same");
        
        //gStyle->SetLabelSize(0.06, "x");
        //gStyle->SetLabelSize(0.06, "y");

        hcal->GetXaxis()->SetTitle("ADC channel");
        hcal->SetStats(0);


        fcnBck->SetLineColor(1);
        fcnBck->SetLineStyle(4);
        fcnBck->SetLineWidth(2);
        fcnShelfMn1->SetLineColor(4);                            // blue for shelf
        fcnShelfMn1->SetLineStyle(4);
        fcnShelfMn1->SetLineWidth(2);
        fcnShelfTi1->SetLineColor(4);                            // blue for shelf
        fcnShelfTi1->SetLineStyle(4);
        fcnShelfTi1->SetLineWidth(2);
        fcnTika1->SetLineColor(6);                                // Magenda line for Ti
        fcnTika1->SetLineWidth(2);               
        fcnTika1->SetNpx(10000);
        fcnMnka1->SetLineColor(3);                                // Green for Mn
        fcnMnka1->SetLineWidth(2);               
        fcnMnka1->SetNpx(10000);
       // fcnCaka1->SetLineColor(5);                                // Yellow line for Ca
       // fcnCaka1->SetLineWidth(2);                                
       // fcnCaka1->SetNpx(10000);
        fcnTTiK1->SetLineColor(7);                                // Cyan for tail 
        fcnTTiK1->SetLineWidth(2);
        fcnTTiK1->SetLineStyle(5);
        fcnTTiK1->SetNpx(10000);
        fcnTTiK2->SetLineColor(7);                                // Cyan for tail 
        fcnTTiK2->SetLineWidth(2);
        fcnTTiK2->SetLineStyle(5);
        fcnTTiK2->SetNpx(10000);
        fcnTMnK1->SetLineColor(7);
        fcnTMnK1->SetLineWidth(2);
        fcnTMnK1->SetLineStyle(5);
        fcnTMnK1->SetNpx(10000);
        fcnTMnK2->SetLineColor(7);
        fcnTMnK2->SetLineWidth(2);
        fcnTMnK2->SetLineStyle(5);
        fcnTMnK2->SetNpx(10000);
       // fcnPTiK1->SetLineColor(4);                                // Blue for pile-up
       // fcnPTiK1->SetLineStyle(5);
       // fcnPTiK1->SetNpx(10000);
       // fcnPMnK1->SetLineColor(4);
       // fcnPMnK1->SetLineStyle(5);
       // fcnPMnK1->SetNpx(10000);
        fcnEscTika1->SetLineColor(6);
        fcnEscTika1->SetLineWidth(2);
        fcnEscTika1->SetLineStyle(5);
        fcnEscTika1->SetNpx(10000);
        fcnEscMnka1->SetLineColor(3);
        fcnEscMnka1->SetLineWidth(2);
        fcnEscMnka1->SetLineStyle(5);
        fcnEscMnka1->SetNpx(10000);

        fcnBck->SetParameters(pfit);
        fcnShelfMn1->SetParameters(pfit[5], pfit[9], pfit[3], pfit[4], pfit[21], pfit[9], pfit[10]); // Mn Ka1 shelf function
        fcnShelfTi1->SetParameters(pfit[5], pfit[9], pfit[3], pfit[4], pfit[21], pfit[5], pfit[6]); // Ti Ka1 shelf function
        fcnTika1->SetParameters(pfit[5], pfit[9], pfit[6], pfit[3], pfit[4]);
        fcnMnka1->SetParameters(pfit[9], pfit[5], pfit[10],pfit[3], pfit[4]);
        //fcnCaka1->SetParameters(pfit[9], pfit[5], pfit[23], pfit[3], pfit[4]);
        fcnTTiK1->SetParameters(pfit[9], pfit[5], pfit[14], pfit[3], pfit[4], pfit[6], pfit[19], pfit[15]);
        fcnTMnK1->SetParameters(pfit[9], pfit[5], pfit[14], pfit[3], pfit[4], pfit[10],pfit[20], pfit[17]);

        fcnTTiK2->SetParameters(pfit[11], pfit[7], pfit[14], pfit[3], pfit[4], pfit[6], pfit[19], pfit[16], pfit[8]); // Ti Kb tail function; 
        fcnTMnK2->SetParameters(pfit[11], pfit[7], pfit[14], pfit[3], pfit[4], pfit[10],pfit[20], pfit[18], pfit[12]); // Mn Kb tail function


     //   fcnPTiK1->SetParameters(pfit[9], pfit[5], pfit[15], pfit[3], pfit[4], pfit[6], pfit[21], pfit[16]);
     //   fcnPMnK1->SetParameters(pfit[9], pfit[5], pfit[15], pfit[3], pfit[4], pfit[10],pfit[21], pfit[16]);
        fcnEscTika1->SetParameters(pfit[9], pfit[5], pfit[6], pfit[3], pfit[4], pfit[13]);
        fcnEscMnka1->SetParameters(pfit[9], pfit[5], pfit[6], pfit[3], pfit[4], pfit[13]);

        lTika1 = new TLine(pfit[5], 0, pfit[5], 1.e5);
        lTikb  = new TLine(pfit[7], 0, pfit[7], 1.e5);
        lMnka1 = new TLine(pfit[9], 0, pfit[9], 1.e5);
        lMnkb  = new TLine(pfit[11],0, pfit[11],1.e5);
       // lCaka1 = new TLine(pfit[22],0, pfit[22],1.e5);
        lTika1->SetLineStyle(2);
        lTika1->SetLineColor(6);
        lTika1->SetLineWidth(2);
        lMnka1->SetLineStyle(2); 
        lMnka1->SetLineColor(3);
        lMnka1->SetLineWidth(2);
       // lCaka1->SetLineStyle(2);
       // lCaka1->SetLineColor(5);
       // lCaka1->SetLineWidth(2);

        fcnBck->Draw("same");
        fcnTika1->Draw("same");
        fcnMnka1->Draw("same");
        fcnShelfMn1->Draw("same");
        fcnShelfTi1->Draw("same");
        fcnTTiK1->Draw("same");
        fcnTMnK1->Draw("same");
        fcnTTiK2->Draw("same");
        fcnTMnK2->Draw("same");
       // fcnPTiK1->Draw("same");
       // fcnPMnK1->Draw("same");
        fcnEscTika1->Draw("same");
        fcnEscMnka1->Draw("same");
        //lCaka1->Draw("same");
        lTika1->Draw("same");
        lMnka1->Draw("same");

        c1->Update();
        cout << "-----------------------------------" << endl;
        cout << "  chi2 = " << chi2 << ";  NDF = " << ndf << ";  chi2/ndf = " 
             << chi2 / ndf << endl;
        cout << "# Plotted fit functions to pads." << endl;
        
 /*       //////// Tail ana part //////////
        lbtl[0] = pfit[5] - 300.;
        ubtl[0] = pfit[5] + 100.;
        lbtl[1] = pfit[9] - 300.;
        ubtl[1] = pfit[9] + 100.;
        lbpk[0] = pfit[5] - 200.;
        ubpk[0] = pfit[5] + 200.;
        lbpk[1] = pfit[9] - 200.;
        ubpk[1] = pfit[9] + 200.;
        ntiT   = fcnTTiK1->Integral( lbtl[0], ubtl[0] );
        ntiTEr = sqrt( ntiT );
        //ntiTEr = fcnTTiK1->IntegralError( lbtl[0], ubtl[0] );
        ntiK   = fcnTika1->Integral( lbpk[0], ubpk[0] );
        //ntiKEr = fcnTika1->IntegralError( lbpk[0], ubpk[0] );
        ntiKEr = sqrt( ntiK );
        nmnT   = fcnTMnK1->Integral( lbtl[1], ubtl[1] );
        //nmnTEr = fcnTMnK1->IntegralError( lbtl[1], ubtl[1] );
        nmnTEr = sqrt( nmnT );
        nmnK   = fcnMnka1->Integral( lbpk[1], ubpk[1] );
        //nmnKEr = fcnMnka1->IntegralError( lbpk[1], ubpk[1] );
        nmnKEr = sqrt( nmnK );
        ntiK   = ntiK * 1.5;
        nmnK   = nmnK * 1.5;
        ntiKEr = ntiKEr * 1.5;
        nmnKEr = nmnKEr * 1.5;
        if( ntiK > 0 && nmnK > 0 ){
            rtiT = ntiT / ntiK;
            rmnT = nmnT / nmnK; 
            rtiTEr = GetError( ntiT,  ntiTEr, ntiTEr, ntiK,  ntiKEr,  ntiKEr, "DIV");
            rmnTEr = GetError( nmnT,  nmnTEr, nmnTEr, nmnK,  nmnKEr,  nmnKEr, "DIV");
        }else{
            rtiT = 0.;
            rmnT = 0.;
            rtiTEr = 0.;
            rmnTEr = 0.;
        }
        cout << "############ Number of ti tail event : " << ntiT << " +- " << ntiTEr << endl;
        cout << "#            Number of ti peak event : " << ntiK << " +- " << ntiKEr << endl;
        cout << "#            Number of mn tail event : " << nmnT << " +- " << nmnTEr << endl;
        cout << "#            Number of mn peak event : " << nmnK << " +- " << nmnKEr << endl;
        cout << "#            tail to peak ratio ti   : " << rtiT << " +- " << rtiTEr << endl;
        cout << "#            tail to peak ratio mn   : " << rmnT << " +- " << rmnTEr << endl;*/
 
        
         
    }

	

/*  --------------here I tried to write fit results to file
  ofstream myfile;
  myfile.open ("FitParameters.txt",ios::app);
  TString output = "From file: " + rootfile;
  myfile << output << endl;
  myfile << "Ti-position: "<< pXtika1 << "Mn-position: "<< pXmnka1 << "FWHM @ Mn Ka: " << fwhm;
  myfile.close();
*/

    return status;
    
}


Int_t SDDclass::FitTiCu( Bool_t saveflag )
{
    cout << "# In SDDclass::FitTiCu() ... " << endl;
    Bool_t   kmnfe;
    status = 0;          // Fit status:Return value of Fit. 0: converged; 4: abnormal termination
    chi2   = 0.;
    ndf    = 0;
    rebin  = 10;
    /// for tail ana /////
    ntiT   = 0.;
    ntiTEr = 0.;
    ntiK   = 0.;
    ntiKEr = 0.;
    ncuT   = 0.;
    ncuK   = 0.;
    ncuTEr = 0.;
    ncuKEr = 0.;
    rtiT   = 0.;
    rtiTEr = 0.;
    rcuT   = 0.;
    rcuTEr = 0.;
    for( Int_t i = 0; i <= 1; i++ ){
    lbtl[i] =  0.;
    ubtl[i] =  0.;
    lbtl[i] =  0.;
    ubtl[i] =  0.;
    lbpk[i] =  0.;
    ubpk[i] =  0.;
    lbpk[i] =  0.;
    ubpk[i] =  0.;
    }
    ////////////////////

    kmnfe = kTRUE;
    if( k2Peaks == true ){
        Double_t  l_bin = peakX[0] - 300.;
        Double_t  u_bin = peakX[1] + 600.;
        Double_t  ti_gain_setvalue = peakY[0] / 2. * rebin;
        Double_t  cu_gain_setvalue = peakY[1] / 2. * rebin;
        Double_t  init_mean_ti     = peakX[0] ;
        Double_t  init_mean_cu     = peakX[1] ;
        Double_t  init_fano        = 0.15;
        Double_t  init_cnst        = 33.;
        //Double_t  init_tikb2ka     = 0.2;
        Int_t     range = (Int_t)u_bin - (Int_t)l_bin + 1;
        range     = ( Int_t )( range / rebin );
        Double_t  init_mean_tikb   = init_mean_ti + 120.;
        Double_t  tikb_lbound      = init_mean_tikb - 80.;
        Double_t  tikb_ubound      = init_mean_tikb + 80.;
        Double_t  init_mean_cukb   = init_mean_cu + 250.;
        Double_t  cukb_lbound      = init_mean_cukb - 150.;
        Double_t  cukb_ubound      = init_mean_cukb + 200.;

        froot = new TFile( rootfile, "READ" );

        TTree *tr_t = (TTree*)froot->Get("tr");
        Int_t ent = (Int_t)tr_t->GetEntries();
        cout << "# entries in tr: " << ent << endl;
        TBranch *adcEvent   = tr_t->GetBranch("padc");
        adcEvent->SetAddress(adc);

        for( Int_t ient=0; ient<ent; ient ++ )
        {
            adcEvent->GetEntry(ient);
            if( adc[padc_ch] > 200 ){
                hcal->Fill((Double_t)(adc[padc_ch]));
            }
        }

        hcal->Rebin( rebin );
        cout << "# Fit histogram in :" << rootfile << endl;
        //fcnFit = new TF1("fcnFit", TiCuSourceMnFeFitFunc, l_bin, u_bin, nPar );
        fcnFit = new TF1("fcnFit", TiCuFullFitFunc, l_bin, u_bin, nPar );
        fcnFit->SetParNames("p0", "p1", "p2", "FANO", "CONST_N",
                            "TiKa1Ch", "TiKa1Gain", "TiKbCh", "RTiKbKa",
                            "CuKa1Ch");
        fcnFit->SetParName( 10, "CuKa1Gain ");
        fcnFit->SetParName( 11, "CuKbCh    ");
        fcnFit->SetParName( 12, "RCuKbKa   ");
        fcnFit->SetParName( 13, "MnGain    ");
        fcnFit->SetParName( 14, "FeGain    ");
        fcnFit->SetParName( 15, "pShift    ");
        fcnFit->SetParName( 16, "pGainRa   ");
        fcnFit->SetParName( 17, "tShift    ");
        fcnFit->SetParName( 18, "tGainTi   ");
        fcnFit->SetParName( 19, "tGainCu   ");
        fcnFit->SetParName( 20, "tBeta     ");
        fcnFit->SetParName( 21, "rSigR     ");

        Double_t  u = 2.; 
        Double_t  l = 0.2;
        fcnFit->SetParameter(0, bpar[0]);
        fcnFit->SetParLimits(0, bpar[0] * l,  bpar[0] * u);
        fcnFit->SetParameter(1, bpar[1]);
        fcnFit->SetParLimits(1, bpar[1] * l,  bpar[1] * u);
        fcnFit->SetParameter(2, bpar[2]);
        fcnFit->SetParLimits(2, bpar[2] * l,  bpar[2] * u);
        fcnFit->SetParameter(3, init_fano );           // FANO
        fcnFit->SetParLimits(3, 0.09, 0.4);
        fcnFit->SetParameter(4, init_cnst );           // Constant noise 
        fcnFit->SetParLimits(4, 0, 100.);           
        fcnFit->SetParameter(5, init_mean_ti);         // TiKa1 mean in ch
        fcnFit->SetParLimits(5, l_bin, 1700.);
        fcnFit->SetParameter(6, ti_gain_setvalue);     // TiKa1 gain
        fcnFit->SetParLimits(6, 100, 1.e7);
        fcnFit->SetParameter(7, init_mean_tikb);       // TiKb mean in ch
        fcnFit->SetParLimits(7, tikb_lbound, tikb_ubound);
        fcnFit->SetParameter(8, 0.3);                  // Ti Kb to Ka gain ratio
        fcnFit->SetParLimits(8, 0.05, 0.8);
        fcnFit->SetParameter(9, init_mean_cu);         // CuKa1 mean in ch
        fcnFit->SetParLimits(9, 1500., u_bin);
        fcnFit->SetParameter(10, cu_gain_setvalue);    // CuKa1 gain
        fcnFit->SetParLimits(10, 10., 1.e7);
        fcnFit->SetParameter(11, init_mean_cukb);      // CuKb mean in ch
        fcnFit->SetParLimits(11, cukb_lbound, cukb_ubound);
        fcnFit->SetParameter(12, 0.3);                 // Cu Kb to Ka gain ratio
        fcnFit->SetParLimits(12, 0.05, 0.7);
        fcnFit->SetParameter(15, 200);                 // Pileup mean shift eV
        fcnFit->SetParLimits(15, 100, 300);              
        //fcnFit->FixParameter(15, 150.);
        //fcnFit->FixParameter(16, 0.);
        fcnFit->SetParameter(16, 0.01);                // Pileup gain ratio 
        fcnFit->SetParLimits(16, 0.001, 0.30);         // upper 0.08 to 0.30 22/02/2010
        fcnFit->SetParameter(17, 150.);                // Tail mean shift eV 
        fcnFit->SetParLimits(17, 30., 250.);           // lower limit 100 to 30. 22/02/2010
        //fcnFit->FixParameter(17, 40.);
        //fcnFit->FixParameter(18,  0.);
        fcnFit->SetParameter(18, 0.04);                 // Tail gain ratio ti
        fcnFit->SetParLimits(18, 0.01, 0.40);           // upper 0.2 to 0.40 
        //fcnFit->FixParameter(19,  0.);
        fcnFit->SetParameter(19, 0.04);                 // Tail gain ratio cu 
        fcnFit->SetParLimits(19, 0.01, 0.60);           // upper 0.2 to 0.40 22/02/2010
        fcnFit->SetParameter(20, 10.);                  // Tail beta slope, from Sato report
        fcnFit->SetParLimits(20, 1., 100.);
        fcnFit->SetParameter(21, 2.);                   // Sigma broadening of pileup
        //fcnFit->SetParLimits(21, 1., 3);
        fcnFit->FixParameter(21, 2.);

        if( kmnfe ){
            fcnFit->SetParameter(13, 1);               // Mn to CuKa1 gain ratio
            fcnFit->SetParLimits(13, -0.000001, 2);
            //fcnFit->FixParameter(13, 0);
            fcnFit->SetParameter(14, 0.001);               // Fe to CuKa1 gain ratio
            fcnFit->SetParLimits(14, -0.000001, 0.1);
        }else{
            fcnFit->FixParameter(13, 0.);
            fcnFit->FixParameter(14, 0.);
        }
        fcnFit->SetNpx(range);
        cout << "# Npx :" << range << endl;
        cout << "# Set parameters for fit function." << endl;

        c1->cd(1);
        gPad->SetTicks();
        gPad->SetGridy();
        gPad->SetGridx();
        gPad->SetLogy();
        hcal->GetXaxis()->SetRangeUser( ll, ul );
        status = hcal->Fit("fcnFit", "R");
        fcnFit->GetParameters(pfit);
        chi2 = fcnFit->GetChisquare();
        ndf  = fcnFit->GetNDF();

        //fcnResult = new TF1("fcnResult", TiCuSourceMnFeFitFunc, l_bin, u_bin, nPar );
        fcnResult = new TF1("fcnResult", TiCuFullFitFunc, l_bin, u_bin, nPar );
        fcnResult->SetLineColor(2);
        fcnResult->SetParameters(pfit);
        fcnResult->Draw("same");

        Double_t  fem = (pfit[9]+pfit[5])/2 + (FeKa1-(TiKa1+CuKa1)/2) * (pfit[9]-pfit[5])/(CuKa1-TiKa1); 
        fcnBck   = new TF1("fcnBck",   backFunc,  l_bin,  u_bin, 3);
        fcnTika1 = new TF1("fcnTika1", tika1Func, l_bin,  u_bin, 5);
        fcnCuka1 = new TF1("fcnCuka1", cuka1Func, l_bin,  u_bin, 5);
        fcnTTiK1 = new TF1("fcnTTiK1", titailFunc,  l_bin,  u_bin, 8);
        fcnTCuK1 = new TF1("fcnTCuK1", cutailFunc,  l_bin,  u_bin, 8);
        fcnPTiK1 = new TF1("fcnPTiK1", tipileFunc,  l_bin,  u_bin, 8);
        fcnPCuK1 = new TF1("fcnPCuK1", cupileFunc,  l_bin,  u_bin, 8);
        fcnBck->SetLineColor(1);
        fcnBck->SetLineWidth(1);
        fcnTika1->SetLineColor(2);                                // Red line for Ti
        fcnTika1->SetNpx(10000);
        fcnCuka1->SetLineColor(3);                                // Green for Cu
        fcnCuka1->SetNpx(10000);
        fcnTTiK1->SetLineColor(7);
        fcnTTiK1->SetLineStyle(5);
        fcnTTiK1->SetNpx(10000);
        fcnTCuK1->SetLineColor(7);
        fcnTCuK1->SetLineStyle(5);
        fcnTCuK1->SetNpx(10000);
        fcnPTiK1->SetLineColor(4);
        fcnPTiK1->SetLineStyle(5);
        fcnPTiK1->SetNpx(10000);
        fcnPCuK1->SetLineColor(4);
        fcnPCuK1->SetLineStyle(5);
        fcnPCuK1->SetNpx(10000);
        fcnBck->SetParameters(pfit);
        fcnTika1->SetParameters(pfit[5], pfit[9], pfit[6], pfit[3], pfit[4]);
        fcnCuka1->SetParameters(pfit[9], pfit[5], pfit[10],pfit[3], pfit[4]);
        fcnTTiK1->SetParameters(pfit[9], pfit[5], pfit[17], pfit[3], pfit[4], pfit[6], pfit[20], pfit[18]);
        fcnTCuK1->SetParameters(pfit[9], pfit[5], pfit[17], pfit[3], pfit[4], pfit[10],pfit[20], pfit[19]);
        fcnPTiK1->SetParameters(pfit[9], pfit[5], pfit[15], pfit[3], pfit[4], pfit[6], pfit[21], pfit[16]);
        fcnPCuK1->SetParameters(pfit[9], pfit[5], pfit[15], pfit[3], pfit[4], pfit[10],pfit[21], pfit[16]);
        lTika1 = new TLine(pfit[5], 0, pfit[5], 1.e5);
        lTikb  = new TLine(pfit[7], 0, pfit[7], 1.e5);
        lCuka1 = new TLine(pfit[9], 0, pfit[9], 1.e5);
        lCukb  = new TLine(pfit[11],0, pfit[11],1.e5);
        lFe    = new TLine( fem    ,0, fem,     1.e5);
        lTika1->SetLineStyle(2);
        lTika1->SetLineColor(2);
        lCuka1->SetLineStyle(2); 
        lCuka1->SetLineColor(3);
        lFe   ->SetLineStyle(3);
        lFe   ->SetLineColor(6);
        fcnBck->Draw("same");
        fcnTika1->Draw("same");
        fcnCuka1->Draw("same");
        fcnTTiK1->Draw("same");
        fcnTCuK1->Draw("same");
        fcnPTiK1->Draw("same");
        fcnPCuK1->Draw("same");
        lTika1->Draw("same");
        lCuka1->Draw("same");
        lFe   ->Draw("same");
        c1->Update();
        cout << "-----------------------------------" << endl;
        cout << "  chi2 = " << chi2 << ";  NDF = " << ndf << ";  chi2/ndf = " 
             << chi2 / ndf << endl;
        cout << "# Plotted fit functions to pads." << endl;
        
        //////// Tail ana part //////////
        lbtl[0] = pfit[5] - 300.;
        ubtl[0] = pfit[5] + 100.;
        lbtl[1] = pfit[9] - 300.;
        ubtl[1] = pfit[9] + 100.;
        lbpk[0] = pfit[5] - 200.;
        ubpk[0] = pfit[5] + 200.;
        lbpk[1] = pfit[9] - 200.;
        ubpk[1] = pfit[9] + 200.;
        ntiT   = fcnTTiK1->Integral( lbtl[0], ubtl[0] );
        ntiTEr = sqrt( ntiT );
        //ntiTEr = fcnTTiK1->IntegralError( lbtl[0], ubtl[0] );
        ntiK   = fcnTika1->Integral( lbpk[0], ubpk[0] );
        //ntiKEr = fcnTika1->IntegralError( lbpk[0], ubpk[0] );
        ntiKEr = sqrt( ntiK );
        ncuT   = fcnTCuK1->Integral( lbtl[1], ubtl[1] );
        //ncuTEr = fcnTCuK1->IntegralError( lbtl[1], ubtl[1] );
        ncuTEr = sqrt( ncuT );
        ncuK   = fcnCuka1->Integral( lbpk[1], ubpk[1] );
        //ncuKEr = fcnCuka1->IntegralError( lbpk[1], ubpk[1] );
        ncuKEr = sqrt( ncuK );
        ntiK   = ntiK * 1.5;
        ncuK   = ncuK * 1.5;
        ntiKEr = ntiKEr * 1.5;
        ncuKEr = ncuKEr * 1.5;
        if( ntiK > 0 && ncuK > 0 ){
            rtiT = ntiT / ntiK;
            rcuT = ncuT / ncuK; 
            rtiTEr = GetError( ntiT,  ntiTEr, ntiTEr, ntiK,  ntiKEr,  ntiKEr, "DIV");
            rmnTEr = GetError( nmnT,  nmnTEr, nmnTEr, nmnK,  nmnKEr,  nmnKEr, "DIV");
        }else{
            rtiT = 0.;
            rcuT = 0.;
            rtiTEr = 0.;
            rcuTEr = 0.;
        }
        cout << "############ Number of ti tail event : " << ntiT << " +- " << ntiTEr << endl;
        cout << "#            Number of ti peak event : " << ntiK << " +- " << ntiKEr << endl;
        cout << "#            Number of cu tail event : " << ncuT << " +- " << ncuTEr << endl;
        cout << "#            Number of cu peak event : " << ncuK << " +- " << ncuKEr << endl;
        cout << "#            tail to peak ratio ti   : " << rtiT << " +- " << rtiTEr << endl;
        cout << "#            tail to peak ratio cu   : " << rcuT << " +- " << rcuTEr << endl;
    }

    return status;
}


Int_t SDDclass::FitTiMnCu( Bool_t saveflag )
{
    cout << "# In SDDclass::FitTiMnCu() ... " << endl;
    Bool_t   kmnfe;
    status = 0;          // Fit status:Return value of Fit. 0: converged; 4: abnormal termination
    chi2   = 0.;
    ndf    = 0;
    rebin  = 2;
    /// for tail ana /////
    ntiT   = 0.;
    ntiTEr = 0.;
    ntiK   = 0.;
    ntiKEr = 0.;
    ncuT   = 0.;
    ncuK   = 0.;
    ncuTEr = 0.;
    ncuKEr = 0.;
    rtiT   = 0.;
    rtiTEr = 0.;
    rcuT   = 0.;
    rcuTEr = 0.;
    for( Int_t i = 0; i <= 1; i++ ){
    lbtl[i] =  0.;
    ubtl[i] =  0.;
    lbtl[i] =  0.;
    ubtl[i] =  0.;
    lbpk[i] =  0.;
    ubpk[i] =  0.;
    lbpk[i] =  0.;
    ubpk[i] =  0.;
    }
    ////////////////////

    kmnfe = kTRUE;
    if( k2Peaks == true ){
        Double_t  ti_gain_setvalue = peakY[0] * rebin;
        Double_t  mn_gain_setvalue = peakY[1] * rebin;
        Double_t  cu_gain_setvalue = mn_gain_setvalue/30;
        Double_t  bg_setvalue      = mn_gain_setvalue/1000;
        Double_t  init_mean_ti     = peakX[0] ;
        Double_t  init_mean_mn     = peakX[1] ;
        Double_t  init_fano        = 0.15;
        Double_t  init_cnst        = 10.;
        Double_t  slope            = ( init_mean_mn - init_mean_ti )/( MnKa1 - TiKa1 ); // in ch/eV
        Double_t  init_mean_cu     = init_mean_mn + ( CuKa1 - MnKa1 ) *  slope;
        Double_t  l_bin = peakX[0] - 300.;
        Double_t  u_bin = init_mean_cu + 300.; // put further up to be able to calculate the background better

        //cout << "l-bin = " << l_bin << "; u_bin = " << u_bin << " " << init_mean_cu << endl;

        Int_t     range = (Int_t)u_bin - (Int_t)l_bin + 1;
        //range     = ( Int_t )( range / rebin );
        Double_t  init_mean_tikb   = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-TiKb1)*slope;
        Double_t  tikb_lbound      = init_mean_tikb - 20.;
        Double_t  tikb_ubound      = init_mean_tikb + 20.;
        Double_t  init_mean_mnkb   = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-MnKb1)*slope;
        Double_t  mnkb_lbound      = init_mean_mnkb - 20.;
        Double_t  mnkb_ubound      = init_mean_mnkb + 20.;
        Double_t  init_mean_cukb   = (init_mean_ti+init_mean_mn)/2 - ((MnKa1+TiKa1)/2-CuKb1)*slope;
        Double_t  cukb_lbound      = init_mean_cukb - 20.;
        Double_t  cukb_ubound      = init_mean_cukb + 20.;

        Double_t  mnka_lbound      = init_mean_mn - 20.;
        Double_t  mnka_ubound      = init_mean_mn + 20.;
        Double_t  tika_lbound      = init_mean_ti - 20.;
        Double_t  tika_ubound      = init_mean_ti + 20.;
        Double_t  cuka_lbound      = init_mean_cu - 20.;
        Double_t  cuka_ubound      = init_mean_cu + 20.;

        froot = new TFile( rootfile, "READ" );
        hcal = new TH1F("hcal", "hcal", 4096, -0.5, 4095.5);

        TTree *tr_t = (TTree*)froot->Get("tr");     
        Int_t ent = (Int_t)tr_t->GetEntries();
        cout << "# entries in tr: " << ent << endl;

        TBranch *adcEvent   = tr_t->GetBranch("adc");
        adcEvent->SetAddress(adc);

        for( Int_t ient=0; ient<ent; ient++ ){
            adcEvent->GetEntry(ient);
            if( adc[padc_ch] > 400 ){
                hcal->Fill((Double_t)(adc[padc_ch]));
            }
        }

        hcal->Rebin( rebin );
        cout << "# Fit histogram in :" << rootfile << endl;
        fcnFit = new TF1("fcnFit", TiMnCuFullFitFunc, l_bin, u_bin, nPar );

        fcnFit->SetParNames("p0", "p1", "p2", "FANO", "CONST_N", // parameters 0 ... 9
                            "TiKa1Ch", "TiKa1Gain", "TiKbCh", "RTiKbKa",
                            "MnKa1Ch");
        fcnFit->SetParName( 10, "MnKa1Gain ");
        fcnFit->SetParName( 11, "MnKbCh    ");
        fcnFit->SetParName( 12, "RMnKbKa   ");
        fcnFit->SetParName( 13, "EscapeGain");        
        fcnFit->SetParName( 14, "tShift    ");
        fcnFit->SetParName( 15, "tGainTiKa   ");
        fcnFit->SetParName( 16, "tGainTiKb   ");
        fcnFit->SetParName( 17, "tGainMnKa   ");
        fcnFit->SetParName( 18, "tGainMnKb   ");
        fcnFit->SetParName( 19, "tTiBeta     ");
        fcnFit->SetParName( 20, "tMnBeta     ");
        fcnFit->SetParName( 21, "shelfGain   ");
        fcnFit->SetParName( 22, "CuKa1Ch   ");
        fcnFit->SetParName( 23, "CuKa1Gain     ");
        fcnFit->SetParName( 24, "CuKbCh   ");
        fcnFit->SetParName( 25, "RCuKbKa      ");
        fcnFit->SetParName( 26, "tGainCu     ");


        fcnFit->SetParameter(0, bg_setvalue);                 // constant background
        fcnFit->SetParLimits(0, 0.,  bg_setvalue*10);
        /*fcnFit->SetParameter(1, 50.);                   // exponential gain
        fcnFit->SetParLimits(1, 1., 500.);                   
        fcnFit->SetParameter(2, 0.0001);                //exponential slope (positive!)
        fcnFit->SetParLimits(2, 0., 0.001); */
        fcnFit->FixParameter(1 , 0.);
        fcnFit->FixParameter(2 , 0.);             
        fcnFit->SetParameter(3, init_fano );           // FANO
        fcnFit->SetParLimits(3, 0.09, 0.4);
        fcnFit->SetParameter(4, init_cnst );           // Constant noise 
        fcnFit->SetParLimits(4, 0, 100.);           
        fcnFit->SetParameter(5, init_mean_ti);         // TiKa1 mean in ch
        fcnFit->SetParLimits(5, tika_lbound, tika_ubound);
        fcnFit->SetParameter(6, ti_gain_setvalue);     // TiKa1 gain
        fcnFit->SetParLimits(6, 100, 1.e7);
        fcnFit->SetParameter(7, init_mean_tikb);       // TiKb mean in ch
        fcnFit->SetParLimits(7, tikb_lbound, tikb_ubound);
        fcnFit->SetParameter(8, 0.26);                  // Ti Kb to Ka gain ratio
        fcnFit->SetParLimits(8, 0.05, 0.8);
        fcnFit->SetParameter(9, init_mean_mn);         // MnKa1 mean in ch
        fcnFit->SetParLimits(9, mnka_lbound, mnka_ubound);
        fcnFit->SetParameter(10, mn_gain_setvalue);    // MnKa1 gain
        fcnFit->SetParLimits(10, 10., 1.e7);
        fcnFit->SetParameter(11, init_mean_mnkb);      // MnKb mean in ch
        fcnFit->SetParLimits(11, mnkb_lbound, mnkb_ubound);
        fcnFit->SetParameter(12, 0.7);                 // Mn Kb to Ka gain ratio -> high for LNGS data
        fcnFit->SetParLimits(12, 0.05, 0.8);        
        fcnFit->SetParameter(14, 50.);                 // Tail mean shift eV !!
        fcnFit->SetParLimits(14, 0., 200.);           // lower limit 100 to 30. 22/02/2010
        fcnFit->SetParameter(13, 0.01);                 // Escape peak gain
        fcnFit->SetParLimits(13, 0.0001, 0.05);


        fcnFit->SetParameter(15, 0.2);                 // Tail gain ratio ti Ka ---- set gains to 0.2 from 0.02
        fcnFit->SetParLimits(15, 0.001, 0.40);           // upper 0.2 to 0.60 
        //fcnFit->FixParameter(19,  0.);
        fcnFit->SetParameter(16, 0.2);                 // Tail gain ratio ti Kb
        fcnFit->SetParLimits(16, 0.001, 0.55);          // upper 0.2 to 0.40 22/02/2010
        fcnFit->SetParameter(17, 0.2);                 // Tail gain ratio mn Ka
        fcnFit->SetParLimits(17, 0.001, 0.40);           // upper 0.2 to 0.60 

        fcnFit->SetParameter(18, 0.2);                 // Tail gain ratio mn Kb 
        fcnFit->SetParLimits(18, 0.001, 0.55);          // upper 0.2 to 0.40 22/02/2010

        fcnFit->SetParameter(19, 10.);                  // Ti Tail beta slope, from Sato report
        fcnFit->SetParLimits(19, 1., 100.);
        fcnFit->SetParameter(20, 10.);                  // Mn Tail beta slope, from Sato report
        fcnFit->SetParLimits(20, 1., 100.);
        fcnFit->SetParameter(21, 0.001);                  // shelf gain
        //fcnFit->SetParLimits(21, 1., 100.);

        fcnFit->SetParameter(22, init_mean_cu);         // CuKa1 mean in ch
        fcnFit->SetParLimits(22, cuka_lbound, cuka_ubound);
        fcnFit->SetParameter(23, cu_gain_setvalue);     // CuKa1 gain
        fcnFit->SetParLimits(23, 1, 1.e6);
        fcnFit->SetParameter(24, init_mean_cukb);       // CuKb mean in ch
        fcnFit->SetParLimits(24, cukb_lbound, cukb_ubound);
        fcnFit->SetParameter(25, 0.26);                  // Cu Kb to Ka gain ratio
        fcnFit->SetParLimits(25, 0.05, 0.8);
        fcnFit->SetParameter(26, 0.2);                 // Tail gain ratio all Cu lines
        fcnFit->SetParLimits(26, 0.001, 0.40);           // upper 0.2 to 0.60 
       // Cu tail slope is taken from mn tail slope for testing
       // cu also has no escape peaks for now and no shelf




        fcnFit->SetNpx(range);
        cout << "# Npx :" << range << endl;
        cout << "# Set parameters for fit function." << endl;

        c1->cd(1);
        gStyle->SetLabelSize(0.06, "x");
        gStyle->SetLabelSize(0.06, "y");
        gPad->SetTicks();
        gPad->SetGridy();
        gPad->SetGridx();
        gPad->SetLogy();
        hcal->GetXaxis()->SetRangeUser( ll, ul);
        c1->Update();

        status = hcal->Fit("fcnFit", "R"); // here the magic (fit) happens --------
	// "R" Use the Range specified in the function range, "E" Perform better Errors estimation using Minos technique, "I" Use integral of function in bin, normalized by the bin volume, instead of value at bin center
        fcnFit->GetParameters(pfit);
        chi2 = fcnFit->GetChisquare();
        ndf  = fcnFit->GetNDF();

        fcnResult = new TF1("fcnResult", TiMnCuFullFitFunc, l_bin, u_bin, nPar );
        fcnResult->SetLineColor(2);
        fcnResult->SetLineWidth(2);
        fcnResult->SetParameters(pfit);
        fcnResult->SetNpx(range);
        fcnResult->Draw("same");

        fcnBck   = new TF1("fcnBck",   backFunc,  l_bin,  u_bin, 3);
        fcnShelfMn1   = new TF1("fcnShelfMn1",   shelfFunc,  l_bin,  u_bin, 7);
        fcnShelfTi1   = new TF1("fcnShelfTi1",   shelfFunc,  l_bin,  u_bin, 7);

        fcnTika1 = new TF1("fcnTika1", tika1Func, l_bin,  u_bin, 6);
        fcnMnka1 = new TF1("fcnMnka1", mnka1Func, l_bin,  u_bin, 5);
        fcnCuka1 = new TF1("fcnCuka1", cuka1Func, l_bin,  u_bin, 5);

        fcnTTiK1 = new TF1("fcnTTiK1", titailFunc,  l_bin,  u_bin, 8);
        fcnTMnK1 = new TF1("fcnTMnK1", mntailFunc,  l_bin,  u_bin, 8);
        fcnTCuK1 = new TF1("fcnTCuK1", cutailFunc,  l_bin,  u_bin, 9);

        fcnTTiK2 = new TF1("fcnTTiK2", tibetatailFunc,  l_bin,  u_bin, 9);
        fcnTMnK2 = new TF1("fcnTMnK2", mnbetatailFunc,  l_bin,  u_bin, 9);

        fcnEscTika1 = new TF1("fcnEscTika1", tiescFunc, l_bin, u_bin, 6);
        fcnEscMnka1 = new TF1("fcnEscMnka1", mnescFunc, l_bin, u_bin, 6);

        // Draw the annotations for the fitted lines
        Double_t  tmn = pfit[9] + 40; 
        Double_t  tti = pfit[5] + 40; 
        Double_t  tcu = pfit[22] + 40; 
        //Double_t  tca = pfit[22] + 40; 
        Double_t  ty = pfit[6] / 20.;
        //Double_t  ty_ti = pfit[6] / 10.;
        //Double_t  ty_cu = pfit[23] / 10.;
 
        TLatex *Tmn = new TLatex( tmn, ty, "MnKa1");
        TLatex *Tti = new TLatex( tti, ty, "TiKa1");
        TLatex *Tcu = new TLatex( tcu, ty, "CuKa1");
     //   TLatex *Tca = new TLatex( tca, ty, "CaKa1");
        Tmn->SetTextSize(0.04); Tmn->SetTextAngle(90);
        Tti->SetTextSize(0.04); Tti->SetTextAngle(90);
        Tcu->SetTextSize(0.04); Tcu->SetTextAngle(90);
     //   Tca->SetTextSize(0.04); Tca->SetTextAngle(90);
        Tmn->Draw("same");
        Tti->Draw("same");
        Tcu->Draw("same");
     //   Tca->Draw("same");

        //gStyle->SetLabelSize(0.06, "x");
        //gStyle->SetLabelSize(0.06, "y");
        
        hcal->GetXaxis()->SetTitle("ADC channel");
        hcal->SetStats(0);


        fcnBck->SetLineColor(1);
        fcnBck->SetLineStyle(4);
        fcnBck->SetLineWidth(2);
        fcnShelfMn1->SetLineColor(4);                            // blue for shelf
        fcnShelfMn1->SetLineStyle(4);
        fcnShelfMn1->SetLineWidth(2);
        fcnShelfTi1->SetLineColor(4);                            // blue for shelf
        fcnShelfTi1->SetLineStyle(4);
        fcnShelfTi1->SetLineWidth(2);

        fcnTika1->SetLineColor(6);                                // Magenda line for Ti
        fcnTika1->SetLineWidth(2);               
        fcnTika1->SetNpx(10000);
        fcnCuka1->SetLineColor(28);                                // brownish line for Cu
        fcnCuka1->SetLineWidth(2);               
        fcnCuka1->SetNpx(10000);
        fcnMnka1->SetLineColor(3);                                // Green for Mn
        fcnMnka1->SetLineWidth(2);               
        fcnMnka1->SetNpx(10000);

        fcnTTiK1->SetLineColor(7);                                // Cyan for tail 
        fcnTTiK1->SetLineWidth(2);
        fcnTTiK1->SetLineStyle(5);
        fcnTTiK1->SetNpx(10000);
        fcnTTiK2->SetLineColor(7);                                // Cyan for tail 
        fcnTTiK2->SetLineWidth(2);
        fcnTTiK2->SetLineStyle(5);
        fcnTTiK2->SetNpx(10000);
        fcnTMnK1->SetLineColor(7);
        fcnTMnK1->SetLineWidth(2);
        fcnTMnK1->SetLineStyle(5);
        fcnTMnK1->SetNpx(10000);
        fcnTMnK2->SetLineColor(7);
        fcnTMnK2->SetLineWidth(2);
        fcnTMnK2->SetLineStyle(5);
        fcnTMnK2->SetNpx(10000);
        fcnTCuK1->SetLineColor(7);
        fcnTCuK1->SetLineWidth(2);
        fcnTCuK1->SetLineStyle(5);
        fcnTCuK1->SetNpx(10000);

        fcnEscTika1->SetLineColor(6);
        fcnEscTika1->SetLineWidth(2);
        fcnEscTika1->SetLineStyle(5);
        fcnEscTika1->SetNpx(10000);
        fcnEscMnka1->SetLineColor(3);
        fcnEscMnka1->SetLineWidth(2);
        fcnEscMnka1->SetLineStyle(5);
        fcnEscMnka1->SetNpx(10000);

        fcnBck->SetParameters(pfit);
        fcnShelfMn1->SetParameters(pfit[5], pfit[9], pfit[3], pfit[4], pfit[21], pfit[9], pfit[10]); // Mn Ka1 shelf function
        fcnShelfTi1->SetParameters(pfit[5], pfit[9], pfit[3], pfit[4], pfit[21], pfit[5], pfit[6]); // Ti Ka1 shelf function

        fcnTika1->SetParameters(pfit[5], pfit[9], pfit[6], pfit[3], pfit[4], 0);
        fcnMnka1->SetParameters(pfit[9], pfit[5], pfit[10],pfit[3], pfit[4]);
        fcnCuka1->SetParameters(pfit[22], pfit[9], pfit[5], pfit[23], pfit[3], pfit[4]); // cu ka1 gaussian peak

        fcnTTiK1->SetParameters(pfit[9], pfit[5], pfit[14], pfit[3], pfit[4], pfit[6], pfit[19], pfit[15]);
        fcnTMnK1->SetParameters(pfit[9], pfit[5], pfit[14], pfit[3], pfit[4], pfit[10],pfit[20], pfit[17]);
        fcnTCuK1->SetParameters(pfit[9], pfit[5], pfit[22], pfit[14], pfit[3], pfit[4], pfit[23], pfit[20], pfit[26]); // cupper ka1 tail function

        fcnTTiK2->SetParameters(pfit[11], pfit[7], pfit[14], pfit[3], pfit[4], pfit[6], pfit[19], pfit[16], pfit[8]); // Ti Kb tail function; 
        fcnTMnK2->SetParameters(pfit[11], pfit[7], pfit[14], pfit[3], pfit[4], pfit[10],pfit[20], pfit[18], pfit[12]); // Mn Kb tail function



        fcnEscTika1->SetParameters(pfit[9], pfit[5], pfit[6], pfit[3], pfit[4], pfit[13]);
        fcnEscMnka1->SetParameters(pfit[9], pfit[5], pfit[6], pfit[3], pfit[4], pfit[13]);

        lTika1 = new TLine(pfit[5], 0, pfit[5], 1.e5);
        //lTikb  = new TLine(pfit[7], 0, pfit[7], 1.e5);
        lMnka1 = new TLine(pfit[9], 0, pfit[9], 1.e5);
        //lMnkb  = new TLine(pfit[11],0, pfit[11],1.e5);
        lCuka1 = new TLine(pfit[22], 0, pfit[22], 1.e5);

        lTika1->SetLineStyle(2);
        lTika1->SetLineColor(6);
        lTika1->SetLineWidth(2);
        lMnka1->SetLineStyle(2); 
        lMnka1->SetLineColor(3);
        lMnka1->SetLineWidth(2);
        lCuka1->SetLineStyle(2); 
        lCuka1->SetLineColor(28);
        lCuka1->SetLineWidth(2);


        fcnBck->Draw("same");
        fcnTika1->Draw("same");
        fcnMnka1->Draw("same");
        fcnCuka1->Draw("same");
        fcnShelfMn1->Draw("same");
        fcnShelfTi1->Draw("same");
        fcnTTiK1->Draw("same");
        fcnTMnK1->Draw("same");
        fcnTTiK2->Draw("same");
        fcnTMnK2->Draw("same");
        fcnTCuK1->Draw("same");
        fcnEscTika1->Draw("same");
        fcnEscMnka1->Draw("same");
        lTika1->Draw("same");
        lMnka1->Draw("same");
        lCuka1->Draw("same");

        c1->Update();
        cout << "-----------------------------------" << endl;
        cout << "  chi2 = " << chi2 << ";  NDF = " << ndf << ";  chi2/ndf = " 
             << chi2 / ndf << endl;
        cout << "# Plotted fit functions to pads." << endl;
        
 /*       //////// Tail ana part //////////
        lbtl[0] = pfit[5] - 300.;
        ubtl[0] = pfit[5] + 100.;
        lbtl[1] = pfit[9] - 300.;
        ubtl[1] = pfit[9] + 100.;
        lbpk[0] = pfit[5] - 200.;
        ubpk[0] = pfit[5] + 200.;
        lbpk[1] = pfit[9] - 200.;
        ubpk[1] = pfit[9] + 200.;
        ntiT   = fcnTTiK1->Integral( lbtl[0], ubtl[0] );
        ntiTEr = sqrt( ntiT );
        //ntiTEr = fcnTTiK1->IntegralError( lbtl[0], ubtl[0] );
        ntiK   = fcnTika1->Integral( lbpk[0], ubpk[0] );
        //ntiKEr = fcnTika1->IntegralError( lbpk[0], ubpk[0] );
        ntiKEr = sqrt( ntiK );
        ncuT   = fcnTCuK1->Integral( lbtl[1], ubtl[1] );
        //ncuTEr = fcnTCuK1->IntegralError( lbtl[1], ubtl[1] );
        ncuTEr = sqrt( ncuT );
        ncuK   = fcnCuka1->Integral( lbpk[1], ubpk[1] );
        //ncuKEr = fcnCuka1->IntegralError( lbpk[1], ubpk[1] );
        ncuKEr = sqrt( ncuK );
        ntiK   = ntiK * 1.5;
        ncuK   = ncuK * 1.5;
        ntiKEr = ntiKEr * 1.5;
        ncuKEr = ncuKEr * 1.5;
        if( ntiK > 0 && ncuK > 0 ){
            rtiT = ntiT / ntiK;
            rcuT = ncuT / ncuK; 
            rtiTEr = GetError( ntiT,  ntiTEr, ntiTEr, ntiK,  ntiKEr,  ntiKEr, "DIV");
            rmnTEr = GetError( nmnT,  nmnTEr, nmnTEr, nmnK,  nmnKEr,  nmnKEr, "DIV");
        }else{
            rtiT = 0.;
            rcuT = 0.;
            rtiTEr = 0.;
            rcuTEr = 0.;
        }
        cout << "############ Number of ti tail event : " << ntiT << " +- " << ntiTEr << endl;
        cout << "#            Number of ti peak event : " << ntiK << " +- " << ntiKEr << endl;
        cout << "#            Number of cu tail event : " << ncuT << " +- " << ncuTEr << endl;
        cout << "#            Number of cu peak event : " << ncuK << " +- " << ncuKEr << endl;
        cout << "#            tail to peak ratio ti   : " << rtiT << " +- " << rtiTEr << endl;
        cout << "#            tail to peak ratio cu   : " << rcuT << " +- " << rcuTEr << endl;
 */   }

    return status;
}


Int_t SDDclass::FitTiCuZr( Bool_t saveflag )
{
    Int_t nPar_zr = 29;
    ll = 400; // setting limits for the hcal histogram (and its residuals)
    ul = 3500;
    cout << "# In SDDclass::FitTiCuZr() ... " << endl;
    Bool_t   kmnfe;
    status = 0;          // Fit status:Return value of Fit. 0: converged; 4: abnormal termination
    chi2   = 0.;
    ndf    = 0;
    rebin  = 2;
    /// for tail ana /////
    ntiT   = 0.;
    ntiTEr = 0.;
    ntiK   = 0.;
    ntiKEr = 0.;
    ncuT   = 0.;
    ncuK   = 0.;
    ncuTEr = 0.;
    ncuKEr = 0.;
    rtiT   = 0.;
    rtiTEr = 0.;
    rcuT   = 0.;
    rcuTEr = 0.;
    for( Int_t i = 0; i <= 1; i++ ){
    lbtl[i] =  0.;
    ubtl[i] =  0.;
    lbtl[i] =  0.;
    ubtl[i] =  0.;
    lbpk[i] =  0.;
    ubpk[i] =  0.;
    lbpk[i] =  0.;
    ubpk[i] =  0.;
    }
    ////////////////////

    kmnfe = kTRUE;
    if( k2Peaks == true ){
        Double_t  ti_gain_setvalue = peakY[0] * rebin;
        Double_t  cu_gain_setvalue = peakY[1] * rebin;
        Double_t  zr_gain_setvalue = cu_gain_setvalue*10;
        Double_t  init_mean_ti     = peakX[0] ;
        Double_t  init_mean_cu     = peakX[1] ;

        Double_t  init_fano        = 0.15;
	Double_t  init_cnst	   = 0.001;
        Double_t  slope            = ( init_mean_cu - init_mean_ti )/( CuKa1 - TiKa1 ); // in ch/eV
        Double_t  init_mean_zr     = init_mean_cu + ( ZrKa1 - CuKa1 ) *  slope;
	Double_t  init_brems       = init_mean_zr + 800;
        Double_t  l_bin = peakX[0] - 300.;
        Double_t  u_bin = init_mean_zr + 1500.; // put further up to be able to calculate the background better

        //cout << "l-bin = " << l_bin << "; u_bin = " << u_bin << " " << init_mean_cu << endl;

        Int_t     range = (Int_t)u_bin - (Int_t)l_bin + 1;
        //range     = ( Int_t )( range / rebin );
        Double_t  init_mean_tikb   = (init_mean_ti+init_mean_cu)/2 - ((CuKa1+TiKa1)/2-TiKb1)*slope;
        Double_t  tikb_lbound      = init_mean_tikb - 20.;
        Double_t  tikb_ubound      = init_mean_tikb + 20.;
        Double_t  init_mean_cukb   = (init_mean_ti+init_mean_cu)/2 + (CuKb1-(CuKa1+TiKa1)/2)*slope;
        Double_t  cukb_lbound      = init_mean_cukb - 20.;
        Double_t  cukb_ubound      = init_mean_cukb + 20.;
        Double_t  init_mean_zrkb   = init_mean_zr + (ZrKb1 - ZrKa1)*slope;
        Double_t  zrkb_lbound      = init_mean_zrkb - 50.;
        Double_t  zrkb_ubound      = init_mean_zrkb + 50.;

        Double_t  zrka_lbound      = init_mean_zr - 50.;
        Double_t  zrka_ubound      = init_mean_zr + 50.;
        Double_t  tika_lbound      = init_mean_ti - 20.;
        Double_t  tika_ubound      = init_mean_ti + 20.;
        Double_t  cuka_lbound      = init_mean_cu - 20.;
        Double_t  cuka_ubound      = init_mean_cu + 20.;

        froot = new TFile( rootfile, "READ" );
        hcal = new TH1F("hcal", "hcal", 4096, -0.5, 4095.5);

        TTree *tr_t = (TTree*)froot->Get("tr");     
        Int_t ent = (Int_t)tr_t->GetEntries();
        cout << "# entries in tr: " << ent << endl;

        TBranch *adcEvent   = tr_t->GetBranch("adc");
        adcEvent->SetAddress(adc);

        for( Int_t ient=0; ient<ent; ient++ ){
            adcEvent->GetEntry(ient);
            if( adc[padc_ch] > 400 ){
                hcal->Fill((Double_t)(adc[padc_ch]));
            }
        }

        hcal->Rebin( rebin );
        cout << "# Fit histogram in :" << rootfile << endl;
        fcnFit = new TF1("fcnFit", TiCuZrFullFitFunc, l_bin, u_bin, nPar_zr );

        fcnFit->SetParNames("FANO", "CONST_N", // parameters 0 ... 5
                            "TiKa1Ch", "TiKa1Gain", "TiKbCh", "RTiKbKa");
        fcnFit->SetParName( 6, "CuKa1Ch ");
        fcnFit->SetParName( 7, "CuKa1Gain ");
        fcnFit->SetParName( 8, "CuKbCh ");
        fcnFit->SetParName( 9, "RCuKbKa ");
        fcnFit->SetParName( 10, "ZrKa1Ch ");
        fcnFit->SetParName( 11, "ZrKa1Gain ");
        fcnFit->SetParName( 12, "ZrKbCh ");
        fcnFit->SetParName( 13, "RZrKbKa "); 
        fcnFit->SetParName( 14, "tShiftZr    ");
        fcnFit->SetParName( 15, "tGainZrKa   ");
        fcnFit->SetParName( 16, "tZrBeta   ");
        fcnFit->SetParName( 17, "tGainZrKb   ");
        fcnFit->SetParName( 18, "tGainTi   ");
        fcnFit->SetParName( 19, "tGainCu   ");
        fcnFit->SetParName( 20, "tTiBeta   ");
        fcnFit->SetParName( 21, "tCuBeta   ");
        fcnFit->SetParName( 22, "tShiftTi    ");
        fcnFit->SetParName( 23, "tShiftCu    ");
        fcnFit->SetParName( 24, "ZrShelfGain    ");
        fcnFit->SetParName( 25, "BremsCh ");
        fcnFit->SetParName( 26, "BremsGain ");
        fcnFit->SetParName( 27, "BremsBeta ");
        fcnFit->SetParName( 28, "BremsSigma ");


         
        fcnFit->SetParameter(0, init_fano );           // FANO
        fcnFit->SetParLimits(0, 0.09, 0.4);
        fcnFit->SetParameter(1, init_cnst );           // Constant noise 
        fcnFit->SetParLimits(1, 0, 100.); 
          
        fcnFit->SetParameter(2, init_mean_ti);         // TiKa1 mean in ch
        fcnFit->SetParLimits(2, tika_lbound, tika_ubound);
        fcnFit->SetParameter(3, ti_gain_setvalue);     // TiKa1 gain
        fcnFit->SetParLimits(3, 100, 1.e7);
        fcnFit->SetParameter(4, init_mean_tikb);       // TiKb mean in ch
        fcnFit->SetParLimits(4, tikb_lbound, tikb_ubound);
        fcnFit->SetParameter(5, 0.26);                  // Ti Kb to Ka gain ratio
        fcnFit->SetParLimits(5, 0.05, 0.8);

        fcnFit->SetParameter(6, init_mean_cu);         // CuKa1 mean in ch
        fcnFit->SetParLimits(6, cuka_lbound, cuka_ubound);
        fcnFit->SetParameter(7, cu_gain_setvalue);    // CuKa1 gain
        fcnFit->SetParLimits(7, 10., 1.e7);
        fcnFit->SetParameter(8, init_mean_cukb);      // CuKb mean in ch
        fcnFit->SetParLimits(8, cukb_lbound, cukb_ubound);
        fcnFit->SetParameter(9, 0.26);                 // Cu Kb to Ka gain ratio
        fcnFit->SetParLimits(9, 0.05, 0.8);  
      
        fcnFit->SetParameter(10, init_mean_zr);         // ZrKa1 mean in ch
        fcnFit->SetParLimits(10, zrka_lbound, zrka_ubound);
        fcnFit->SetParameter(11, zr_gain_setvalue);    // ZrKa1 gain
        fcnFit->SetParLimits(11, 10., 1.e7);
        fcnFit->SetParameter(12, init_mean_zrkb);      // ZrKb mean in ch
        fcnFit->SetParLimits(12, zrkb_lbound, zrkb_ubound);
        fcnFit->SetParameter(13, 0.2);                 // Zr Kb to Ka gain ratio
        fcnFit->SetParLimits(13, 0.05, 0.8);  
        fcnFit->SetParameter(14, 100.);                 // tail shift for zr in eV
        fcnFit->SetParLimits(14, 1., 500.);           
        fcnFit->SetParameter(15, 0.2);                 // Tail gain ratio zr Ka
        fcnFit->SetParLimits(15, 0.001, 0.55);          
        fcnFit->SetParameter(16, 10.);                 // Tail beta for zr
        fcnFit->SetParLimits(16, 1., 100.);  
        fcnFit->SetParameter(17, 0.85);                 // Tail gain ratio zr Kb
        fcnFit->SetParLimits(17, 0.1, 1); 
        fcnFit->SetParameter(18, 0.35);                 // Tail gain ratio Ti Ka1
        fcnFit->SetParLimits(18, 0.001, 0.8);
        fcnFit->SetParameter(19, 0.25);                 // Tail gain ratio Cu Ka1
        fcnFit->SetParLimits(19, 0.001, 0.8);
        fcnFit->SetParameter(20, 20.);                 // Tail beta for Ti
        fcnFit->SetParLimits(20, 1., 100.); 
        fcnFit->SetParameter(21, 20.);                 // Tail beta for Cu
        fcnFit->SetParLimits(21, 1., 100.); 
        fcnFit->SetParameter(22, 100.);                 // tail shift for ti in eV
        fcnFit->SetParLimits(22, 1., 300.);
        fcnFit->SetParameter(23, 100.);                 // tail shift for Cu in eV
        fcnFit->SetParLimits(23, 1., 300.);
        fcnFit->SetParameter(24, 0.0001);                 // zr ka1 shelf 
        fcnFit->SetParLimits(24, 0., 0.1);
        fcnFit->SetParameter(25, init_brems);                 // Bremsstrahlung starting channel
        fcnFit->SetParLimits(25, init_brems - 300., init_brems + 300); 
        fcnFit->SetParameter(26, 1000.);                 // Bremsstrahlungsgain
        fcnFit->SetParLimits(26, 100., 50000.);
        fcnFit->SetParameter(27, 10.);                 // Bremsstrahlungs beta slope
        fcnFit->SetParLimits(27, 1., 100.);
        fcnFit->SetParameter(28, 100.);                 //Bremsstrahlungs sigma
        fcnFit->SetParLimits(28, 0., 1000.);


/*        fcnFit->SetParameter(18, 2500);                 // Brems mean channel
        fcnFit->SetParLimits(18, 2000, 3000);  
        fcnFit->SetParameter(19, 10000.);                 // Brems gain (absolute gain not relative to zr ka1 or so)
        fcnFit->SetParLimits(19, 1., 100000.);  
        fcnFit->SetParameter(20, 5.);                 // Brems Beta
        fcnFit->SetParLimits(20, 0.1, 100.);            
*/

        fcnFit->SetNpx(range);
        cout << "# Npx :" << range << endl;
        cout << "# Set parameters for fit function." << endl;

        c1->cd(1);
        gStyle->SetLabelSize(0.06, "x");
        gStyle->SetLabelSize(0.06, "y");
        gPad->SetTicks();
        gPad->SetGridy();
        gPad->SetGridx();
        gPad->SetLogy();
        hcal->GetXaxis()->SetRangeUser( ll, ul);
        c1->Update();

        status = hcal->Fit("fcnFit", "R"); // here the magic (fit) happens --------
   
        fcnFit->GetParameters(pfit);
        chi2 = fcnFit->GetChisquare();
        ndf  = fcnFit->GetNDF();


        fcnResult = new TF1("fcnResult", TiCuZrFullFitFunc, l_bin, u_bin, nPar_zr );
        fcnResult->SetLineColor(2);
        fcnResult->SetLineWidth(2);
        fcnResult->SetParameters(pfit);
        fcnResult->SetNpx(range);
        fcnResult->Draw("same");

        fcnShelfZr1   = new TF1("fcnShelfZr1",   shelfFunc,  l_bin,  u_bin, 8);

	cout << 1 << endl;

        fcnTika1 = new TF1("fcnTika1", tika1Func, l_bin,  u_bin, 6);
        fcnZrka1 = new TF1("fcnZrka1", zrka1Func, l_bin,  u_bin, 5); // here only 5 par!
        fcnCuka1 = new TF1("fcnCuka1", cuka1Func, l_bin,  u_bin, 7);


        fcnTTiK1 = new TF1("fcnTTiK1", titailFunc,  l_bin,  u_bin, 9);
        fcnTZrK1 = new TF1("fcnTZrK1", zrtailFunc,  l_bin,  u_bin, 9);
        fcnTCuK1 = new TF1("fcnTCuK1", cutailFunc,  l_bin,  u_bin, 10);

        fcnBrems = new TF1("fcnBrems", bremsFunc,  l_bin,  u_bin, 6);


/*
        fcnTTiK2 = new TF1("fcnTTiK2", tibetatailFunc,  l_bin,  u_bin, 9);
        fcnTMnK2 = new TF1("fcnTMnK2", mnbetatailFunc,  l_bin,  u_bin, 9);

        fcnEscTika1 = new TF1("fcnEscTika1", tiescFunc, l_bin, u_bin, 6);
        fcnEscMnka1 = new TF1("fcnEscMnka1", mnescFunc, l_bin, u_bin, 6);
*/

        //gStyle->SetLabelSize(0.06, "x");
        //gStyle->SetLabelSize(0.06, "y");
        
        hcal->GetXaxis()->SetTitle("ADC channel");
        hcal->SetStats(0);


	cout << 2 << endl;
        fcnShelfZr1->SetLineColor(4);                            // blue for shelf
        fcnShelfZr1->SetLineStyle(4);
        fcnShelfZr1->SetLineWidth(2);

	fcnBrems->SetLineWidth(2);				// black for bremsstrahlung
	fcnBrems->SetLineColor(1);


        fcnTika1->SetLineColor(6);                                // Magenda line for Ti
        fcnTika1->SetLineWidth(2);               
        fcnTika1->SetNpx(10000);
        fcnCuka1->SetLineColor(28);                                // brownish line for Cu
        fcnCuka1->SetLineWidth(2);               
        fcnCuka1->SetNpx(10000);
        fcnZrka1->SetLineColor(3);                                // Green for Mn
        fcnZrka1->SetLineWidth(2);               
        fcnZrka1->SetNpx(10000);

	cout << 3 << endl;

        fcnTTiK1->SetLineColor(7);                                // Cyan for tail 
        fcnTTiK1->SetLineWidth(2);
        fcnTTiK1->SetLineStyle(5);
        fcnTTiK1->SetNpx(10000);
/*        fcnTTiK2->SetLineColor(7);                                 
        fcnTTiK2->SetLineWidth(2);
        fcnTTiK2->SetLineStyle(5);
        fcnTTiK2->SetNpx(10000); */
        fcnTZrK1->SetLineColor(7);
        fcnTZrK1->SetLineWidth(2);
        fcnTZrK1->SetLineStyle(5);
        fcnTZrK1->SetNpx(10000);

        fcnTCuK1->SetLineColor(7);
        fcnTCuK1->SetLineWidth(2);
        fcnTCuK1->SetLineStyle(5);
        fcnTCuK1->SetNpx(10000);

	cout << 4 << endl;

/*
        fcnEscTika1->SetLineColor(6);
        fcnEscTika1->SetLineWidth(2);
        fcnEscTika1->SetLineStyle(5);
        fcnEscTika1->SetNpx(10000);
        fcnEscMnka1->SetLineColor(3);
        fcnEscMnka1->SetLineWidth(2);
        fcnEscMnka1->SetLineStyle(5);
        fcnEscMnka1->SetNpx(10000);
*/

        fcnShelfZr1->SetParameters(pfit[2], pfit[10], pfit[0], pfit[1], pfit[24], pfit[10], pfit[11], 1); // Mn Ka1 shelf function

	cout << 5 << endl;



        fcnTika1->SetParameters(pfit[2], pfit[10], pfit[3], pfit[0], pfit[1], 1);
        fcnZrka1->SetParameters(pfit[10], pfit[2], pfit[11],pfit[0], pfit[1]);
        fcnCuka1->SetParameters(pfit[6], pfit[10], pfit[2], pfit[7], pfit[0], pfit[1], 1); // cu ka1 gaussian peak

        fcnTTiK1->SetParameters(pfit[10], pfit[2], pfit[22], pfit[0], pfit[1], pfit[3], pfit[20], pfit[18], 1);
       fcnTZrK1->SetParameters(pfit[10], pfit[2], pfit[14], pfit[0], pfit[1], pfit[11], pfit[16], pfit[15], 1);
        fcnTCuK1->SetParameters(pfit[10], pfit[2], pfit[6], pfit[23], pfit[0], pfit[1], pfit[7], pfit[21], pfit[19], 1); // cupper ka1 tail function

        fcnBrems->SetParameters(pfit[10], pfit[2], pfit[25], pfit[28], pfit[27], pfit[26]);

	cout << 6 << endl;
/* 
        fcnTTiK2->SetParameters(pfit[11], pfit[7], pfit[14], pfit[3], pfit[4], pfit[6], pfit[19], pfit[16], pfit[8]); // Ti Kb tail function; 
        fcnTMnK2->SetParameters(pfit[11], pfit[7], pfit[14], pfit[3], pfit[4], pfit[10],pfit[20], pfit[18], pfit[12]); // Mn Kb tail function



        fcnEscTika1->SetParameters(pfit[9], pfit[5], pfit[6], pfit[3], pfit[4], pfit[13]);
        fcnEscMnka1->SetParameters(pfit[9], pfit[5], pfit[6], pfit[3], pfit[4], pfit[13]);
*/
// ---------  Make the lines and the annotations to the lines ----------

        slope = ( pfit[10] - pfit[2] ) / ( ZrKa1 - TiKa1);
    
        Double_t pXEsctika1, pXEsccuka1, pXEsczrka1;
        pXEsctika1 = pfit[2] - SiKa * slope;
        pXEsczrka1 = pfit[10] - SiKa * slope;
        pXEsccuka1 = pfit[6] - SiKa * slope;
       

        lTika1 = new TLine(pfit[2], 0, pfit[2], 2 * pfit[3]);
        //lTikb  = new TLine(pfit[7], 0, pfit[7], 1.e5);
        lZrka1 = new TLine(pfit[10], 0, pfit[10], 2 * pfit[3]);
        //lMnkb  = new TLine(pfit[11],0, pfit[11],1.e5);
        lCuka1 = new TLine(pfit[6], 0, pfit[6], 2 * pfit[3]);

        lEscTika1 = new TLine(pXEsctika1, 0, pXEsctika1, 2 * pfit[3]);
        lEscZrka1 = new TLine(pXEsczrka1, 0, pXEsczrka1, 2 * pfit[3]);
        lEscCuka1 = new TLine(pXEsccuka1, 0, pXEsccuka1, 2 * pfit[3]);

	cout << 7 << endl;

        lTika1->SetLineStyle(2);
        lTika1->SetLineColor(6);
        lTika1->SetLineWidth(2);
        lZrka1->SetLineStyle(2); 
        lZrka1->SetLineColor(3);
        lZrka1->SetLineWidth(2);
        lCuka1->SetLineStyle(2); 
        lCuka1->SetLineColor(28);
        lCuka1->SetLineWidth(2);

        // Draw the annotations for the fitted lines
        Double_t  tzr = pfit[10] + 60; 
        Double_t  tti = pfit[2] + 40; 
        Double_t  tcu = pfit[6] + 40; 
        //Double_t  tca = pfit[22] + 40; 
        Double_t  ty = pfit[11] / 20.;
        //Double_t  ty_ti = pfit[6] / 10.;
        //Double_t  ty_cu = pfit[23] / 10.;
 
        TLatex *Tzr = new TLatex( tzr, ty, "ZrKa1");
        TLatex *Tti = new TLatex( tti, ty, "TiKa1");
        TLatex *Tcu = new TLatex( tcu, ty, "CuKa1");

        Tzr->SetTextSize(0.04); Tzr->SetTextAngle(90);
        Tti->SetTextSize(0.04); Tti->SetTextAngle(90);
        Tcu->SetTextSize(0.04); Tcu->SetTextAngle(90);

        Tzr->Draw("same");
        Tti->Draw("same");
        Tcu->Draw("same");

	//fcnTMnK2->Draw("same");

       // fcnBck->Draw("same");
        fcnTika1->Draw("same");
        fcnZrka1->Draw("same");
        fcnCuka1->Draw("same");
 //       fcnShelfMn1->Draw("same");
        fcnShelfZr1->Draw("same");
        fcnTTiK1->Draw("same");
        fcnTZrK1->Draw("same");
        fcnTCuK1->Draw("same");
	fcnBrems->Draw("same");
        
 /*     
        fcnEscTika1->Draw("same");
        fcnEscMnka1->Draw("same");*/
        lTika1->Draw("same");
        lZrka1->Draw("same");
        lCuka1->Draw("same");
	lEscTika1->Draw("same");
	lEscCuka1->Draw("same");
	lEscZrka1->Draw("same");


        c1->Update();
        cout << "-----------------------------------" << endl;
        cout << "  chi2 = " << chi2 << ";  NDF = " << ndf << ";  chi2/ndf = " 
             << chi2 / ndf << endl;
        cout << "# Plotted fit functions to pads." << endl;
        
 /*       //////// Tail ana part //////////
        lbtl[0] = pfit[5] - 300.;
        ubtl[0] = pfit[5] + 100.;
        lbtl[1] = pfit[9] - 300.;
        ubtl[1] = pfit[9] + 100.;
        lbpk[0] = pfit[5] - 200.;
        ubpk[0] = pfit[5] + 200.;
        lbpk[1] = pfit[9] - 200.;
        ubpk[1] = pfit[9] + 200.;
        ntiT   = fcnTTiK1->Integral( lbtl[0], ubtl[0] );
        ntiTEr = sqrt( ntiT );
        //ntiTEr = fcnTTiK1->IntegralError( lbtl[0], ubtl[0] );
        ntiK   = fcnTika1->Integral( lbpk[0], ubpk[0] );
        //ntiKEr = fcnTika1->IntegralError( lbpk[0], ubpk[0] );
        ntiKEr = sqrt( ntiK );
        ncuT   = fcnTCuK1->Integral( lbtl[1], ubtl[1] );
        //ncuTEr = fcnTCuK1->IntegralError( lbtl[1], ubtl[1] );
        ncuTEr = sqrt( ncuT );
        ncuK   = fcnCuka1->Integral( lbpk[1], ubpk[1] );
        //ncuKEr = fcnCuka1->IntegralError( lbpk[1], ubpk[1] );
        ncuKEr = sqrt( ncuK );
        ntiK   = ntiK * 1.5;
        ncuK   = ncuK * 1.5;
        ntiKEr = ntiKEr * 1.5;
        ncuKEr = ncuKEr * 1.5;
        if( ntiK > 0 && ncuK > 0 ){
            rtiT = ntiT / ntiK;
            rcuT = ncuT / ncuK; 
            rtiTEr = GetError( ntiT,  ntiTEr, ntiTEr, ntiK,  ntiKEr,  ntiKEr, "DIV");
            rmnTEr = GetError( nmnT,  nmnTEr, nmnTEr, nmnK,  nmnKEr,  nmnKEr, "DIV");
        }else{
            rtiT = 0.;
            rcuT = 0.;
            rtiTEr = 0.;
            rcuTEr = 0.;
        }
        cout << "############ Number of ti tail event : " << ntiT << " +- " << ntiTEr << endl;
        cout << "#            Number of ti peak event : " << ntiK << " +- " << ntiKEr << endl;
        cout << "#            Number of cu tail event : " << ncuT << " +- " << ncuTEr << endl;
        cout << "#            Number of cu peak event : " << ncuK << " +- " << ncuKEr << endl;
        cout << "#            tail to peak ratio ti   : " << rtiT << " +- " << rtiTEr << endl;
        cout << "#            tail to peak ratio cu   : " << rcuT << " +- " << rcuTEr << endl;
 */   }

    return status;
}

void  SDDclass::GetFitParameters( TString source )
{ 

    if( source != "TiCuZr" ){ // start for the case of ti mn source calibration

    cout << "# In SDDclass::GetFitParameters() ... " << endl;
    fano = pfit[3];
    cstn = pfit[4];
    fanoEr = fcnFit->GetParError(3);
    cstnEr = fcnFit->GetParError(4);

    // Default fit function is the TiMn source function
    pXtika1 = pfit[5];
    pXtika1Er = fcnFit->GetParError(5);
    pXmnka1 = pfit[9];
    pXmnka1Er = fcnFit->GetParError(9);
    pXtikb1 = pfit[7];
    pXtikb1Er  = fcnFit->GetParError(7);
//    pXcaka1    = pfit[22];
//    pXcaka1Er  = fcnFit->GetParError(22);
//    pXcakb13   = pfit[24];
//    pXcakb13Er = fcnFit->GetParError(24);
    pXmnkb1 = pfit[11];
    pXmnkb1Er = fcnFit->GetParError(11);

    tikaG   = pfit[6];
    mnkaG   = pfit[10];
    //ttiG    = pfit[18];
    //tcuG    = pfit[19];
    //pileG   = pfit[16];
    //pileGEr = fcnFit->GetParError(16);
    tshift  = pfit[14];
    shelfG  = pfit[21];
    //tBeta   = pfit[];
    
    //pSigR   = pfit[21];
    //pSigREr = fcnFit->GetParError(21);

    if( source == "TiMnCu" )
    {   
         
        
        pXcuka1 = pfit[22];              // For Ti (Mn) Cu source
        pXcuka1Er = fcnFit->GetParError(22);
	pXcukb1 = pfit[24];
	pXcukb1Er = fcnFit->GetParError(24);
        //pXmnka1 = pfit[22];
        //pXmnka1Er = fcnFit->GetParError(22);
        //pXcukb1 = pfit[11];
        //pXcukb1Er = fcnFit->GetParError(11);
        //pXmnkb1 = pfit[23];
        //pXmnkb1Er = fcnFit->GetParError(23);
        //cukaG   = pfit[10];
        //mnkaG   = pfit[13];
        //feG     = pfit[24];
    }

   } // end of the if for both calibrations with ti and mn 

   else { // for the X ray tube calibration - ti and zr as fit fixed points

    cout << "# In SDDclass::GetFitParameters() for Ti + Zr calibration " << endl;
    fano = pfit[0];
    cstn = pfit[1];
    fanoEr = fcnFit->GetParError(0);
    cstnEr = fcnFit->GetParError(1);
   
    pXtika1 = pfit[2];
    pXtika1Er = fcnFit->GetParError(2);
    pXcuka1 = pfit[6];
    pXcuka1Er = fcnFit->GetParError(6);
    pXzrka1 = pfit[10];
    pXzrka1Er = fcnFit->GetParError(10);

    pXtikb1 = pfit[4];
    pXtikb1Er  = fcnFit->GetParError(4);
    pXcukb1 = pfit[8];
    pXcukb1Er = fcnFit->GetParError(8);
    pXzrkb1 = pfit[12];
    pXzrkb1Er = fcnFit->GetParError(12);

  }
    
    return;
}


void  SDDclass::PlotResidue( TString source )
{
    cout << "# In SDDclass::PlotResidue() ... " << endl;
    Double_t lb, ub;
    if( source != "TiCuZr"){ lb = pXtika1 - 200.; ub = pXcuka1 + 200.; }
    if( source == "TiCuZr"){ lb = pXtika1 - 200.; ub = pXzrka1 + 1000.; }    
    Int_t lowBin = hcal->FindBin(lb);
    Int_t upBin  = hcal->FindBin(ub);
    Int_t nBin   = upBin - lowBin + 1;
    Double_t  x[nBin];
    Double_t  y[nBin];
    Double_t  fy[nBin];
    Double_t  res[nBin];
    Double_t  resEr[nBin];
    //Double_t  barP[nBin];
    //Double_t  barM[nBin];
    Double_t  sig[nBin];
    cout << " nBin = " << nBin << endl;
    for( Int_t j = 0; j < nBin; j ++ ){
        x[j]     = hcal->GetBinCenter(lowBin + j);
        y[j]     = hcal->GetBinContent(lowBin + j);
        sig[j]   = sqrt( y[j] );
        fy[j]    = fcnFit->Eval(x[j]); // might be it would be better to take the std deviation of this to calculate the sigma! (bc it is kind of "the deviation of the measured value from the error of the fit"); also possible: the sigma of the residual!!
        if( sig[j] != 0 ){
            res[j]   = ( y[j] - fy[j] ) / sig[j]; // this "compares" the residual to the standard error of the fit -> res = 1 -> 1 standard deviation difference; 68 % of the values should be smaller than +-1 if the model is correct ... http://www.science20.com/quantum_diaries_survivor/those_deceiving_error_bars-85735
        }else{
            res[j]   = 0.;
        }
        resEr[j] = hcal->GetBinError(lowBin + j);
        //barP[j]   = resEr[j] * ( + 1. );
        //barM[j]   = resEr[j] * ( - 1. );
    }

    TGraph *grres = new TGraph(nBin, x, res);
    //TGraph *grbarP= new TGraph(nBin, x, barP);
    //TGraph *grbarM= new TGraph(nBin, x, barM);
    c1->cd(2);
    gStyle->SetLabelSize(0.06, "x");
    gStyle->SetLabelSize(0.06, "y");
    gPad->SetTicks();
    gPad->SetGridx();
    gPad->SetGridy();
    grres->SetTitle("Residual plot");
    grres->GetXaxis()->SetTitle("ADC channel");
    grres->GetXaxis()->SetTitleOffset(0.8);
    grres->GetXaxis()->SetTitleSize(0.06);
    //grres->GetXaxis()->CenterTitle();
    grres->GetYaxis()->SetTitle("Residuals/Sigma");
    grres->GetYaxis()->SetTitleOffset(0.6);
    grres->GetYaxis()->SetTitleSize(0.06);
    grres->GetYaxis()->CenterTitle();
    grres->SetMarkerStyle(7);
    grres->SetMarkerColor(1);
    grres->Draw("AP");
    grres->GetXaxis()->SetRange( 700, 3000 );
    grres->GetXaxis()->SetRangeUser( ll, ul);
    //grbarP->Draw("C");
    //grbarM->Draw("C");
    c1->Update();

    return ;
}


void SDDclass::CalcFwhmMn( TString source )
{
    cout << "# In SDDclass::CalcFwhmMn() ... " << endl;

    // Do not consider about error of sig 
    //Double_t  slope = ( pXcuka1 - pXtika1 ) / ( CuKa1 - TiKa1 );
    //Double_t  slopeEr = sqrt(pXcuka1Er*pXcuka1Er+pXtika1Er*pXtika1Er)/(CuKa1-TiKa1);
    Double_t  slope, slopeEr;
    if( source != "TiCuZr"){ 
	slope = ( pXmnka1 - pXtika1 ) / ( MnKa1 - TiKa1 ); // slope in ch/eV
        slopeEr = sqrt(pXmnka1Er*pXmnka1Er+pXtika1Er*pXtika1Er)/(MnKa1-TiKa1); // slope error propagated from peak position errors (ca 1/100 channel)
        
    }
    else{
	slope = ( pXzrka1 - pXtika1 ) / ( ZrKa1 - TiKa1 ); // slope in ch/eV
	slopeEr = sqrt(pXzrka1Er*pXzrka1Er+pXzrka1Er*pXzrka1Er)/(ZrKa1-TiKa1); // slope error propagated from peak position errors from fit (ca 1/100 channel)
        pXmnka1 = pXtika1 + slope * ( MnKa1 - TiKa1 );
    }

    //sig_tika   = sqrt( pXtika1 * SiW * fano * slope + cstn * cstn );          //[ch]
    sig_mnka   = sqrt( pXmnka1 * SiW * fano * slope + cstn * cstn );          //[ch]
    //sig_tikaEr = 0.1;                   // Only for initialization 
    //sig_mnkaEr = 0.1;                   // Only for initialization 
    fwhm       = SIG2FWHM * ( sig_mnka / slope );     // FWHM for MnKa region [eV] - from fitted fano noise mostly
    fwhmEr     = SIG2FWHM*slopeEr*sig_mnka/slope/slope; //Only propagate slope error, not the error of the sigma (which includes the error of the fano and the const noise)

    //cout << "# Slope, cstn, SiW , pXtika1, Sig_tika, fano , fwhm : " 
    cout << "# Slope: " << slope << ", cstn: " << cstn << ", SiW: " <<  SiW 
         << ", MnKaCh: " << pXmnka1 << ", SigMnKa: " << sig_mnka << ", fano: " 
         <<  fano << ", fwhm (Mn): " << fwhm << ", fwhmEr: " <<  fwhmEr << endl; 
  
    return ;
}


void SDDclass::FitE2ChLine( TString source) // energy to channel line
{
    cout << "# In SDDclass::FitE2ChLine() ... " << endl;
 // s.. standard deviation: s = sqrt((df/dx)^2*s(x)^2+....) 
    // Also plot the linearity of the six ka, kb peaks


    if( source != "TiCuZr" ){// For Ti Mn calibration


    Double_t  slope = ( pXmnka1 - pXtika1 ) / ( MnKa1 - TiKa1 ); // slope in ch/eV
    Double_t  slopeEr = sqrt(pXmnka1Er*pXmnka1Er+pXtika1Er*pXtika1Er)/(MnKa1-TiKa1); // error of slope calculated as sqrt of sum of squares of standard deviations from FIT

    Double_t tikb     = (pXmnka1+pXtika1)/2 + (TiKb1-(TiKa1+MnKa1)/2) * slope; // calculated channel number of ti kb relative to mnka1 and tika1
    Double_t tikbEr   = 1/2*sqrt(sq(pXmnka1Er)+sq(pXtika1Er)+4*sq(TiKb1-(TiKa1+MnKa1)/2)*sq(slopeEr)); // changed here to slopeEr^2 !! sq(x) means squared - defined in common.h
 // apart from that it should be ok
    Double_t tikbChEr = sqrt( sq(tikbEr) + sq(pXtikb1Er) ); // pXtikb1Er is given by fit - when you subtract the 2 errors of the 2 values from each other (for the plot) - the re
// sulting value has this uncertainty
    Double_t mnkb     = (pXmnka1+pXtika1)/2 + (MnKb1-(TiKa1+MnKa1)/2) * slope; 
    Double_t mnkbEr   = 1/2*sqrt(sq(pXmnka1Er)+sq(pXtika1Er)+4*sq(MnKb1-(TiKa1+MnKa1)/2)*sq(slopeEr));
    Double_t mnkbChEr = sqrt( sq(mnkbEr) + sq(pXmnkb1Er) );
    Double_t tika = 0.,     cukb = 0.;
    Double_t tikaEr = 0.,   cukbEr   = 0.;
    Double_t tikaChEr = 0., cukbChEr = 0.;


    pX[0] = TiKa1; pX[1] = TiKb1; pX[2] = MnKa1; pX[3] = MnKb1; pX[4] = 0; pX[5] = 0;
    pY[0] = pXtika1; pY[1] = pXtikb1; pY[2] = pXmnka1; pY[3] = pXmnkb1; pY[4] = 0; pY[5] = 0;
    //pXEr = {0} in SDDclass.h
    pYEr[0] = pXtika1Er; pYEr[1] = pXtikb1Er; pYEr[2] = pXmnka1Er; pYEr[3] = pXmnkb1Er; pYEr[4] = 0; pYEr[5] = 0;
     
    pXmn[0] = TiKa1; pXmn[1] = TiKb1; pXmn[2] = MnKa1; pXmn[3] = MnKb1; pXmn[4] = 0; pXmn[5] =  0;
    pYmn[0] = 0.; pYmn[1] = pXtikb1-tikb; pYmn[2] = 0; pYmn[3] = pXmnkb1-mnkb; pYmn[4] = 0.; pYmn[5] = 0;// difference of fitted channel to where the peak should be according to MnKa1 and TiKa1 positions; value of Kb difference is positive when the fit calculated a value that is too high!!

    // in SDDclass.h Double_t pXmnEr[6] = {0.};    
    pYmnEr[0] = pXtika1Er; pYmnEr[1] = tikbChEr; pYmnEr[2] = pXmnka1Er; pYmnEr[3] = mnkbChEr; pYmnEr[4] = 0.; pYmnEr[5] = 0.;// Ka channel errors given by fit (pXtika1Er) and Kb error given by combined error of fit and kb channel calculation

    pX_line[0] = TiKa1; pX_line[1] = MnKa1;
    pY_line[0] = pXtika1; pY_line[1] = pXmnka1;	

    
   if(source == "TiMnCu" ){
        //slope = ( pXcuka1 - pXmnka1 ) / ( CuKa1 - MnKa1 );
        //slopeEr = sqrt(pXcuka1Er*pXcuka1Er+pXmnka1Er*pXmnka1Er)/(CuKa1-MnKa1);

        //pX[0] = CuKa1;      pY[0] = pXcuka1;     pYEr[0] = pXcuka1Er;
       

        Double_t cuka     = (pXtika1+pXmnka1)/2 + (CuKa1-(MnKa1+TiKa1)/2) * slope; 
        Double_t cukaEr   = 1/2*sqrt(sq(pXtika1Er)+sq(pXmnka1Er)+4*sq(CuKa1-(MnKa1+TiKa1)/2)*sq(slopeEr));
        Double_t cukaChEr = sqrt( sq(cukaEr) + sq(pXcuka1Er) );
	Double_t cukb     = (pXmnka1+pXtika1)/2 + (CuKb1-(TiKa1+MnKa1)/2) * slope;
	Double_t cukbEr   = 1/2*sqrt(sq(pXmnka1Er)+sq(pXtika1Er)+4*sq(CuKb1-(TiKa1+MnKa1)/2)*sq(slopeEr));
	Double_t cukbChEr = sqrt( sq(cukbEr) + sq(pXcukb1Er) );

        pXmn[4] =  CuKa1; pXmn[5] = CuKb1;
        pYmn[4] = pXcuka1 - cuka; pYmn[5] = pXcukb1 - cukb; 
        pYmnEr[4] = cukaChEr; pYmnEr[5] = cukbChEr; 

        pX[4] = CuKa1; pX[5] = CuKb1;
	pY[4] = pXcuka1; pY[5] = pXcukb1;
	pYEr[4] = pXcuka1Er; pYEr[5] = pXcukb1Er;


    }

   } // end of the if statement for the ti mn source cases

  else{ // beginning of ti cu zr x ray tube calibration


    Double_t  slope = ( pXzrka1 - pXtika1 ) / ( ZrKa1 - TiKa1 ); // slope in ch/eV

    Double_t  slopeEr = sqrt(pXzrka1Er*pXzrka1Er+pXtika1Er*pXtika1Er)/(ZrKa1-TiKa1); // error of slope calculated as sqrt of sum of squares of standard deviations from FIT

    Double_t tikb     = (pXzrka1+pXtika1)/2 + (TiKb1-(TiKa1+ZrKa1)/2) * slope; // calculated channel number of ti kb1 relative to zrka1 and tika1
    Double_t tikbEr   = 1/2*sqrt(sq(pXzrka1Er)+sq(pXtika1Er)+4*sq(TiKb1-(TiKa1+ZrKa1)/2)*sq(slopeEr)); // changed here to slopeEr^2 !! sq(x) means squared - defined in common.h
 // apart from that it should be ok
    Double_t tikbChEr = sqrt( sq(tikbEr) + sq(pXtikb1Er) ); // pXtikb1Er is given by fit - when you subtract the 2 errors of the 2 values from each other (for the plot) - the re
// sulting value has this uncertainty

    Double_t zrkb     = (pXzrka1+pXtika1)/2 + (ZrKb1-(TiKa1+ZrKa1)/2) * slope; 
    Double_t zrkbEr   = 1/2*sqrt(sq(pXzrka1Er)+sq(pXtika1Er)+4*sq(ZrKb1-(TiKa1+ZrKa1)/2)*sq(slopeEr));
    Double_t zrkbChEr = sqrt( sq(zrkbEr) + sq(pXzrkb1Er) );

    Double_t cuka     = (pXtika1+pXzrka1)/2 + (CuKa1-(ZrKa1+TiKa1)/2) * slope; //calculated channel number of cu ka1 relative to zr and ti
    Double_t cukaEr   = 1/2*sqrt(sq(pXtika1Er)+sq(pXzrka1Er)+4*sq(CuKa1-(ZrKa1+TiKa1)/2)*sq(slopeEr));
    Double_t cukaChEr = sqrt( sq(cukaEr) + sq(pXcuka1Er) );
    Double_t cukb     = (pXzrka1+pXtika1)/2 + (CuKb1-(TiKa1+ZrKa1)/2) * slope;
    Double_t cukbEr   = 1/2*sqrt(sq(pXzrka1Er)+sq(pXtika1Er)+4*sq(CuKb1-(TiKa1+ZrKa1)/2)*sq(slopeEr));
    Double_t cukbChEr = sqrt( sq(cukbEr) + sq(pXcukb1Er) );

    Double_t tika = 0.;
    Double_t tikaEr = 0.;
    Double_t tikaChEr = 0.;

    

    pX[0] = TiKa1; pX[1] = TiKb1; pX[2] = CuKa1; pX[3] = CuKb1; pX[4] = ZrKa1; pX[5] = ZrKb1;
    pY[0] = pXtika1; pY[1] = pXtikb1; pY[2] = pXcuka1; pY[3] = pXcukb1; pY[4] = pXzrka1; pY[5] = pXzrkb1;
    //pXEr = {0} in SDDclass.h
    pYEr[0] = pXtika1Er; pYEr[1] = pXtikb1Er; pYEr[2] = pXcuka1Er; pYEr[3] = pXcukb1Er; pYEr[4] = pXzrka1Er; pYEr[5] = pXzrkb1Er;
     
    pXmn[0] = TiKa1; pXmn[1] = TiKb1; pXmn[2] = CuKa1; pXmn[3] = CuKb1; pXmn[4] = ZrKa1; pXmn[5] =  ZrKb1;
    pYmn[0] = 0.; pYmn[1] = pXtikb1-tikb; pYmn[2] = pXcuka1-cuka; pYmn[3] = pXcukb1-cukb; pYmn[4] = 0.; pYmn[5] = pXzrkb1-zrkb;// difference of fitted channel to where the peak should be according to MnKa1 and TiKa1 positions; value of Kb difference is positive when the fit calculated a value that is too high!!

    // in SDDclass.h Double_t pXmnEr[6] = {0.};    
    pYmnEr[0] = pXtika1Er; pYmnEr[1] = tikbChEr; pYmnEr[2] = 0.; pYmnEr[3] = 0.; pYmnEr[4] = pXzrka1Er; pYmnEr[5] = zrkbChEr;// Ka channel errors given by fit (pXtika1Er) and Kb error given by combined error of fit and kb channel calculation

    pX_line[0] = TiKa1; pX_line[1] = ZrKa1;
    pY_line[0] = pXtika1; pY_line[1] = pXzrka1;	

  } // end of the if statement for the X ray tube ti zr calibration


    TGraphErrors *gre2cerr = new TGraphErrors(6, pX, pY, pXEr, pYEr ); // graph eV 2 channel with errors
    TGraphErrors *grmn  = new TGraphErrors(6, pXmn, pYmn, pXmnEr, pYmnEr ); // graph with how far from the expected channels some peaks are
    TGraph *gre2c = new TGraph(2,pX_line,pY_line); // graph eV 2 channel

    fe2c = new TF1("fe2c", peakLine, 500, 20000, 2); //  eV 2 ch line
    fe2c->SetNpx(1000);

    ctmp->cd(2);
    gStyle->SetLabelSize(0.05, "x");
    gStyle->SetLabelSize(0.05, "y");
    gPad->SetTicks();
    gPad->SetGridy();
    gPad->SetGridx();
    gre2cerr->SetMarkerStyle(20);
    gre2cerr->SetMarkerColor(4);
    gre2cerr->Draw("AP");
    gre2cerr->SetTitle("ev to channel");
    gre2cerr->GetXaxis()->SetTitle("Energy [eV]");
    gre2cerr->GetXaxis()->SetTitleSize(0.05);
    gre2cerr->GetYaxis()->SetTitle("ADC Channel");
    gre2cerr->GetYaxis()->SetTitleSize(0.05);
    gre2cerr->Draw("AP");
    gre2c->SetMarkerStyle(20);
    gre2c->SetMarkerColor(4);
    gre2c->Fit("fe2c", "R");
    gre2c->Draw("same");

    offset = fe2c->GetParameter(0); // here the offset and the ev2ch ratio are calculated
    ev2ch  = fe2c->GetParameter(1); // in ch/eV
    offsetEr = fe2c->GetParError(0);
    ev2chEr  = fe2c->GetParError(1);
    
    ctmp->cd(1);
    gStyle->SetLabelSize(0.05, "x");
    gStyle->SetLabelSize(0.05, "y");
    gPad->SetTicks();
    gPad->SetGridy();
    gPad->SetGridx();
    grmn ->SetMarkerStyle(20);
    grmn ->SetMarkerColor(4);
    grmn ->SetMarkerSize(1);
    grmn ->GetXaxis()->SetTitle("Energy [eV]");
    grmn ->GetXaxis()->SetTitleSize(0.05);
    grmn ->GetYaxis()->SetTitle("ADC Channel");
    grmn ->GetYaxis()->SetTitleSize(0.05);
    grmn ->SetTitle("Two sources calib linearity");
    grmn ->Draw("AP");
    ctmp->Update();

    //Draw the energy scale histogram
    Double_t bin_ev = 0.;

    froot = new TFile( rootfile, "READ" );

    Double_t bin_number = 3000;
    Double_t low_edge = -0.5;
    Double_t high_edge = -0.5 + bin_number * 1./ev2ch;
    
    TH1F *hev = new TH1F("hev", "Energy spectrum", bin_number, low_edge, high_edge);

    TTree *tr_t = (TTree*)froot->Get("tr");
    Int_t ent = (Int_t)tr_t->GetEntries();

    TBranch *adcEvent   = tr_t->GetBranch("adc");
    adcEvent->SetAddress(adc);


    for( Int_t ient=0; ient<ent; ient++ ){
        adcEvent->GetEntry(ient);
        if( adc[padc_ch] > 200 ){
            bin_ev = ( adc[padc_ch] - offset ) / ev2ch; 
            hev->Fill(bin_ev);
        }
    }
    //hev->Rebin(10);

    cev->cd();
    gStyle->SetLabelSize(0.04, "x");
    gStyle->SetLabelSize(0.04, "y");
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetLogy(); 
    cev->Update();

    hev->GetXaxis()->SetTitle("Energy [eV]");
    hev->GetXaxis()->SetLabelSize(0.04);
    hev->SetStats(0);
    hev->Draw(); //here the energy histogram should be drawn to the canvas
    

   // TFile hist_file("/home/andreas/vip2/reports/1608_VIPReportLNF/sdd1LNGS.root", "RECREATE");
   // hev->Write();  
    //hist_file.close();

    return ;
}


void SDDclass::OpenCanvas()
{
    c1   = new TCanvas("c1",Form("bus%did%02d", bus, id), 
                                       1, 1, 1200, 900); 
    ctmp = new TCanvas("ctmp", "ctmp", 1, 1, 900, 1000);
    cev  = new TCanvas("cev", "energy histogram", 1, 1, 800, 600);

    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
    c1->Divide(1,2);
    ctmp->Divide(1,2);
    return ;
}


void SDDclass::CloseCanvas(Bool_t saveplot)
{
    cout << "# Close canvas c1. " << endl;
    if( saveplot ){                
        c1  ->Print(Form(GRAPH_PATH + "/SDD/calib/" + "%04d%02d%02dCalibFit%02d.pdf",       year, month, day, padc_ch) );
        ctmp->Print(Form(GRAPH_PATH + "/SDD/calib/" + "%04d%02d%02dCalibLinearity%02d.pdf", year, month, day, padc_ch) );
        cev ->Print(Form(GRAPH_PATH + "/SDD/calib/" + "%04d%02d%02dEnergySpectrum%02d.pdf", year, month, day, padc_ch) );
    }
    //c1->Close();
    //ctmp->Close();
    //cev->Close();
    return ;
}

void SDDclass::CloseRootFile()
{
    cout << "# Close TFile froot. " << endl;
    froot->Close();
    return ;
}

void SDDclass::GetHistoRange()
{
    cout << "# id,  padc : " << id << ", "  << padc_ch << endl;
    ll = 350., ul = 1500.; 

    cout << "# Range of histogram set to be " << ll << " - " << ul << endl;
  return;
}


Time   SDDclass::GetRunTime()
{
    cout << "# In GetRunTime() ..." << froot << endl;
    cout << "# Time span of the run " << tdiff << " sec." << endl;
    return tdiff;
}

/*
AnaPar SDDclass::GetAverageRate()
{ 
    return rate; 
}




ID SDDclass::GetBus()  const { return bus; }
ID SDDclass::GetID()   const { return id;  }
ID SDDclass::GetChip() const { return chip; }
ID SDDclass::GetGrp()  const { return grp; }
ID SDDclass::GetPos()  const { return pos; }
ID SDDclass::GetRow()  const { return row; }
ID SDDclass::GetCol()  const { return col; }
ID SDDclass::GetRnk()  const { return rnk; }
*/
Time   SDDclass::GetTlabel() const { return tlabel; }
Time   SDDclass::GetUnixT()  const { return ut; }

Int_t  SDDclass::GetNdf()       const { return ndf; }
AnaPar SDDclass::GetChi2()      const { return chi2; }
AnaPar SDDclass::GetFano()      const { return fano; }
AnaPar SDDclass::GetFanoEr()    const { return fanoEr; }
AnaPar SDDclass::GetCstn()      const { return cstn; }
AnaPar SDDclass::GetCstnEr()    const { return cstnEr; }
AnaPar SDDclass::GetTika1()     const { return pXtika1; }
AnaPar SDDclass::GetTika1Er()   const { return pXtika1Er; }
AnaPar SDDclass::GetCuka1()     const { return pXcuka1; }
AnaPar SDDclass::GetCuka1Er()   const { return pXcuka1Er; }
AnaPar SDDclass::GetMnka1()     const { return pXmnka1; }
AnaPar SDDclass::GetMnka1Er()   const { return pXmnka1Er; }
AnaPar SDDclass::GetMnD()       const { return mnD; }
AnaPar SDDclass::GetMnDEr()     const { return mnDEr; }
AnaPar SDDclass::GetCuKaG()     const { return cukaG; }
AnaPar SDDclass::GetTiTailG()   const { return ttiG; }
AnaPar SDDclass::GetCuTailG()   const { return tcuG; }
AnaPar SDDclass::GetPileG()     const { return pileG; }
AnaPar SDDclass::GetPileGEr()   const { return pileGEr; }
AnaPar SDDclass::GetTMnBeta()     const { return tMnBeta; }
AnaPar SDDclass::GetTTiBeta()     const { return tTiBeta; }
AnaPar SDDclass::GetTShift()    const { return tshift; }
AnaPar SDDclass::GetTiKaG()     const { return tikaG; }
AnaPar SDDclass::GetMnKaG()     const { return mnkaG; }
AnaPar SDDclass::GetFeG()       const { return feG; }
AnaPar SDDclass::GetSigTika()   const { return sig_tika; }
AnaPar SDDclass::GetSigTikaEr() const { return sig_tikaEr; }
AnaPar SDDclass::GetFwhm()      const { return fwhm; }
AnaPar SDDclass::GetFwhmEr()    const { return fwhmEr; }
AnaPar SDDclass::GetEv2ch()     const { return ev2ch; }
AnaPar SDDclass::GetEv2chEr()   const { return ev2chEr; }
AnaPar SDDclass::GetE2ChOffset()   const { return offset; }
AnaPar SDDclass::GetE2ChOffsetEr() const { return offsetEr; }
AnaPar SDDclass::GetNtiK()       const { return ntiK; }
AnaPar SDDclass::GetNtiKEr()     const { return ntiKEr; }
AnaPar SDDclass::GetNcuK()       const { return ncuK; }
AnaPar SDDclass::GetNcuKEr()     const { return ncuKEr; }
AnaPar SDDclass::GetNtiT()       const { return ntiT; }
AnaPar SDDclass::GetNtiTEr()     const { return ntiTEr; }
AnaPar SDDclass::GetNcuT()       const { return ncuT; }
AnaPar SDDclass::GetNcuTEr()     const { return ncuTEr; }
AnaPar SDDclass::GetPSigR()      const { return pSigR;  }
AnaPar SDDclass::GetPSigREr()    const { return pSigREr;}

/*
SlowPar SDDclass::GetVrefPOS()   const { return vrefPOS; }
SlowPar SDDclass::GetVrefPULSE() const { return vrefPULSE; }
SlowPar SDDclass::GetDacOffset() const { return dacOffset; }
SlowPar SDDclass::GetVsa()  const { return vsa; }
SlowPar SDDclass::GetVsss() const { return vsss; }
SlowPar SDDclass::GetRxA()  const { return rxA; }
SlowPar SDDclass::GetRxB()  const { return rxB; }
SlowPar SDDclass::GetBfA()  const { return bfA; }
SlowPar SDDclass::GetBfB()  const { return bfB; }
SlowPar SDDclass::GetVb1A() const { return vb1A;}
SlowPar SDDclass::GetVb1B() const { return vb1B;}
SlowPar SDDclass::GetVb2A() const { return vb2A;}
SlowPar SDDclass::GetVb2B() const { return vb2B;}
SlowPar SDDclass::GetVb3A() const { return vb3A;}
SlowPar SDDclass::GetVb3B() const { return vb3B;}
SlowPar SDDclass::GetR1A()  const { return r1A; }
SlowPar SDDclass::GetR1B()  const { return r1B; }
SlowPar SDDclass::GetP12VA() const { return p12VA; }
SlowPar SDDclass::GetVdd()   const { return vdd; }
SlowPar SDDclass::GetM12VA() const { return m12VA; }
SlowPar SDDclass::GetVss()   const { return vss; }
SlowPar SDDclass::GetP12VB() const { return p12VB; }
SlowPar SDDclass::GetHv()    const { return hv;  }
SlowPar SDDclass::GetM12VB() const { return m12VB; }
SlowPar SDDclass::GetM5V()   const { return m5V; }
SlowPar SDDclass::GetTsenseA()   const { return tsenseA; }
SlowPar SDDclass::GetTsenseB()   const { return tsenseB; }
SlowPar SDDclass::GetPt100()     const { return pt100; }
SlowPar SDDclass::GetTADC()      const { return tADC; }
SlowPar SDDclass::GetXtubeAm()   const { return xtubeAm; }
*/
