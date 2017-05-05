/************************************
  ReadFullDaqBinData.C
  For the NI daq at LNF. 
  2013/10/27
  For CAEN V792 QDC at this moment. 
  2013/12/16 
  modified for binary data reading.
  2015/05/29 
  modified for 32 scinti 6 sdd binary data.
  H.Shi
************************************/
#include  <stdlib.h>
#include  <stdio.h>
#include  <iostream>

#include  "TStyle.h"
#include  "TObject.h"
#include  "TDatime.h"
#include  "TFile.h"
#include  "TTree.h"
#include  "TBranch.h"
#include  "TCanvas.h"
#include  "TH1.h"
#include  "TH2.h"
#include  "TAxis.h"

#include  "common_Kaku.h"
#include  "PhaseOneScintiConnectionTable.h"
#include  "ReadFullDaqBinData_Kaku.h" 
//#include  "SDDcalib.h"
//#include  "qdc_ana.C"


using namespace std;


void MakeTree( TString filelist, TString rootfilename ) 
{

    Int_t  read_status = - 999;

    Int_t  qdc_th[NScinti] = { 0 };

    EventStruct  evt; 
    TString rootfile; 
    
    FILE  *flist;
    if( (flist = fopen(filelist, "r")) == NULL ){
        cout << "Cannot open file: " << filelist << endl;
    }
    cout << filelist << endl;

    char  listline[MAXCHAR];
    char  file_name[MAXCHAR];

    rootfile = ROOT_PATH + "/" + rootfilename; 
    TFile *f = new TFile( rootfile, "RECREATE");
    f->cd();

    TTree *tr = new TTree("tr", "SiPM and SDD data");
    tr->Branch("evid",   &evt.evid,   "Event number/I");
    tr->Branch("evid1",  &evt.evid1,  "Event number qdc/S");
    tr->Branch("evid2",  &evt.evid2,  "Event number padc/S");
    tr->Branch("evid3",  &evt.evid3,  "Event number tdc/S");
    tr->Branch("ut",     &evt.ut,     "Unix time tag/I");
    tr->Branch("rate",   &evt.rate,   "Trigger rate/S");
    tr->Branch("mul",    &evt.mul,    "Event multiplicity/S");
    tr->Branch("qdc",    evt.qdc,     "QDC Channel[32]/S");
    tr->Branch("tdc",    evt.tdc,     "TDC Channel[16]/S");
    tr->Branch("adc",    evt.padc,    "PeakADC[16]/S");
    tr->Branch("layer",  &evt.layer,  "Scinti layer/S");
    tr->Branch("col",    &evt.col,    "Scinti col/S");
    tr->Branch("trgid",  &evt.trgid,  "Trigger type/S");
    tr->Branch("hr",     &evt.hr,     "High Rate/S");
    tr->Branch("clk",    &evt.clk,    "LV ms clock/I");

    TH1F* ht[ENABLED_TDC_CH];   // TDC 
    TH1F* hp[ENABLED_PADC_CH];  // PADC
    TH1F* he[ENABLED_PADC_CH];  // Energy histogram
    TH1F* hq[NScinti];          // QDC

    TH2F* hpos;                 // Hit position

    for( Int_t i=0; i < ENABLED_TDC_CH; i++ ){
        ht[i] = new TH1F(Form("ht[%d]", i), "Timing spectra", 2048, -200.5, 1847.5);
    }

    for( Int_t i=0; i < ENABLED_PADC_CH; i++ ){
        hp[i] = new TH1F(Form("hp[%d]", i), "PADC spectra", 512, -0.5, 4095.5);
    }

    for( Int_t i=0; i < ENABLED_PADC_CH; i++ ){
        he[i] = new TH1F(Form("he[%d]", i), "SDD Energy spectra", 32768, -0.5, 32767.5);
    }

    for( Int_t i=0; i < NSiPM; i++ ){
        hq[i] = new TH1F(Form("hq[%d]", i), "QDC spectra", 512, -0.5, 4095.5);
    }

    hpos = new TH2F( "hpos", "Hit position at Veto counters", 8, -0.5, 7.5, 6, -0.5, 5.5);

    f->Write();
    f->Close();

    while( fgets(listline, MAXCHAR, flist) != NULL)
    {
        cout << file_name << endl;
        sscanf( listline, "%s\n", file_name );

        read_status = ReadData( file_name, rootfile, qdc_th);

        cout << "Return value : " << read_status << endl;
    }
    
    MakePlots( rootfile, qdc_th );

    fclose(flist);

    return;
}



Int_t ReadData( TString datafilename, TString root_file, Int_t qdc_th[NScinti] )
{

    Int_t     rtn = -999; 

    Short_t   mul = 0;
    Short_t   dummyword = -1;
    Short_t   tdc_tmp[3] = {0};

    UShort_t   clk_db = 0,  clk_ub = 0;

    EventStruct  evt_r; 

    evt_r.evid  = 0,     evt_r.evid1 = 0,  evt_r.evid2 = 0,  evt_r.evid3 = 0;
    evt_r.ut    = 0,     evt_r.rate  = -99;   
    evt_r.year  = 2014,  evt_r.month = 0;  
    evt_r.day   = 0,     evt_r.hour  = 0,   evt_r.min = 0;
    evt_r.layer = 0,     evt_r.col   = 0;  
    evt_r.trgid = -1,    evt_r.hr    = 0;
    evt_r.clk   = 0;

    // Histograms for timing
    TH1F* ht;

    // Histograms for SiPMs
    TH1F* hvq;
    TH2F* hpos;

    // Histograms for SDDs
    TH1F* hsp;

    TFile *ff = new TFile( root_file, "UPDATE");
    ff->cd();

    TTree *tree = (TTree*)(ff->Get("tr"));

    TBranch *evidEvent  = tree->GetBranch("evid");
    TBranch *evid1Event = tree->GetBranch("evid1");
    TBranch *evid2Event = tree->GetBranch("evid2");
    TBranch *evid3Event = tree->GetBranch("evid3");
    TBranch *utEvent    = tree->GetBranch("ut");
    TBranch *rateEvent  = tree->GetBranch("rate");
    TBranch *mulEvent   = tree->GetBranch("mul");
    TBranch *qdcEvent   = tree->GetBranch("qdc");
    TBranch *tdcEvent   = tree->GetBranch("tdc");
    TBranch *padcEvent  = tree->GetBranch("adc");
    TBranch *layerEvent = tree->GetBranch("layer");
    TBranch *colEvent   = tree->GetBranch("col");
    TBranch *trgidEvent = tree->GetBranch("trgid");
    TBranch *hrEvent    = tree->GetBranch("hr");
    TBranch *clkEvent   = tree->GetBranch("clk");


    evidEvent  ->SetAddress(&evt_r.evid);
    evid1Event ->SetAddress(&evt_r.evid1);
    evid2Event ->SetAddress(&evt_r.evid2);
    evid3Event ->SetAddress(&evt_r.evid3);
    utEvent    ->SetAddress(&evt_r.ut);
    rateEvent  ->SetAddress(&evt_r.rate);
    mulEvent   ->SetAddress(&evt_r.mul);
    qdcEvent   ->SetAddress(evt_r.qdc);
    tdcEvent   ->SetAddress(evt_r.tdc);
    padcEvent  ->SetAddress(evt_r.padc);
    layerEvent ->SetAddress(&evt_r.layer);
    colEvent   ->SetAddress(&evt_r.col);
    trgidEvent ->SetAddress(&evt_r.trgid);
    hrEvent    ->SetAddress(&evt_r.hr);
    clkEvent   ->SetAddress(&evt_r.clk);


    FILE *datafile;

    // binary data for Phaes One setup in the lab: 
    TString binary_file = BIN_PATH + "/PhaseOne/" + datafilename; 
    ///TString binary_file = BIN_PATH + "/" + datafilename; 

    if((datafile = fopen(binary_file, "r"))==NULL){
        cout << "Cannot open file: " << binary_file << endl;
    }

    cout << binary_file << endl;

    // Get the time tag of the run
    sscanf ( datafilename, "%04hd%02hd%02hd_%02hd%02hd", &evt_r.year, &evt_r.month, 
                           &evt_r.day, &evt_r.hour, &evt_r.min );

    TDatime *datime = new TDatime( evt_r.year, evt_r.month, evt_r.day, 
                                   evt_r.hour, evt_r.min,   0 );
    evt_r.ut   = datime->Convert();
    //datime->Print();
    //evt_r.time = evt_r.year*1e10 + evt_r.month*1e8 + evt_r.day*1e6 
    //             + evt_r.hour*1e4 + evt_r.min*1e2; 
    delete datime;


    // Read the binary data
    while(1)
    {
        // new line begins with xx 00 in hex 
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 )
        {
            if ( feof (datafile) ){
                cout << "Reached end of file " << endl;
                break;
            }else {
                cout << "Unexpected termination ! " << endl;
                break;
            }
        }else{
            //cout << " begin of line : " << dummyword << endl;
        }

        evt_r.evid ++;

        // header 0
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;

        // daq rate 
        if ( fread(&evt_r.rate, sizeof(evt_r.rate), 1, datafile) != 1 ) break;

        // evid from qdc  
        if ( fread(&evt_r.evid1, sizeof(evt_r.evid1), 1, datafile) != 1 ) break;     // comment out if only PADC data
        //cout << " evid of qdc : " << evt_r.evid1  << endl;

        // V792_QDC_CH coloumns for qdc
        if ( fread(&evt_r.qdc, sizeof(evt_r.qdc[0]), V792_QDC_CH, datafile) != V792_QDC_CH ) break;
        // for Phase One 
        ///if ( fread(&evt_r.qdc, sizeof(evt_r.qdc[0]), NSiPM, datafile) != NSiPM) break;   // comment out if only PADC data

        // 0 column for 32 ch qdc and -1 division
        ///if ( fread(&dummyword, sizeof(dummyword), 2, datafile) != 2 ) break;      // if disabled 
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;      // only -1 column if enabled

        // evid from padc  
        if ( fread(&evt_r.evid2, sizeof(evt_r.evid2), 1, datafile) != 1 ) break;     // comment out if no PADC data
        ///cout << " evid of padc : " << evt_r.evid2  << endl;

        for( Int_t i=0; i < V785_ADC_CH; i++){
            evt_r.padc[i] = 0;  
        }

        // V785_ADC_CH coloumns for ph adc for Phase One binary data 
        // number of enabled channels NPADC_CH is controlled from the DAQ  
        if ( fread(&evt_r.padc, sizeof(evt_r.padc[0]), ENABLED_PADC_CH, datafile) != ENABLED_PADC_CH ) break;
        if ( fread(&dummyword, sizeof(dummyword), 2, datafile) != 2 ) break;  // 0, and -1 division coloumn afterwards;

        // evid from tdc 
        if ( fread(&evt_r.evid3, sizeof(evt_r.evid3), 1, datafile) != 1 ) break;
        //cout << " evid of tdc : " << evt_r.evid3  << endl;

        /// /*         enable this comment out if only PADC data
        for( Int_t i=0; i < ENABLED_TDC_CH; i++){
            evt_r.tdc[i] = -10;
        }

        // read the tdc data, empty channel has an entry of -31073 
        for( Int_t i=0; i < ENABLED_TDC_CH ; i++){
            if( fread ( &tdc_tmp[0], sizeof(tdc_tmp[0]), 1, datafile ) != 1 ) break;
            if( fread ( &tdc_tmp[1], sizeof(tdc_tmp[1]), 1, datafile ) != 1 ) break;
            if( fread ( &tdc_tmp[2], sizeof(tdc_tmp[2]), 1, datafile ) != 1 ) break;
            if( tdc_tmp[0] < -30000){      // after 10.04.2014 the out of range events are the limit of Short_t
                tdc_tmp[0] = -100;
                tdc_tmp[2] = -10;
            }else{
                mul ++ ;
            }

            if( tdc_tmp[1] == i ){
                evt_r.tdc[tdc_tmp[1]] = tdc_tmp[0]; 
            }
            evt_r.mul = mul;
        }

        // Fill the histograms
        Float_t  tdc    = .0;

        // Histograms for TDC 
        // Reference TDC channel tdc[14],  SDD OR ch[12],  scinti trigger ch[13].
        // the SDD timing is given by tdc[12] - tdc[14]; 
        for( Int_t i=0; i < ENABLED_TDC_CH ; i++){
            ht = (TH1F*)ff->Get( Form( "ht[%d]", i) );
            if( evt_r.tdc[i] > 0 ) {
                tdc = ( evt_r.tdc[TDC_REF_CH] - evt_r.tdc[i] ) * NSECTDCCH ;
                //cout << "tdc ch : " << i << ";  tdc value : "  << tdc << "; " << endl;
                ht->Fill( tdc ); 
            }
        }

        // Fill the histograms for SiPMs
        // qdc
        hpos = (TH2F*)ff->Get( "hpos" );
        for( Int_t j=0; j < NSiPM; j++ )
        {
            hvq = (TH1F*)ff->Get( Form("hq[%d]", j) );
            hvq->Fill( evt_r.qdc[j] );

            // Fill hit position
            // Threshold value 
            if( qdc_th[j] == 0 ){
                qdc_th[j] = 700;
            }
            if( evt_r.qdc[j] > qdc_th[j] ){
                evt_r.layer = Qdc2Layer[j];
                evt_r.col   = Qdc2Column[j];
                hpos->Fill( evt_r.col, evt_r.layer );
            }
        }

        /////////////////////////////////
        /// */  // enable this in case of only PADC data

        // Fill PADC histograms
        for( Int_t i=0; i < ENABLED_PADC_CH ; i++)
        {
            hsp = (TH1F*)ff->Get( Form("hp[%d]", i) );
            if( evt_r.padc[i] > 50 ){
                hsp->Fill( evt_r.padc[i] );
            }
        }

        //  0 columns, tdc word count, and -1 division 
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        //if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        if ( fread(&clk_ub, sizeof(clk_ub), 1, datafile) != 1 ) break;
        if ( fread(&clk_db, sizeof(clk_db), 1, datafile) != 1 ) break;

        evt_r.clk = clk_ub << 16 ^ clk_db;

        //cout << clk_ub << ",    " << clk_db << endl; 
        //cout << evt_r.clk << endl;

        // Fill Tree
        tree->Fill();

        mul = 0;
    }


    if ( fclose(datafile) == 0){
        rtn = 0;
    }else{
        rtn = -1;
    }

    ff->Write("", TObject::kOverwrite);
    ff->Close();

    cout << "Finished run : " << binary_file << endl;

    return rtn; 
}


////////////////////////////////////////////
// Can be directly called if root file exists
void MakePlots(TString rootfile, Int_t qdc_th[NScinti] )
{
    // Plot the qdc, tdc, padc 

    TFile *f = new TFile(rootfile, "READ");
    f->cd();

    TH1F* ht;
    TH1F* hmul; 
    TH1F* hadc;
    //TH1F* hev;
    TH1F* hvq;
    TH2F* hpos;

    TString graph_title;

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    cc->cd();
    cc->SetGridx();
    cc->SetGridy();

    for( Int_t i=0; i<ENABLED_TDC_CH; i++){
        ht = (TH1F*)f->Get( Form("ht[%d]", i) );
        ht->GetXaxis()->SetTitle("TDC [ns]");
        ht->Draw();
        cc->SetLogy();
        cc->Print( Form( GRAPH_PATH + "/TDC_CH%02d.pdf", i ) );
    }

    for( Int_t i=0; i<ENABLED_PADC_CH; i++){
        hadc = (TH1F*)f->Get( Form("hp[%d]", i) );
        hadc->GetXaxis()->SetTitle("PADC [ch]");
        hadc->GetXaxis()->SetRangeUser(200, 2500);
        hadc->Draw();
        cc->SetLogy(0);
        cc->Print( Form( GRAPH_PATH + "/PADC_CH%02d.pdf", i ) );
    }

    for( Int_t i=0; i<NSiPM; i++ ){
        hvq  = (TH1F*)f->Get( Form("hq[%d]", i ) );
        hvq->GetXaxis()->SetTitle("QDC [ch]");
        //hvq->Draw();
        //cc->SetLogy();
       
        qdc_th[i] = FitQdcSpectrum(hvq, cc);
        /*
        if( i < 10 || i == 11 || i >23 ){     // For 201505 data set
            qdc_th[i] = FitQdcSpectrum(hvq, cc);
        }else{
            cc->cd();
            hvq->Draw();
            cc->SetLogy();
        }
        */
        cc->Print( Form( GRAPH_PATH + "/QDC_SiPM_CH%02d.pdf", i) );
    }


    // Fill and plot again the hit position 2d plot
    EventStruct  evt_p; 
    evt_p.mul = 0;

    TTree *tree = (TTree*)(f->Get("tr"));

    TBranch *evidEvent  = tree->GetBranch("evid");
    TBranch *qdcEvent   = tree->GetBranch("qdc");
    TBranch *adcEvent   = tree->GetBranch("adc");
    TBranch *tdcEvent   = tree->GetBranch("tdc");
    TBranch *layerEvent = tree->GetBranch("layer");
    TBranch *colEvent   = tree->GetBranch("col");

    evidEvent  ->SetAddress(&evt_p.evid);
    qdcEvent   ->SetAddress(evt_p.qdc);
    adcEvent   ->SetAddress(evt_p.padc);
    tdcEvent   ->SetAddress(evt_p.tdc);
    layerEvent ->SetAddress(&evt_p.layer);
    colEvent   ->SetAddress(&evt_p.col);
    

    // Reset the bin content of the hit position histogram;
    Int_t entries = tree->GetEntries();
    hpos= (TH2F*)f->Get( "hpos" );

    for( Int_t icol = 0; icol < 8; icol ++ ){
        for( Int_t ilayer = 0; ilayer < 6; ilayer ++ ){
            hpos->SetBinContent( icol, ilayer, 0.);
        }
    }

    hmul = new TH1F("hmul", "scinti multiplicity", 100, -0.5, 99.5);

    // Refill the hit postiion histogram
    for( Int_t ient=0; ient < entries; ient++ ){
        for( Int_t j=0; j < NSiPM; j++ )
        {
            qdcEvent    ->GetEntry( ient );
            layerEvent  ->GetEntry( ient );
            colEvent    ->GetEntry( ient );

            // Fill hit position
            // Threshold value 
            if( qdc_th[j] == 0 ){
                qdc_th[j] = 1000;
            }
            if( evt_p.qdc[j] > qdc_th[j] ){
                evt_p.layer = Qdc2Layer[j];
                evt_p.col   = Qdc2Column[j];
                hpos->Fill( evt_p.col, evt_p.layer );

                evt_p.mul ++ ;
            }
        }
        hmul->Fill( evt_p.mul);
        evt_p.mul = 0;
    }


    /// event mismatch of plus minus one event id check
    TH1F* hevt  = new TH1F("hevt", "adc all", 512, 0.5, 4096.5);
    TH1F* hac   = new TH1F("hac", "adc adc cut", 512, 0.5, 4096.5);
    TH1F* htc   = new TH1F("htc", "adc tdc cut", 512, 0.5, 4096.5);
    TH1F* hatc  = new TH1F("hatc", "adc tdc mp one event cut", 512, 0.5, 4096.5);

    Bool_t    adc_flag  = 0;
    Bool_t    adc_flag1 = 0;
    Bool_t    adc_flag2 = 0;
    Bool_t    tdc_flag  = 0;
    Bool_t    tdc_flag1 = 0;
    Bool_t    tdc_flag2 = 0;

    for( Int_t ient=0; ient < entries; ient++ ){
        adcEvent    ->GetEntry( ient );
        tdcEvent    ->GetEntry( ient );

        for( Int_t k=0; k < 1; k++ )
        {

            if( evt_p.padc[k] > 300 ){

                adc_flag = 1;

                hevt->Fill( (Double_t)(evt_p.padc[k]) );
                hac->Fill( (Double_t)(evt_p.padc[k]) );

                if( evt_p.tdc[12] > 0 ){
                    htc->Fill( (Double_t)(evt_p.padc[k]) );
                    hatc->Fill( (Double_t)(evt_p.padc[k]) );

                    tdc_flag = 1;
                }
            }
        }

        if( adc_flag == 1 && tdc_flag == 0 && ient > 0 && ient < entries )
        {

            // check previous event 
            adcEvent    ->GetEntry( ient - 1 );
            tdcEvent    ->GetEntry( ient - 1 );

            if( evt_p.tdc[12] > 0 )tdc_flag1 = 1; 

            for( Int_t k=0; k < 1; k++ )
            {
                if( evt_p.padc[k] > 300 ){
                    adc_flag1 = 1;
                }
            }

            // check following event 
            adcEvent    ->GetEntry( ient + 1 );
            tdcEvent    ->GetEntry( ient + 1 );

            if( evt_p.tdc[12] > 0 )tdc_flag2 = 1; 

            for( Int_t k=0; k < 1; k++ )
            {
                if( evt_p.padc[k] > 300 ){
                    adc_flag2 = 1;
                }
            }
        }

        // Fill the hatc histo if pre or next event has no adc entry but tdc entry 
        if( (adc_flag1 == 0 && tdc_flag1 == 1 && adc_flag == 1) || 
                (adc_flag2 == 0 && tdc_flag2 == 1 && adc_flag == 1) ){

            adcEvent    ->GetEntry( ient );
            for( Int_t k=0; k < 1; k++ )
            {
                if( evt_p.padc[k] > 300 ){
                    hatc->Fill( (Double_t)(evt_p.padc[k]) );
                }
            }
        }
        
        adc_flag  = 0;
        adc_flag1 = 0;
        adc_flag2 = 0;
        tdc_flag  = 0;
        tdc_flag1 = 0;
        tdc_flag2 = 0;
    }

    cc->cd();
    hpos->GetXaxis()->SetTitle("Column");
    hpos->GetYaxis()->SetTitle("Layer");
    hpos->Draw("BOX,TEXT");
    cc->SetLogy(0);
    graph_title = GRAPH_PATH + "/Veto_Counter_hit_pos_distribution.pdf";
    cc->Print( graph_title );

    cc->cd();
    hmul->GetXaxis()->SetTitle("scinti multiplicity");
    hmul->GetXaxis()->SetRangeUser(0, 32);
    hmul->Draw();
    cc->SetLogy();
    graph_title = GRAPH_PATH + "/Veto_Counter_multiplicity.pdf";
    cc->Print( graph_title );

    cc->cd();
    cc->SetLogy(0);
    hevt->GetXaxis()->SetRangeUser(200, 2000);
    hevt->Draw();
    hac->Draw("same");
    htc->SetLineColor(2);
    htc->Draw("same");
    graph_title = GRAPH_PATH + "/Evid_matching.pdf";
    cc->Print( graph_title );
    hatc->GetXaxis()->SetRangeUser(200, 2000);
    hatc->Draw();
    hatc->SetLineColor(2);
    hevt->Draw("same");
    graph_title = GRAPH_PATH + "/Evid_matching_pm_events.pdf";
    cc->Print( graph_title );

    cc->Close();

    f->Close();

    return; 
}

