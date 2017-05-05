/************************************
  ReadFullDaqBinData.C

  modified for binary data reading.
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

#include  "common.h"
#include  "ReadFullDaqBinData.h" 
#include  "PhaseOneScintiConnectionTable.h"


using namespace std;
//gROOT->Reset();

//ofstream hannes;

void MakeTree( TString filelist, TString rootfilename ) 
{

    Int_t  read_status = - 999;

    EventStruct  evt; // defined in header
    TString rootfile; 
    
    FILE  *flist;
    if( (flist = fopen(BIN_PATH + "/" +filelist, "r")) == NULL ){
        cout << "Cannot open file: " << filelist << endl;
    }
    cout << BIN_PATH + "/" +filelist << endl;

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
    tr->Branch("ut",     &evt.ut,     "Unix time tag/i");
    tr->Branch("rate",   &evt.rate,   "Trigger rate/S");
    tr->Branch("mul",    &evt.mul,    "Event multiplicity/S");
    tr->Branch("qdc",    evt.qdc,     "QDC Channel[32]/S");
    tr->Branch("tdc",    evt.tdc,     "TDC Channel[64]/S");
    tr->Branch("tot",    evt.tot,     "TimeOverThreshold[64]/S");
    tr->Branch("padc",   evt.padc,    "PeakADC[32]/S");
    //tr->Branch("layer",  evt.layer,   "Scinti layer/S");
    //tr->Branch("col",    evt.col,     "Scinti col/S");
    tr->Branch("trgid",     evt.trgid,      "Trigger type/S");
    tr->Branch("hr",        evt.hr,         "High Rate/S");

    TH1F* ht[V1190_TDC_CH];   // TDC 

    TH1F* hp[ENABLED_NPADC_CH];       // PADC
    TH1F* hq[NScinti];    // QDC
    TH1F* hqc[NScinti];   // QDC with tot cut

    TH2F* hqt[NScinti];   // QDC - ToT correlation

    TH2F* hpos;         // Hit position

    for( Int_t i=0; i < V1190_TDC_CH; i++ ){
        ht[i] = new TH1F(Form("ht[%d]", i), "Timing spectra", 2048, -40000, 20000);
    }

    for( Int_t i=0; i<ENABLED_NPADC_CH; i++ ){
        hp[i] = new TH1F(Form("hp[%d]", i), "PADC spectra", 2048, -0.5, 4095.5);
    }

    for( Int_t i=0; i<NScinti; i++ ){
        hq[i] = new TH1F(Form("hq[%d]", i), "QDC spectra", 512, -0.5, 4095.5);
    }

    for( Int_t i=0; i<NScinti; i++ ){
        hqc[i] = new TH1F(Form("hqc[%d]", i), "QDC spectra with tot cut", 512, -0.5, 4095.5);
    }

    for( Int_t i=0; i<NScinti; i++ ){
        hqt[i] = new TH2F( Form("hqt[%d]", i), "QDC to ToT correlation", 
                                               256, 0.5, 4096.5, 128, -100.5, 923.5);
    }

    hpos = new TH2F( "hpos", "Hit position at Veto counters", 10, -0.5, 9.5, 6, -0.5, 5.5);

    f->Write();
    f->Close();
     
    
//    hannes.open("20151003_SDD.txt");


    // here the data is read and stored to the root file
    while( fgets(listline, MAXCHAR, flist) != NULL)//check what name there is in the filelist and store this name under listline
    {
          
        cout << file_name << endl;
        // read name of files from listline and store it under file_name
        sscanf( listline, "%s\n", file_name );

        read_status = ReadData( file_name, rootfile);
      
    }
    
//    fclose(flist);
//    hannes.close();

   // MakePlots( rootfile );

    return;
}



Int_t ReadData( TString datafilename, TString root_file )
{

    Int_t     rtn = -999; 

    Short_t   mul = 0;
    Short_t   dummyword = -1;
    Short_t   tdc_tmp[3] = {0};
    Short_t   rate_buf = 0;
    Short_t   rate_tag = 0;
    Short_t   tdc = 0;

    EventStruct  evt_r; // defined in header


    evt_r.evid = 0,   evt_r.evid1 = 0,  evt_r.evid2 = 0; evt_r.evid3 = 0;
    evt_r.ut    = 0,  evt_r.rate  = -99;   
    evt_r.year  = 2014,  evt_r.month = 0;  
    evt_r.day   = 0,    evt_r.hour = 0,   evt_r.min = 0;
    evt_r.layer = 0,    evt_r.col  = 0;  
    evt_r.trgid  = -1,  evt_r.hr   = 0;

    // Histograms for SiPMs
    TH1F* ht[V1190_TDC_CH];
    TH1F* hvq[V792_QDC_CH];
//    TH1F* hvqc;
//    TH2F* hqt;
    TH2F* hpos;

    // Histograms for SDDs
    TH1F* hsp[NPADC_CH];
    TH1F* hst[7];


    TFile *ff = new TFile( root_file, "UPDATE");
    ff->cd();

    // very strange that using this loop to get the histograms gives segmentation error ...
    // probably due to the limit of memory space ..
    /*
    for( Int_t i=0; i<NPADC_CH; i++ ){
        hsp[i] = (TH1F*)ff->Get( Form("hp[%d]", i) );
        hsp[i]->GetXaxis()->SetTitle("SDD ADC [channel]");
    }
    */
    // Instead only possible to get the histograms from root file like this ...
    hsp[0]  = (TH1F*)ff->Get(Form("hp[%d]", 0));
    hsp[1]  = (TH1F*)ff->Get(Form("hp[%d]", 1));
    hsp[2]  = (TH1F*)ff->Get(Form("hp[%d]", 2));
    hsp[3]  = (TH1F*)ff->Get(Form("hp[%d]", 3));
    hsp[4]  = (TH1F*)ff->Get(Form("hp[%d]", 4));
    hsp[5]  = (TH1F*)ff->Get(Form("hp[%d]", 5));
    hsp[6]  = (TH1F*)ff->Get(Form("hp[%d]", 6));
    hsp[7]  = (TH1F*)ff->Get(Form("hp[%d]", 7));
    hsp[8]  = (TH1F*)ff->Get(Form("hp[%d]", 8));
    hsp[9]  = (TH1F*)ff->Get(Form("hp[%d]", 9));
    hsp[10] = (TH1F*)ff->Get(Form("hp[%d]", 10));
    hsp[11] = (TH1F*)ff->Get(Form("hp[%d]", 11));
    hsp[12] = (TH1F*)ff->Get(Form("hp[%d]", 12));
    hsp[13] = (TH1F*)ff->Get(Form("hp[%d]", 13));
    hsp[14] = (TH1F*)ff->Get(Form("hp[%d]", 14));
    hsp[15] = (TH1F*)ff->Get(Form("hp[%d]", 15));


/*
    for( Int_t i=0; i<7; i++ ){
        Int_t  ihist = i + 48;
        hst[i] = (TH1F*)ff->Get( Form("ht[%d]", ihist) );
        hst[i]->GetXaxis()->SetTitle("SDD timing [ns]");
    }
*/

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
    TBranch *totEvent   = tree->GetBranch("tot");
    TBranch *padcEvent  = tree->GetBranch("padc");
    //TBranch *layerEvent = tree->GetBranch("layer");
    //TBranch *colEvent   = tree->GetBranch("col");
    TBranch *trgidEvent = tree->GetBranch("trgid");
    TBranch *hrEvent    = tree->GetBranch("hr");

    evidEvent  ->SetAddress(&evt_r.evid);
    evid1Event ->SetAddress(&evt_r.evid1);
    evid2Event ->SetAddress(&evt_r.evid2);
    evid3Event ->SetAddress(&evt_r.evid3);
    utEvent    ->SetAddress(&evt_r.ut);
    rateEvent  ->SetAddress(&evt_r.rate);
    mulEvent   ->SetAddress(&evt_r.mul);
    qdcEvent   ->SetAddress(evt_r.qdc);
    tdcEvent   ->SetAddress(evt_r.tdc);
    totEvent   ->SetAddress(evt_r.tot);
    padcEvent  ->SetAddress(evt_r.padc);
    //layerEvent ->SetAddress(&evt_r.layer);
    //colEvent   ->SetAddress(&evt_r.col);
    trgidEvent    ->SetAddress(&evt_r.trgid);
    hrEvent       ->SetAddress(&evt_r.hr);


    FILE *datafile;

    // binary data for Phase One setup in the lab: 
    TString binary_file = BIN_PATH + "/" + datafilename; 

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
        // new line begins with xx 00 00 00 in hex 
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 )
        {
            if ( feof (datafile) ){
                cout << " Reached end of file " << endl;
                break;
            }else {
                cout << " Unexpected termination ! " << endl;
                break;
            }
        }else{
            //cout << " begin of line : " << dummyword << endl;
        }

        evt_r.evid ++;

        // header 0
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        // daq rate 
//        if ( fread(&evt_r.rate, sizeof(evt_r.rate), 1, datafile) != 1 ) break;---------------------------------------------------

         //read another 0 ----  es sind glaub ich jedenfalls 4 wörter vor den daten
     //   if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;//-------------------------------hier evtl ändern wenn qdc anzahl sich ändert

        // set the rate information
        if( evt_r.rate >= rate_buf ){
            if( evt_r.rate > rate_buf ){
                rate_tag = 1; 
            }
            rate_buf = evt_r.rate; 
            evt_r.hr = rate_tag; 
        }else{ 
            rate_tag = -1; 
            evt_r.hr = rate_tag;
            rate_buf = evt_r.rate;
        }

 
        // evid from qdc  
         if ( fread(&evt_r.evid1, sizeof(evt_r.evid1), 1, datafile) != 1 ) break;     // comment out if only PADC data

        // V792_QDC_CH coloumns for qdc
        if ( fread(&evt_r.qdc, sizeof(evt_r.qdc[0]), V792_QDC_CH, datafile) != V792_QDC_CH ) break;
        

        // 0 column for 32 ch qdc and -1 division
	if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break; // -1 coloumn only if disabled
        //if ( fread(&dummyword, sizeof(dummyword), 2, datafile) != 2 ) break;
  //      if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break; // -1 coloumn only if disabled

        for( Int_t i=0; i < V785_ADC_CH; i++){
            evt_r.padc[i] = 0;  
        }


         // evid from padc 
     //   if ( fread(&evt_r.evid2, sizeof(evt_r.evid2), 1, datafile) != 1 ) break;--------------------------------------

        // V785_ADC_CH coloumns for ph adc for Phase One binary data 
        // number of enabled channels NPADC_CH is controlled from the DAQ  
        if ( fread(&evt_r.padc, sizeof(evt_r.padc[0]), ENABLED_NPADC_CH, datafile) != ENABLED_NPADC_CH ) break;
        if ( fread(&dummyword, sizeof(dummyword), 2, datafile) != 2 ) break;  // 0, and -1 division coloumn afterwards;
        //if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break // -1 coloumn only if disabled
       
       // if ( fread(&dummyword, sizeof(dummyword), 4, datafile) != 4 ) break;  // read 3 times 0 and -1
       //                enable this if only PADC data -----------------------------
        for( Int_t i=0; i < V1190_TDC_CH; i++){
            evt_r.tdc[i] = -10,  evt_r.tot[i] = -10; 
        }

    
	// read the tdc event id
	if ( fread(&evt_r.evid3, sizeof(dummyword), 1, datafile) != 1 ) break;



        // read the tdc data, empty channels has an entry of 9999 
        for( Int_t i=0; i < V1190_TDC_CH; i++){
            if( fread ( &tdc_tmp, sizeof(tdc_tmp[0]), 3, datafile ) != 3 ) break;
            
            if( tdc_tmp[0] < -10000){      // after 10.04.2014 the out of range events are the limit of Short_t
                tdc_tmp[0] = -10;
                tdc_tmp[2] = -10;
            }
	    else{
                mul ++ ;
            }

            //   cout << tdc_tmp[0] << " " << tdc_tmp[1] << " " << tdc_tmp[2] << " "  << endl;

            evt_r.tdc[tdc_tmp[1]] = tdc_tmp[0]; 
            evt_r.tot[tdc_tmp[1]] = tdc_tmp[2];
            //cout << "ch : " << tdc_tmp[1] << ";  tot : " << evt_r.tot[tdc_tmp[1]] << endl;
            evt_r.mul = mul;
            

	 }

   

/*            ///////////////////////////////////
            // Histograms for SDDs 
            // Reference TDC channel changed to tdc[61] after the self trigger of SDDs is included
            // from 10.04.2014, however the SDD timing is still given by tdc[48] - tdc[63]; 
            // 
            if( evt_r.tdc[48] > 0 ){     // SDD OR channel
                if( evt_r.tdc[48] >= 1150 && evt_r.tdc[48]<1250 ){
                    hst[0]->Fill( ( evt_r.tdc[48] - evt_r.tdc[63] - 1075 ) * NSECTDCCH );
                }else if( evt_r.tdc[48] > 0 && evt_r.tdc[63] > 0 ){ 
                    hst[0]->Fill( ( evt_r.tdc[48] - evt_r.tdc[63] ) * NSECTDCCH );
                }
            }
            if( evt_r.tdc[49] > 0 ){     // SDD 4
                hst[1]->Fill( ( evt_r.tdc[49] - evt_r.tdc[63] ) * NSECTDCCH );
                hsp[9]->Fill( evt_r.padc[9] );
                hsp[15]->Fill( evt_r.padc[9] );
            }
            if( evt_r.tdc[50] > 0 ){   // SDD 3 
                hst[2]->Fill( ( evt_r.tdc[50] - evt_r.tdc[63] ) * NSECTDCCH );
                hsp[7]->Fill( evt_r.padc[7] );
                hsp[15]->Fill( evt_r.padc[7] );
            }
            if( evt_r.tdc[51] > 0 ){   // SDD 5
                hst[3]->Fill( ( evt_r.tdc[51] - evt_r.tdc[63] ) * NSECTDCCH );
                hsp[11]->Fill( evt_r.padc[11] );
                hsp[15]->Fill( evt_r.padc[11] );
            }
            if( evt_r.tdc[52] > 0 ){  // Open channel at TDC
                //hst[4]->Fill( ( evt_r.tdc[52] - evt_r.tdc[63] ) * NSECTDCCH );
                //hsp[]->Fill( evt_r.padc[] );
            }
            if( evt_r.tdc[53] > 0 ){   // SDD 1
                if( evt_r.tdc[53] >= 1100 && evt_r.tdc[53]<1250 ){
                    hst[5]->Fill( ( evt_r.tdc[53] - evt_r.tdc[63] - 1075 ) * NSECTDCCH );
                }else{
                    hst[5]->Fill( ( evt_r.tdc[53] - evt_r.tdc[63] ) * NSECTDCCH );
                }
                hsp[3]->Fill( evt_r.padc[3] );
                hsp[15]->Fill( evt_r.padc[3] );
            }

            // Fill SDD 6 without TDC data recorded
            if( evt_r.padc[13] > 200 ){
                hsp[12]->Fill( evt_r.padc[13] );
                hsp[13]->Fill( evt_r.padc[13] );
                hsp[14]->Fill( evt_r.padc[13] );
            }
            // Other SDD adc without timing data selection. 
            if( evt_r.padc[3] > 200 ){
                hsp[2]->Fill( evt_r.padc[3] );
                hsp[14]->Fill( evt_r.padc[3] );
            }
            if( evt_r.padc[5] > 200 ){
                hsp[4]->Fill( evt_r.padc[5] );
                hsp[14]->Fill( evt_r.padc[5] );
            }
            if( evt_r.padc[7] > 200 ){
                hsp[6]->Fill( evt_r.padc[7] );
                hsp[14]->Fill( evt_r.padc[7] );
            }
            if( evt_r.padc[9] > 200 ){
                hsp[8]->Fill( evt_r.padc[9] );
                hsp[14]->Fill( evt_r.padc[9] );
            }
            if( evt_r.padc[11] > 200 ){
                hsp[10]->Fill( evt_r.padc[11] );
                hsp[14]->Fill( evt_r.padc[11] );
            }
        }
        /////////////////////////////////
    //         enable this in case of only PADC data -----------------------------
*/

        // Fill the SDD histograms
       for( Int_t j=0; j < ENABLED_NPADC_CH; j++ ){

            hsp[j] -> Fill(evt_r.padc[j]);
	        
	}

  	// Fill in qdc histograms
	for( Int_t j=0; j < NScinti; j++ ){

	    hvq[j] = (TH1F*)ff->Get( Form("hq[%d]", j) );

        }


      for ( Int_t j=0; j < NScinti; j++ ){

           hvq[j] -> Fill( evt_r.qdc[j] );

      }

       // Fill in tdc histograms
       for( Int_t j = 0; j < V1190_TDC_CH; j++ ){	 

            ht[j] = (TH1F*)ff->Get( Form( "ht[%d]", j ) );
                
	}

       for ( Int_t j=0; j < V1190_TDC_CH; j++ ){

           tdc = evt_r.tdc[j] - evt_r.tdc[TDC_REF_CH];

           ht[j] -> Fill( tdc );

      }

         


 /*       //write data to textfile
        for( Int_t j=0; j < 3; j++ ){

            if( j == 2 ) { hannes << evt_r.padc[j] << endl; }
            else  {  hannes << evt_r.padc[j] << "\t"; }

        }
*/
/*
        // Fill the histograms for SiPMs
        Int_t    tdc_ch =  0;
        Float_t  tdc    = .0;

        for( Int_t j=0; j < NScinti; j++ )
        {

            tdc_ch = Qdc2Tdc[j];
            tdc = evt_r.tdc[TDC_REF_CH] - evt_r.tdc[tdc_ch] ;

            // tdc
            if( tdc > -200 && tdc < 2000 ) {
                ht = (TH1F*)ff->Get( Form( "ht[%d]", tdc_ch ) );
                ht->Fill( tdc ); 
            }

            // qdc -> qdc data into histogram
            hvq = (TH1F*)ff->Get( Form("hq[%d]", j) );
            hvq->Fill( evt_r.qdc[j] );

            hvqc = (TH1F*)ff->Get( Form("hqc[%d]", j) );
            if( evt_r.tot[tdc_ch] > 0 ){
                hvqc->Fill( evt_r.qdc[j] );
            }

            // tot qdc correlation
            if( j == 13 || j == 14 ){
                if( evt_r.qdc[j] > 300 ){
                    hqt = (TH2F*)ff->Get( Form( "hqt[%d]", j) );
                    hqt->Fill( evt_r.qdc[j], evt_r.tot[tdc_ch] );
                }
            }else{
                if( evt_r.qdc[j] > 400 ){
                    hqt = (TH2F*)ff->Get( Form( "hqt[%d]", j) );
                    hqt->Fill( evt_r.qdc[j], evt_r.tot[tdc_ch] );
                }
            }

            // Fill hit position
            if( evt_r.tot[tdc_ch] > 0 && evt_r.tot[Qdc2Tdc[j]] > 0 ){
                hpos = (TH2F*)ff->Get( "hpos" );
                evt_r.layer = Qdc2Layer[j];
                evt_r.col   = Qdc2Column[j];
                hpos->Fill( evt_r.col, evt_r.layer );
            }

        }
*/
/*
        ///////////////////////////////
        // Sort the hit pattern of trigger
        // S17 & S22 : 0;  S18 & S24 : 1;  S19 & S25 : 2;  S20 & S26 : 3;  S21 & S27 : 4;
        // S17 & S24 :10;  S18 & S25 :12;  S19 & S26 :14;  S20 & S27 :16;  
        //                 S18 & S22 :13;  S19 & S24 :15;  S20 & S25 :17;  S21 & S26 :19;
        if( mul >= 7 && mul < 10 )
        {
            evt_r.trgid = 5;

            for( Int_t i=0; i<5; i++ ) 
            {
                // straight two layers 
                if( evt_r.tot[Qdc2Tdc[i*2]] > 0 && evt_r.tot[Qdc2Tdc[i*2+1]] > 0
                        &&  evt_r.tot[Qdc2Tdc[i*2+10]] > 0 && evt_r.tot[Qdc2Tdc[i*2+11]] > 0 )
                {
                    evt_r.trgid = i;
                }
                // cross two layers
                if( i < 4 ){
                    if( evt_r.tot[Qdc2Tdc[i*2]] > 0 && evt_r.tot[Qdc2Tdc[i*2+1]] > 0
                            &&  evt_r.tot[Qdc2Tdc[i*2+12]] > 0 && evt_r.tot[Qdc2Tdc[i*2+13]] > 0 )
                    {
                        evt_r.trgid = i*2 + 10;
                    }
                }
                if( i > 0 ){
                    if( evt_r.tot[Qdc2Tdc[i*2]] > 0 && evt_r.tot[Qdc2Tdc[i*2+1]] > 0
                            &&  evt_r.tot[Qdc2Tdc[i*2+8]] > 0 && evt_r.tot[Qdc2Tdc[i*2+9]] > 0 )
                    {
                        evt_r.trgid = i*2 + 11;
                    }
                }
            }
        }else if( mul >= 10 && mul < 13 ){
            evt_r.trgid = 29;
            // straight two layers 
            for( Int_t i=0; i<5; i++ ) 
            {
                if( evt_r.tot[Qdc2Tdc[i*2]] > 0 && evt_r.tot[Qdc2Tdc[i*2+1]] > 0
                        &&  evt_r.tot[Qdc2Tdc[i*2+10]] > 0 && evt_r.tot[Qdc2Tdc[i*2+11]] > 0 )
                {
                    evt_r.trgid = i + 20 ;
                }
            }
        }else if( mul >= 13 && mul < 16 ){
            evt_r.trgid = 39;
            for( Int_t i=0; i<5; i++ ) 
            {
                if( evt_r.tot[Qdc2Tdc[i*2]] > 0 && evt_r.tot[Qdc2Tdc[i*2+1]] > 0
                        &&  evt_r.tot[Qdc2Tdc[i*2+10]] > 0 && evt_r.tot[Qdc2Tdc[i*2+11]] > 0 )
                {
                    evt_r.trgid = i + 30 ;
                }
            }
        }else{
            evt_r.trgid = 49;
            for( Int_t i=0; i<5; i++ ) 
            {
                if( evt_r.tot[Qdc2Tdc[i*2]] > 0 && evt_r.tot[Qdc2Tdc[i*2+1]] > 0
                        &&  evt_r.tot[Qdc2Tdc[i*2+10]] > 0 && evt_r.tot[Qdc2Tdc[i*2+11]] > 0 )
                {
                    evt_r.trgid = i + 40 ;
                }
            }
        }

*/
        // Fill Tree
        tree->Fill();

        mul = 0;

        //  0 column, tdc word count, and -1 division, rate coloumn  after 
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;
        if ( fread(&dummyword, sizeof(dummyword), 1, datafile) != 1 ) break;

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
void MakePlots(TString rootfile )
{
    // Plot the qdc, tdc, tot to qdc for SiPMs 02 Mar. 2014

    TFile *f = new TFile(rootfile, "READ");
    f->cd();

    TH1F* ht;
    TH1F* hsddt;
    TH1F* hvq;
    TH1F* hvqc;
    TH2F* hqt;
    TH2F* hpos;

    TString graph_title;

    Int_t tdc_ch = 0;

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    cc->cd();
    cc->SetGridx();
    cc->SetGridy();

    for( Int_t i=0; i<NScinti; i++ ){

        tdc_ch = Qdc2Tdc[i];

        hqt = (TH2F*)f->Get( Form("hqt[%d]", i) );
        hqt->GetXaxis()->SetTitle("QDC [channel]");
        hqt->GetYaxis()->SetTitle("Time over threshold [ns]");
        hqt->Draw("BOX");
        cc->SetLogy(0);
        cc->Print( Form( GRAPH_PATH + "/ToT_QDC_Correlation_SiPM%02d.pdf", i ) );

        ht = (TH1F*)f->Get( Form("ht[%d]", tdc_ch ) );
        ht->GetXaxis()->SetTitle("TDC [ns]");
        ht->Draw();
        cc->SetLogy(0);
        cc->Print( Form( GRAPH_PATH + "/TDC_SiPM%02d_TDC_CH%02d.pdf", i, tdc_ch ) );

        hvq  = (TH1F*)f->Get( Form("hq[%d]", i ) );
        //FitSpectrum( hvq );
        hvqc = (TH1F*)f->Get( Form("hqc[%d]", i ) );
        hvq->GetXaxis()->SetTitle("QDC [ch]");
        hvq->Draw();
        hvqc->SetLineColor(2);
        hvqc->Draw("same");
        cc->SetLogy();
        cc->Print( Form( GRAPH_PATH + "/QDC_SiPM_CH%02d.pdf", i) );
    }

    hpos= (TH2F*)f->Get( "hpos" );
    hpos->GetXaxis()->SetTitle("Column");
    hpos->GetYaxis()->SetTitle("Layer");
    hpos->Draw("BOX,TEXT");
    cc->SetLogy(0);
    graph_title = GRAPH_PATH + "/Veto_Counter_hit_pos_distribution.pdf";
    cc->Print( graph_title );


    hsddt = (TH1F*)f->Get( Form("ht[48]"));
    hsddt->GetXaxis()->SetTitle("SDD timing [ns]");
    hsddt->Draw();
    graph_title = GRAPH_PATH + "/SDD_Timing.pdf";
    cc->Print( graph_title );


    cc->Close();
    f->Close();

    return; 
}

void MakeOneQDCPlot(TString rootfile)
{
    TString rootfilename = ROOT_PATH + "/" + rootfile;
    TFile *f = new TFile(rootfilename, "READ");
    f->cd();

    TH1F* hvq;
   // Int_t list [4] = {15,16,19,21};

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    cc->Divide(4,8);
    cc->cd();
    
    cc->SetGridx();
    cc->SetGridy();

    for( Int_t i = 0 ; i < 32 ; i++){

       cc->cd(i+1);
       gPad->SetLogy();
       hvq = (TH1F*)f->Get( Form("hq[%d]", i ) );
       hvq->GetXaxis()->SetTitle( Form("QDC [%d]", i) );
       
       hvq->Draw();    
       
    }

}


void MakeOnePADCPlot(TString rootfile)
{
    TString rootfilename = ROOT_PATH + "/" + rootfile;
    TFile *f = new TFile(rootfilename, "READ");
    f->cd();

    TH1F* hsp;
    Int_t list [6] = {0,2,3,5,6,7};
  //  Int_t list [1] = {2};

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    cc->Divide(2,3);
    cc->cd();
    
    cc->SetGridx();
    cc->SetGridy();

    for( Int_t i = 0 ; i < 6 ; i++){

       
       cc->cd(i+1);
       gPad->SetLogy();
       
       hsp = (TH1F*)f->Get( Form("hp[%d]", list[i] ) );
       //hsp -> SetMaximum(5000);
       //hsp -> SetLineColor(2 for red);
       hsp -> GetXaxis()->SetRangeUser(0,1200);
       hsp->GetXaxis()->SetTitle( Form("PADC [%d]" , i));
       
       hsp->Draw();    
       
    }

}


void MakeOnePADCPlot2Files(TString rootfile1, TString rootfile2)
{
   Double_t scaling_factor;
   TH1F* hsp1;
   TH1F* hsp2;

   TString rootfilename1 = ROOT_PATH + "/" + rootfile1;
   TString rootfilename2 = ROOT_PATH + "/" + rootfile2;

   TFile *f1 = new TFile(rootfilename1, "READ");
   TFile *f2 = new TFile(rootfilename2, "READ");

   TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
   cc->cd();

   hsp1 = (TH1F*)f1->Get("hp[7]");
   hsp2 = (TH1F*)f2->Get("hp[7]");

   Double_t counts_one = hsp1->GetEntries();
   Double_t counts_two = hsp2->GetEntries();

   scaling_factor = counts_one/counts_two;

   hsp2->SetLineColor(2);
   hsp2 -> GetXaxis()->SetRangeUser(500,1200);

   hsp2->Scale(scaling_factor);

   hsp1->Fit("gaus","","",670,730);
   hsp1->GetFunction("gaus")->SetLineColor(4);

   hsp2->Fit("gaus","","",670,730);
   hsp2->GetFunction("gaus")->SetLineColor(2);

   hsp2->Draw();
   hsp1->Draw("same");
  
   
}

void GetMnPeakPos(){

   TString rootfilename1 = ROOT_PATH + "/" + "0615_without.root";
  
   Double_t thr = 0.9;
   Double_t NumberOfFiles = 1;

   TFile *f1 = new TFile(rootfilename1, "READ");
   TH1F* htmp;
   htmp = (TH1F*)f1->Get("hp[7]");

   htmp -> GetXaxis()->SetRangeUser(500,900);

   TSpectrum *s = new TSpectrum();
   s->Search(htmp, 1, "", thr );                

   Double_t peakN = s->GetNPeaks();

   Double_t *pX = s->GetPositionX();
   Double_t *pY = s->GetPositionY();

   htmp->Fit("gaus","","",670,730);

   cout << "Number of peaks: " << peakN << endl;
   cout << "Mn peak at: " << pX[0] << " " << pX[1] << endl;

}




void MakeOneTDCPlot(TString rootfile)
{
    TString rootfilename = ROOT_PATH + "/" + rootfile;
    TFile *f = new TFile(rootfilename, "READ");
    f->cd();

    TH1F* htt;
   // Int_t list [6] = {0,2,3,5,6,7};
  

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    cc->Divide(3,5);
    cc->cd();
    
    cc->SetGridx();
    cc->SetGridy();

    for( Int_t i = 0 ; i < 15 ; i++){

       
       cc->cd(i+1);
       gPad->SetLogy();
       
       htt = (TH1F*)f->Get( Form("ht[%d]", i ) );
       //hsp -> SetMaximum(5000);
       htt -> GetXaxis()->SetRangeUser(-15000,5000);
       htt->GetXaxis()->SetTitle( Form("TDC [%d]" , i));
       
       htt->Draw();    
       
    }

}














