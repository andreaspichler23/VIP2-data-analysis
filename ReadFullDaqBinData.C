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
 * 
 * final form: Andreas Pichler
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
#include  "PhaseOneScintiConnectionTable.h"
#include  "PhaseOneSDDConnectionTable.h"
#include  "ReadFullDaqBinData.h" 
//#include  "SDDcalib.h"
//#include  "qdc_ana.C"


using namespace std;



void MakeTree( TString filelist, TString rootfilename, TString place ) 
{

    Int_t  read_status = - 999;
    Int_t root_file_event_counter = 0;

    Int_t  qdc_th[NScinti] = { 700 };
    
    //TObjString bin_file;

    EventStruct  evt;  
    Int_t size = place.Sizeof();
    //cout << size << endl;
    
    FILE  *flist;

    if(size==5){if( (flist = fopen(LIST_PATH_LNGS + "/" +filelist, "r")) == NULL ){
        cout << "Cannot open file: " << filelist << endl;
    }}

    if(size==4){if( (flist = fopen(LIST_PATH_SMI + "/" +filelist, "r")) == NULL ){
        cout << "Cannot open file: " << filelist << endl;
    }}

    cout << filelist << endl;
    //cout << rootfilename << endl;
    

    char  listline[MAXCHAR];
    char  file_name[MAXCHAR];


    TString rootfile;
    if(size==4)rootfile = ROOT_PATH_SMI + "/" + rootfilename;
    if(size==5)rootfile = ROOT_PATH_LNGS + "/" + rootfilename;
    
    //cout << rootfile << endl;
    
    TFile *f = new TFile( rootfile, "RECREATE");
    f->cd();

    TTree *tr = new TTree("tr", "SiPM and SDD data");
    tr->Branch("evid",   &evt.evid,   "Event number/I");
    tr->Branch("evid1",  &evt.evid1,  "Event number qdc/S");
    tr->Branch("evid2",  &evt.evid2,  "Event number padc/S");
    tr->Branch("evid3",  &evt.evid3,  "Event number tdc/S");
    tr->Branch("evnr",   &evt.evnr,   "Event number Rootfile/I");
    tr->Branch("ut",     &evt.ut,     "Unix time tag/I");
    tr->Branch("time",   evt.time,    "Array with date and time[6]/S");
    tr->Branch("rate",   &evt.rate,   "Trigger rate/S");
    //tr->Branch("mul",    &evt.mul,    "Event multiplicity/S");
    tr->Branch("qdc",    evt.qdc,     "QDC Channel[32]/S");
    tr->Branch("tdc",    evt.tdc,     "TDC Channel[16]/S");
    tr->Branch("adc",    evt.padc,    "PeakADC[16]/S");
    //tr->Branch("layer",  &evt.layer,  "Scinti layer/S");
    //tr->Branch("col",    &evt.col,    "Scinti col/S");
    tr->Branch("trgid",  &evt.trgid,  "Trigger type/S");
    //tr->Branch("hr",     &evt.hr,     "High Rate/S");
    tr->Branch("clk",    &evt.clk,    "LV ms clock/I");
    tr->Branch("slow",   evt.slow,    "Slow Parameter[19]/D");
    tr->Branch("tdc_dig",evt.tdc_dig, "TDC Digtal array[9]/S");
    tr->Branch("qdc_dig",evt.qdc_dig, "QDC Digtal array[37]/S");
    tr->Branch("adc_dig",evt.adc_dig, "ADC Digtal array[8]/S");
    //tr->Branch("binFile",&bin_file);

    TH1F* ht[ENABLED_TDC_CH];   // TDC 
    TH1F* hp[ENABLED_PADC_CH];  // PADC
    TH1F* he[ENABLED_PADC_CH];  // Energy histogram
    TH1F* hq[NScinti];          // QDC

    TH2F* hpos;                 // Hit position

    for( Int_t i=0; i < ENABLED_TDC_CH; i++ ){
        ht[i] = new TH1F(Form("ht[%d]", i), "Timing spectra", 2048, -200.5, 1847.5);
    }

    for( Int_t i=0; i < ENABLED_PADC_CH; i++ ){
        hp[i] = new TH1F(Form("hp[%d]", i), "PADC spectra", 2048, -0.5, 4095.5);
    }

    for( Int_t i=0; i < ENABLED_PADC_CH; i++ ){
        he[i] = new TH1F(Form("he[%d]", i), "SDD Energy spectra", 32768, -0.5, 32767.5);
    }

    for( Int_t i=0; i < NSiPM; i++ ){
        hq[i] = new TH1F(Form("hq[%d]", i), "QDC spectra", 512, -0.5, 4095.5);
    }

    hpos = new TH2F( "hpos", "Hit position at Veto counters", 8, -0.5, 7.5, 6, -0.5, 5.5); // 8 columns, 6 rows

    f->Write();
    f->Close();
    

    while( fgets(listline, MAXCHAR, flist) != NULL) // reads MAXCHAR = 512 characters (or until newline in our case) from flist (pointer to filestream from file with names of all binary files) and stores it to listline
    { 

        
        sscanf( listline, "%s\n", file_name ); // stores listline to file_name -> necessary bc listline contains newlines
        cout << file_name << endl;
        
        read_status = ReadData( file_name, rootfile, qdc_th, place, root_file_event_counter);

        cout << "Return value : " << read_status << endl << endl;
    }
    
 //   MakePlots( rootfile, qdc_th );

    fclose(flist);

    return;
}

void MakeTreeLoop( TString filelistFile, TString place){
    
    
    FILE  *flist;
    Int_t size = place.Sizeof();

    if( (flist = fopen(LIST_PATH + "/" +filelistFile, "r")) == NULL ){
        cout << "Cannot open file: " << filelistFile << endl;
    }

    cout << filelistFile << endl; // filelistFile is the name of the file which contains a list of strings, which correspond with a ".list" ending to the filenames of the files which include
   //a list of binary files (already existing), and with a ".root" ending to the root files which are to be made

    char  listline[MAXCHAR];
    char  file_name[MAXCHAR];
    //string listline;
    
    TString rootfilename;
    TString filelist;
    TString filename_tstring;

    
    while( fgets(listline, MAXCHAR, flist) != NULL){ // reads MAXCHAR = 512 characters (or until newline in our case) from flist (pointer to filestream from file with names of all file lists) and stores it to listline
        // listline is now a one string in the filelistFile with a newline!
     
        
        
        sscanf( listline, "%s\n", file_name ); // stores listline to file_name -> necessary bc listline contains newlines
        cout << file_name << endl;
        
        filename_tstring = TString(file_name);
        
        rootfilename = filename_tstring + ".root";
        filelist     = filename_tstring + ".list";
        
        //cout << rootfilename << " " << filelist << endl;
        
        MakeTree( filelist, rootfilename, place);

        //cout << "Return value : " << read_status << endl << endl;
        cout << endl << endl << endl;
    }
    
    fclose(flist);
    return;
    
}

Int_t ReadData( TString datafilename, TString root_file, Int_t qdc_th[NScinti], TString place, int& root_file_event_counter )
{

    
    long ms_start, ms_current, ms_diff, ms_previous, ms_diffToPrevious;
    Int_t ms_gap = 0;
    Int_t sec_start, sec_diff;
    Int_t success = 0;
    Short_t ref_channel_val = 0;
    Short_t layer, col;
    
    Int_t bin_file_event_count = 0;
    
    Int_t     rtn = -999; 
    
    ofstream logfile;   
    TString logfileName = ANALYSIS_PATH + "/logfile.txt";
    
    logfile.open(logfileName,std::ofstream::app);

    Short_t   mul = 0;
    Short_t   dummyword = -1;
    Short_t   tdc_tmp[3] = {0};

    UShort_t   clk_db = 0,  clk_ub = 0;

    EventStruct  evt_r; 

    evt_r.evid  = 0,     evt_r.evid1 = 0,  evt_r.evid2 = 0,  evt_r.evid3 = 0; evt_r.evnr = 0;
    evt_r.ut    = 0,     evt_r.rate  = -99;   
    evt_r.year  = 2014,  evt_r.month = 0;  
    evt_r.day   = 0,     evt_r.hour  = 0,   evt_r.min = 0; evt_r.sec = 0;
    //evt_r.layer = 0,     evt_r.col   = 0;  
    //evt_r.trgid = -1,    evt_r.hr    = 0;
    evt_r.clk   = 0;
    //const char *binFilename = datafilename.Data();
    //evt_r.binFile = binFilename;
    
    //evt_r.binFile = TText(0,0,datafilename.Data());

    //TObjString bin_file;
      

 
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
    TBranch *evnrEvent  = tree->GetBranch("evnr");
    TBranch *timeEvent  = tree->GetBranch("time");
    TBranch *utEvent    = tree->GetBranch("ut"); // time of the current binary file
    TBranch *rateEvent  = tree->GetBranch("rate");
    //TBranch *mulEvent   = tree->GetBranch("mul");
    TBranch *qdcEvent   = tree->GetBranch("qdc");
    TBranch *tdcEvent   = tree->GetBranch("tdc");
    TBranch *padcEvent  = tree->GetBranch("adc");
    TBranch *slowEvent  = tree->GetBranch("slow");
    //TBranch *layerEvent = tree->GetBranch("layer");
    //TBranch *colEvent   = tree->GetBranch("col");
    TBranch *trgidEvent = tree->GetBranch("trgid");
    //TBranch *hrEvent    = tree->GetBranch("hr");
    TBranch *clkEvent   = tree->GetBranch("clk");
    TBranch *tdc_digEvent = tree->GetBranch("tdc_dig");
    TBranch *qdc_digEvent = tree->GetBranch("qdc_dig");
    TBranch *adc_digEvent = tree->GetBranch("adc_dig");
    //TBranch *binFileEvent = tree->GetBranch("binFile");
  
    
    
    evidEvent  ->SetAddress(&evt_r.evid);
    evid1Event ->SetAddress(&evt_r.evid1);
    evid2Event ->SetAddress(&evt_r.evid2);
    evid3Event ->SetAddress(&evt_r.evid3);
    evnrEvent  ->SetAddress(&evt_r.evnr); 
    timeEvent  ->SetAddress(evt_r.time);
    utEvent    ->SetAddress(&evt_r.ut);
    rateEvent  ->SetAddress(&evt_r.rate);
    //mulEvent   ->SetAddress(&evt_r.mul);
    qdcEvent   ->SetAddress(evt_r.qdc);
    tdcEvent   ->SetAddress(evt_r.tdc);
    padcEvent  ->SetAddress(evt_r.padc);
    slowEvent  ->SetAddress(evt_r.slow); 
    //layerEvent ->SetAddress(&evt_r.layer);
    //colEvent   ->SetAddress(&evt_r.col);
    trgidEvent ->SetAddress(&evt_r.trgid);
    //hrEvent    ->SetAddress(&evt_r.hr);
    clkEvent   ->SetAddress(&evt_r.clk);
    tdc_digEvent  ->SetAddress(evt_r.tdc_dig);
    qdc_digEvent  ->SetAddress(evt_r.qdc_dig);
    adc_digEvent  ->SetAddress(evt_r.adc_dig);
    //binFileEvent ->SetAddress(&bin_file);

    

    //bin_file = TObjString(datafilename);
    //TString test = bin_file.GetString();
    
    //cout << "AATTENTION: " << test << endl;
   
    FILE *datafile;

    // binary data for Phaes One setup in the lab: 
  
    TString binary_file;
    Int_t size = place.Sizeof();
    if(size==5){binary_file = BIN_PATH_LNGS + "/" + datafilename;}
    if(size==4){binary_file = BIN_PATH_SMI + "/" + datafilename; }

    if((datafile = fopen(binary_file, "r"))==NULL){
        cout << "Cannot open file: " << binary_file << endl;
        logfile << "Cannot open file: " << binary_file << endl;
    }

    cout << binary_file << endl;

    // Get the time tag of the run
    sscanf ( datafilename, "%04hd%02hd%02hd_%02hd%02hd", &evt_r.year, &evt_r.month, // ?!?!?!??!??!!?!?!? if longer filename!! -> seems to be fine!!
                           &evt_r.day, &evt_r.hour, &evt_r.min ); // %hd specifies short integer

    TDatime *datime_start = new TDatime( evt_r.year, evt_r.month, evt_r.day, 
                                   evt_r.hour, evt_r.min,   0 );
    evt_r.ut   = datime_start->Convert();
    sec_start = evt_r.ut;
    //datime_start->Print();
    //cout << "year: " << evt_r.year << ",month: " << evt_r.month << ",day: " << evt_r.day << ",hour: " << evt_r.hour << ",minute: " << evt_r.min << endl;
    
    //evt_r.time = evt_r.year*1e10 + evt_r.month*1e8 + evt_r.day*1e6 
    //             + evt_r.hour*1e4 + evt_r.min*1e2; 
    //delete datime;
    
    
    // Read the binary data
    while(1)  
    //for(Int_t bla = 1; bla < 10; bla ++)
    {
       
        //cout << root_file_event_counter << " " << bin_file_event_count  << " " << evt_r.evnr << endl;
        
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
        
        root_file_event_counter += 1;
        bin_file_event_count += 1;
        evt_r.evnr = root_file_event_counter;
        
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

        ///          enable this comment out if only PADC data
        for( Int_t i=0; i < ENABLED_TDC_CH; i++){
            evt_r.tdc[i] = -10;
        }

        // read the tdc data, empty channel has an entry of -31073 
        for( Int_t i=0; i < ENABLED_TDC_CH ; i++){
            if( fread ( &tdc_tmp[0], sizeof(tdc_tmp[0]), 1, datafile ) != 1 ) break;
            if( fread ( &tdc_tmp[1], sizeof(tdc_tmp[1]), 1, datafile ) != 1 ) break;
            if( fread ( &tdc_tmp[2], sizeof(tdc_tmp[2]), 1, datafile ) != 1 ) break;
            if( tdc_tmp[0] < -30000){      // after 10.04.2014 the out of range events are the limit of Short_t
                tdc_tmp[0] = -100; // no data set to -100 in rootfile
                tdc_tmp[2] = -10;  // does not make any sense; would be 0 otherwise; is not used afterwards
            }

            if( tdc_tmp[1] == i ){
                evt_r.tdc[tdc_tmp[1]] = tdc_tmp[0]; 
            }
            
        }
        //evt_r.mul = mul;
        // Fill the digital branches
        
        //TDC; only using bottom and top
          
        ref_channel_val = evt_r.tdc[14];
        evt_r.tdc_dig[8] = 0; // tdc multiplicity
        
        if( evt_r.tdc[4] > -200 && evt_r.tdc[4] < 0){ evt_r.tdc_dig[0] = 0; } // top outer tdc
        else if( evt_r.tdc[4] > 0 && evt_r.tdc[4] < 7001 ){ evt_r.tdc_dig[0] = 1; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else if( evt_r.tdc[4] > 7000 && evt_r.tdc[4] < 16000){ evt_r.tdc_dig[0] = 2; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else{ evt_r.tdc_dig[0] = -1; }

        if( evt_r.tdc[5] > -200 && evt_r.tdc[5] < 0){ evt_r.tdc_dig[1] = 0; } // top inner tdc
        else if( evt_r.tdc[5] > 0 && evt_r.tdc[5] < 7001 ){ evt_r.tdc_dig[1] = 1; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else if( evt_r.tdc[5] > 7000 && evt_r.tdc[5] < 16000){ evt_r.tdc_dig[1] = 2; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else{ evt_r.tdc_dig[1] = -1; }
        
        if( evt_r.tdc[1] > -200 && evt_r.tdc[1] < 0){ evt_r.tdc_dig[2] = 0; } // bottom inner tdc
        else if( evt_r.tdc[1] > 0 && evt_r.tdc[1] < 7001 ){ evt_r.tdc_dig[2] = 1; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else if( evt_r.tdc[1] > 7000 && evt_r.tdc[1] < 16000){ evt_r.tdc_dig[2] = 2; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else{ evt_r.tdc_dig[2] = -1; }
        
        if( evt_r.tdc[8] > -200 && evt_r.tdc[8] < 0){ evt_r.tdc_dig[3] = 0; } // bottom outer tdc
        else if( evt_r.tdc[8] > 0 && evt_r.tdc[8] < 7001 ){ evt_r.tdc_dig[3] = 1; evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1;}
        else if( evt_r.tdc[8] > 7000 && evt_r.tdc[8] < 16000){ evt_r.tdc_dig[3] = 2;evt_r.tdc_dig[8] = evt_r.tdc_dig[8] + 1; }
        else{ evt_r.tdc_dig[3] = -1; }   
        
        if( evt_r.tdc_dig[0] > 0 || evt_r.tdc_dig[3] > 0){ evt_r.tdc_dig[4] = 1; } // tdc outer flag
        else{ evt_r.tdc_dig[4] = 0; }
       
        if( evt_r.tdc_dig[1] > 0 || evt_r.tdc_dig[2] > 0){ evt_r.tdc_dig[5] = 1; } // tdc inner flag
        else{ evt_r.tdc_dig[5] = 0; }
        
        if( evt_r.tdc_dig[4] > 0 || evt_r.tdc_dig[5] > 0 ){ evt_r.tdc_dig[6] = 1; } // outer OR inner
        else( evt_r.tdc_dig[6] = 0 );
        
        if( evt_r.tdc_dig[4] > 0 && evt_r.tdc_dig[5] > 0 ){ evt_r.tdc_dig[7] = 1; } // outer AND inner
        else( evt_r.tdc_dig[7] = 0 );
        
        // QDC 
        evt_r.qdc_dig[32] = 0; // qdc outer flag
        evt_r.qdc_dig[33] = 0; // qdc inner flag
        
        evt_r.qdc_dig[36] = 0; // qdc multiplicity
        hpos = (TH2F*)ff->Get( "hpos" );
        
        for( Int_t i = 0 ; i < 24 ; i++ ){
            //if( i == 3 ){ continue; }
            if( evt_r.qdc[i] > QdcCuts_100nsGate[i] ){ 
                
                evt_r.qdc_dig[i] = 1;
                layer = Qdc2Layer[i];
                col = Qdc2Column[i];
                evt_r.qdc_dig[36] = evt_r.qdc_dig[36] + 1; // qdc multiplicity
                
                hpos -> Fill( col, layer );
                
                if( Qdc2Outer[i] == 1 ){ evt_r.qdc_dig[32] = 1; } // qdc outer flag
                else{ evt_r.qdc_dig[33] = 1; } // qdc inner flag
                           
            }
            else{ evt_r.qdc_dig[i] = 0; }
            
            
        }
        
        if( evt_r.qdc_dig[32] == 1 || evt_r.qdc_dig[33] == 1 ){ evt_r.qdc_dig[34] = 1; } // inner OR outer
        else{ evt_r.qdc_dig[34] = 0; }

        if( evt_r.qdc_dig[32] == 1 && evt_r.qdc_dig[33] == 1 ){ evt_r.qdc_dig[35] = 1; } // inner AND outer
        else{ evt_r.qdc_dig[35] = 0; }  
        
        //ADC
        
        evt_r.adc_dig[6] = 0; // adc OR
        evt_r.adc_dig[7] = 0; // adc multiplicity
        for( Int_t sdd = 1 ; sdd < 7 ; sdd++ ){
            
            if( evt_r.padc[SDDToPadc[sdd]] > AdcCuts[sdd-1] ){
                evt_r.adc_dig[sdd-1] = 1; 
                evt_r.adc_dig[6] = 1; // adc OR
                evt_r.adc_dig[7] = evt_r.adc_dig[7] + 1; // adc multiplicity
            }
            else{ evt_r.adc_dig[sdd-1] = 0; }
            
            
        }
        
        // TRIGID
        
        evt_r.trgid = -1;
        
        if( evt_r.tdc[12] > 12000 && evt_r.tdc_dig[6] == 0 && evt_r.tdc[13] < 0 ){ evt_r.trgid = 1; } // SDD only
        else if( evt_r.tdc[12] < 0 && (evt_r.tdc_dig[6] == 1 || evt_r.tdc[13] > 0) ){ evt_r.trgid = 2; } // scintillator only
        else if( evt_r.tdc[12] > 12000 && (evt_r.tdc_dig[7] == 1 || evt_r.tdc[13] > 0 ) ){ evt_r.trgid = 3; } // sdd + inner AND outer tdc layer
        else if( evt_r.tdc[12] > 12000 && evt_r.tdc_dig[6] == 1 ){ evt_r.trgid = 4; } // sdd + 1 scintillator layer
        else if( evt_r.tdc[12] < 0 && evt_r.tdc[13] < 0 ){ evt_r.trgid = 0; } // no trigger

        // Histograms for TDC 
        // Reference TDC channel tdc[14],  SDD OR ch[12],  scinti trigger ch[13].
        // the SDD timing is given by tdc[12] - tdc[14]; 

        Float_t  tdc    = .0;
        for( Int_t i=0; i < ENABLED_TDC_CH ; i++){
            ht = (TH1F*)ff->Get( Form( "ht[%d]", i) );
            if( evt_r.tdc[i] > 0 ) {
                tdc = ( evt_r.tdc[TDC_REF_CH] - evt_r.tdc[i] ) * NSECTDCCH ;
                //cout << "tdc ch : " << i << ";  tdc value : "  << tdc << "; " << endl;
                ht->Fill( tdc ); 
            }
        }
/*

        for( Int_t i=0; i < ENABLED_TDC_CH ; i++){// -------------------------------- change here again!!!!!!!!!!!!!!!!!!!!!
            ht = (TH1F*)ff->Get( Form( "ht[%d]", i) );
            if( evt_r.tdc[i] > 0 ) {
                tdc = ( evt_r.tdc[i] ) * NSECTDCCH ;
                //cout << "tdc ch : " << i << ";  tdc value : "  << tdc << "; " << endl;
                ht->Fill( tdc ); 
            }
        }

 */      

        // Fill the histograms for SiPMs
        // qdc

        
        
        for( Int_t j=0; j < NSiPM; j++ )
        {
            hvq = (TH1F*)ff->Get( Form("hq[%d]", j) );
            hvq->Fill( evt_r.qdc[j] );
	
            // Fill hit position
            // Threshold value 
           // if( qdc_th[j] == 0 ){
           //     qdc_th[j] = 700;
           // }

         
//            if( evt_r.qdc[j] > QdcCuts[j] ){// -------------------------------------------this has to be changed again!!!!!!!!!!
//                evt_r.layer = Qdc2Layer[j];
//                evt_r.col   = Qdc2Column[j];
//                hpos->Fill( evt_r.col, evt_r.layer );
//            }
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
        if ( fread(&clk_ub, sizeof(clk_ub), 1, datafile) != 1 ) break;
        if ( fread(&clk_db, sizeof(clk_db), 1, datafile) != 1 ) break;

        evt_r.clk = clk_ub << 16 ^ clk_db; // = shift clk_ub 16 places to the left -> towards the msb; ^ is the exclusive or: = 1 if either one bit is 1; but 0 if both are 0 or both are 1
        // -> as the lower 16 bits are 0 after the shift, this results in: [16 bits clk_ub] [16 bits clk_db]
        
        if(bin_file_event_count == 1){ ms_start = evt_r.clk;}
        
        ms_current = evt_r.clk;
        
        if( ms_gap == 1 ){ ms_current = ms_current + 4294967295; }
        
        ms_diff = ms_current - ms_start;
        ms_diffToPrevious = ms_current - ms_previous;
        
        if( ms_diff < 0 ){ // sometimes the ms_current is at -2036400115 for single events
            
            cout << "ms_diff: " << ms_diff << " ms diff 2 previous: " << ms_diffToPrevious <<
                    " ms gap = " << ms_gap << " ms current: " << ms_current << " ms_previous: " << ms_previous << " ms start: " << ms_start << endl;
            if( ms_current > -2147480000 ){ ms_current = ms_previous; ms_diff = ms_current - ms_start; } 
            else{ cout << "THERE WAS A RESET OF THE MS CLOCK: ms current: " << ms_current << endl; ms_gap = 1; }
        }

        
        
        sec_diff = (Int_t)ms_diff/1000;
        evt_r.ut = sec_start + sec_diff;
        
        TDatime datime_curr(evt_r.ut);
        
        evt_r.time[5]  = datime_curr.GetSecond();
        evt_r.time[4]  = datime_curr.GetMinute();
        evt_r.time[3] =  datime_curr.GetHour();
        evt_r.time[2]  = datime_curr.GetDay();
        evt_r.time[1]  = datime_curr.GetMonth();
        evt_r.time[0]  = datime_curr.GetYear();
        

//       cout << "slow time before: " << evt_r.slow[0] << " " << evt_r.slow[1] << " " << evt_r.slow[2] << " " << evt_r.slow[3] << " " << evt_r.slow[4] << " " << evt_r.slow[5] << " " << evt_r.slow[6] << " " 
//                << evt_r.slow[7] << " " << evt_r.slow[8] << " " << evt_r.slow[9] << " " << evt_r.slow[10] << " " << evt_r.slow[11] << " " << evt_r.slow[12] << " " << evt_r.slow[13] << " " 
//                << evt_r.slow[14] << " " << evt_r.slow[15] << " " << evt_r.slow[16] << " " << evt_r.slow[17] << " " << evt_r.slow[18] << endl;
        success = GetSlowParameter(evt_r.ut,evt_r.slow);

        if(success == 0){
            
            logfile << "PROBLEM!!!, file: " << datafilename << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << " day: " << evt_r.time[2] << " month: " << evt_r.time[1] << " year: " << evt_r.time[0] << endl;
            logfile << "no slow control data found: stopping run of current binary file: " << binary_file << endl << endl;
            break;
            
        }
        ms_previous = ms_current;
        
        //cout << ms_start << " " << ms_diff << " " << ms_current <<  " " << sec_start << " " << sec_diff << " "<<  evt_r.ut << endl;
        
        
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
    
    logfile.close();
    delete ff;
    delete datime_start;

    cout << "Finished run : " << binary_file << endl;

    return rtn; 
}

Int_t GetSlowParameter(UInt_t time, Double_t *slow){
    
    ifstream slow_file;
    string slow_line;
    
    
    
    Int_t success = 0;
    Int_t loop_count = 0;
    
    for(Int_t i = 0; i < 19; i++){
        
        slow[i] = -1000;
        
    }
    
    start: // start for trying the loop again for a time 60 seconds later
    
    TDatime *datime = new TDatime(time);
    TDatime datime_s;
    
    UInt_t year_s, month_s, day_s, hour_s, min_s, sec_s;
    UInt_t time_s;
    
    UInt_t year = datime->GetYear();
    UInt_t month = datime->GetMonth();
    UInt_t day = datime->GetDay();
    
    UInt_t date = datime->GetDate();
    
    UInt_t hour = datime->GetHour();
    UInt_t min = datime->GetMinute();
    UInt_t sec = datime->GetSecond();
    
    TString date_string = Form("/%d.log",date);
    TString slow_filename = SLOW_TEXT_PATH + date_string;
    
    //cout << slow_filename << endl;
    
    TString comp_string;// comp_string_cut;
    TString slow_line_tstring, slow_line_final;
    TString am_string = Form("AM");
    TString pm_string = Form("PM");
    char dummie[2];
    
    if(min < 10){
        if(hour == 0){comp_string = Form("%d/%d/%d %d:0%d",month,day,year,12,min);}// hour in the format:00 - 23 for TDatime objects
        if(hour > 0 && hour < 12){comp_string = Form("%d/%d/%d %d:0%d",month,day,year,hour,min);}
        if(hour == 12){comp_string = Form("%d/%d/%d %d:0%d",month,day,year,12,min);}
        if(hour > 12){comp_string = Form("%d/%d/%d %d:0%d",month,day,year,hour-12,min);}
    } else{
        if(hour == 0){comp_string = Form("%d/%d/%d %d:%d",month,day,year,12,min);}// hour in the format:00 - 23 for TDatime objects
        if(hour > 0 && hour < 12){comp_string = Form("%d/%d/%d %d:%d",month,day,year,hour,min);}
        if(hour == 12){comp_string = Form("%d/%d/%d %d:%d",month,day,year,12,min);}
        if(hour > 12){comp_string = Form("%d/%d/%d %d:%d",month,day,year,hour-12,min);}        
    }
    // 5/2/2017 12:00:26 AM   
    
    //cout << comp_string << endl;
    
    slow_file.open(slow_filename);
    
    if( slow_file.fail() ){ slow_file.close(); goto end; }
    
    while( getline(slow_file,slow_line) ){

        slow_line_tstring = slow_line;
        //cout << slow_line_tstring << endl;
    
        if(hour < 12){ 
            
            if(slow_line_tstring.Contains(comp_string) && slow_line_tstring.Contains(am_string)){
            
                slow_line_final = slow_line_tstring;
                success = 1;
                
                //cout << slow_line_final;
                break;
            
            }
            
        } else{
        
            if(slow_line_tstring.Contains(comp_string) && slow_line_tstring.Contains(pm_string)){
            
                slow_line_final = slow_line_tstring;
                success = 1;
                
                //cout << slow_line_final;
                break;
            
            }
        
        }
    }  
    
    //cout << "success: " << success << " loop count: " << loop_count << endl; 
    
    if( success == 0){
        
        loop_count += 1;
        time += TMath::Power(-1,loop_count+1) * loop_count * 60;
        slow_file.close();
        
        if( loop_count == 5 ) {cout << "comp string = " << comp_string << endl; goto end;}
        goto start;

    }   
       
    slow_file.close();
    
    sscanf(slow_line_final, "%d/%d/%d %d:%d:%d",&month_s,&day_s,&year_s,&hour_s,&min_s,&sec_s);   // getting the time of the slow file entry and writing it to slow array
    if( hour == 0 && hour_s == 12 ){ hour_s =  0; }
    if( hour > 12 ){ hour_s = hour_s + 12; }
    
    datime_s = TDatime(year_s, month_s, day_s, hour_s, min_s, sec_s);
    time_s = datime_s.Convert();    
    slow[0] = time_s;

    sscanf(slow_line_final, "%d/%d/%d %d:%d:%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",&month_s,&day_s,&year_s,&hour_s,&min_s,&sec_s, dummie, &slow[1], &slow[2],
             &slow[3], &slow[4], &slow[5], &slow[6], &slow[7], &slow[8], &slow[9], &slow[10], &slow[11], &slow[12], &slow[13], &slow[14], &slow[15], &slow[16], &slow[17], &slow[18]);
    
    end:
             
    delete datime;
    return success;
    
}
////////////////////////////////////////////
// Can be directly called if root file exists
//void MakePlots(TString rootfile, Int_t qdc_th[NScinti] )
//{
//    // Plot the qdc, tdc, padc 
//
//    TFile *f = new TFile(rootfile, "READ");
//    f->cd();
//
//    TH1F* ht;
//    TH1F* hmul; 
//    TH1F* hadc;
//    //TH1F* hev;
//    TH1F* hvq;
//    TH2F* hpos;
//
//    TString graph_title;
//
//    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
//    cc->cd();
//    cc->SetGridx();
//    cc->SetGridy();
//
//    for( Int_t i=0; i<ENABLED_TDC_CH; i++){
//        ht = (TH1F*)f->Get( Form("ht[%d]", i) );
//        ht->GetXaxis()->SetTitle("TDC [ns]");
//        ht->Draw();
//        cc->SetLogy();
//        cc->Print( Form( GRAPH_PATH + "/TDC_CH%02d.pdf", i ) );
//    }
//
//    for( Int_t i=0; i<ENABLED_PADC_CH; i++){
//        hadc = (TH1F*)f->Get( Form("hp[%d]", i) );
//        hadc->GetXaxis()->SetTitle("PADC [ch]");
//        hadc->GetXaxis()->SetRangeUser(200, 2500);
//        hadc->Draw();
//        cc->SetLogy(0);
//        cc->Print( Form( GRAPH_PATH + "/PADC_CH%02d.pdf", i ) );
//    }
//
//    for( Int_t i=0; i<NSiPM; i++ ){
//        hvq  = (TH1F*)f->Get( Form("hq[%d]", i ) );
//        hvq->GetXaxis()->SetTitle("QDC [ch]");
//        //hvq->Draw();
//        //cc->SetLogy();
//       
//       // qdc_th[i] = FitQdcSpectrum(hvq, cc);
//        /*
//        if( i < 10 || i == 11 || i >23 ){     // For 201505 data set
//            qdc_th[i] = FitQdcSpectrum(hvq, cc);
//        }else{
//            cc->cd();
//            hvq->Draw();
//            cc->SetLogy();
//        }
//        */
//        cc->Print( Form( GRAPH_PATH + "/QDC_SiPM_CH%02d.pdf", i) );
//    }
//
//
//    // Fill and plot again the hit position 2d plot
//    EventStruct  evt_p; 
//    evt_p.mul = 0;
//
//    TTree *tree = (TTree*)(f->Get("tr"));
//
//    TBranch *evidEvent  = tree->GetBranch("evid");
//    TBranch *qdcEvent   = tree->GetBranch("qdc");
//    TBranch *adcEvent   = tree->GetBranch("adc");
//    TBranch *tdcEvent   = tree->GetBranch("tdc");
//    TBranch *layerEvent = tree->GetBranch("layer");
//    TBranch *colEvent   = tree->GetBranch("col");
//
//    evidEvent  ->SetAddress(&evt_p.evid);
//    qdcEvent   ->SetAddress(evt_p.qdc);
//    adcEvent   ->SetAddress(evt_p.padc);
//    tdcEvent   ->SetAddress(evt_p.tdc);
//    layerEvent ->SetAddress(&evt_p.layer);
//    colEvent   ->SetAddress(&evt_p.col);
//    
//
//    // Reset the bin content of the hit position histogram;
//    Int_t entries = tree->GetEntries();
//    hpos= (TH2F*)f->Get( "hpos" );
//
//    for( Int_t icol = 0; icol < 8; icol ++ ){
//        for( Int_t ilayer = 0; ilayer < 6; ilayer ++ ){
//            hpos->SetBinContent( icol, ilayer, 0.);
//        }
//    }
//
//    hmul = new TH1F("hmul", "scinti multiplicity", 100, -0.5, 99.5);
//
//    // Refill the hit postiion histogram
//    for( Int_t ient=0; ient < entries; ient++ ){
//        for( Int_t j=0; j < NSiPM; j++ )
//        {
//            qdcEvent    ->GetEntry( ient );
//            layerEvent  ->GetEntry( ient );
//            colEvent    ->GetEntry( ient );
//
//            // Fill hit position
//            // Threshold value 
//            if( qdc_th[j] == 0 ){
//                qdc_th[j] = 1000;
//            }
//            if( evt_p.qdc[j] > qdc_th[j] ){
//                evt_p.layer = Qdc2Layer[j];
//                evt_p.col   = Qdc2Column[j];
//                hpos->Fill( evt_p.col, evt_p.layer );
//
//                evt_p.mul ++ ;
//            }
//        }
//        hmul->Fill( evt_p.mul);
//        evt_p.mul = 0;
//    }
//
//
//    /// event mismatch of plus minus one event id check
//    TH1F* hevt  = new TH1F("hevt", "adc all", 512, 0.5, 4096.5);
//    TH1F* hac   = new TH1F("hac", "adc adc cut", 512, 0.5, 4096.5);
//    TH1F* htc   = new TH1F("htc", "adc tdc cut", 512, 0.5, 4096.5);
//    TH1F* hatc  = new TH1F("hatc", "adc tdc mp one event cut", 512, 0.5, 4096.5);
//
//    Bool_t    adc_flag  = 0;
//    Bool_t    adc_flag1 = 0;
//    Bool_t    adc_flag2 = 0;
//    Bool_t    tdc_flag  = 0;
//    Bool_t    tdc_flag1 = 0;
//    Bool_t    tdc_flag2 = 0;
//
//    for( Int_t ient=0; ient < entries; ient++ ){
//        adcEvent    ->GetEntry( ient );
//        tdcEvent    ->GetEntry( ient );
//
//        for( Int_t k=0; k < 1; k++ )
//        {
//
//            if( evt_p.padc[k] > 300 ){
//
//                adc_flag = 1;
//
//                hevt->Fill( (Double_t)(evt_p.padc[k]) );
//                hac->Fill( (Double_t)(evt_p.padc[k]) );
//
//                if( evt_p.tdc[12] > 0 ){
//                    htc->Fill( (Double_t)(evt_p.padc[k]) );
//                    hatc->Fill( (Double_t)(evt_p.padc[k]) );
//
//                    tdc_flag = 1;
//                }
//            }
//        }
//
//        if( adc_flag == 1 && tdc_flag == 0 && ient > 0 && ient < entries )
//        {
//
//            // check previous event 
//            adcEvent    ->GetEntry( ient - 1 );
//            tdcEvent    ->GetEntry( ient - 1 );
//
//            if( evt_p.tdc[12] > 0 )tdc_flag1 = 1; 
//
//            for( Int_t k=0; k < 1; k++ )
//            {
//                if( evt_p.padc[k] > 300 ){
//                    adc_flag1 = 1;
//                }
//            }
//
//            // check following event 
//            adcEvent    ->GetEntry( ient + 1 );
//            tdcEvent    ->GetEntry( ient + 1 );
//
//            if( evt_p.tdc[12] > 0 )tdc_flag2 = 1; 
//
//            for( Int_t k=0; k < 1; k++ )
//            {
//                if( evt_p.padc[k] > 300 ){
//                    adc_flag2 = 1;
//                }
//            }
//        }
//
//        // Fill the hatc histo if pre or next event has no adc entry but tdc entry 
//        if( (adc_flag1 == 0 && tdc_flag1 == 1 && adc_flag == 1) || 
//                (adc_flag2 == 0 && tdc_flag2 == 1 && adc_flag == 1) ){
//
//            adcEvent    ->GetEntry( ient );
//            for( Int_t k=0; k < 1; k++ )
//            {
//                if( evt_p.padc[k] > 300 ){
//                    hatc->Fill( (Double_t)(evt_p.padc[k]) );
//                }
//            }
//        }
//        
//        adc_flag  = 0;
//        adc_flag1 = 0;
//        adc_flag2 = 0;
//        tdc_flag  = 0;
//        tdc_flag1 = 0;
//        tdc_flag2 = 0;
//    }
//
//    cc->cd();
//    hpos->GetXaxis()->SetTitle("Column");
//    hpos->GetYaxis()->SetTitle("Layer");
//    hpos->Draw("BOX,TEXT");
//    cc->SetLogy(0);
//    graph_title = GRAPH_PATH + "/Veto_Counter_hit_pos_distribution.pdf";
//    cc->Print( graph_title );
//
//    cc->cd();
//    hmul->GetXaxis()->SetTitle("scinti multiplicity");
//    hmul->GetXaxis()->SetRangeUser(0, 32);
//    hmul->Draw();
//    cc->SetLogy();
//    graph_title = GRAPH_PATH + "/Veto_Counter_multiplicity.pdf";
//    cc->Print( graph_title );
//
//    cc->cd();
//    cc->SetLogy(0);
//    hevt->GetXaxis()->SetRangeUser(200, 2000);
//    hevt->Draw();
//    hac->Draw("same");
//    htc->SetLineColor(2);
//    htc->Draw("same");
//    graph_title = GRAPH_PATH + "/Evid_matching.pdf";
//    cc->Print( graph_title );
//    hatc->GetXaxis()->SetRangeUser(200, 2000);
//    hatc->Draw();
//    hatc->SetLineColor(2);
//    hevt->Draw("same");
//    graph_title = GRAPH_PATH + "/Evid_matching_pm_events.pdf";
//    cc->Print( graph_title );
//
//    cc->Close();
//
//    f->Close();
//
//    return; 
//}

