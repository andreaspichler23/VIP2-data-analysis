/*------------------------
/
/ written by Andreas Pichler
/
/       in a hurry
/
------------------------*/

#include  <stdlib.h>
#include  <stdio.h>
#include  <iostream>
#include  <fstream>

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
#include  "PhaseOneSDDConnectionTable.h"
#include  "XrayLines.h"

#include  "MyRootFunctionLib_Basic.C"


void WriteParts2Rootfile(TString rootfile, TString place, Int_t adc_channel, Int_t divider, Int_t part){

   // writes parts of a given rootfile (like the 1st tenth of it) to another rootfile - currently only the tree and the ADC Branch filled with data only from sdd 1
 
    EventStruct  evt_r;
    EventStruct  evt;

    Int_t* ind_array;
    ind_array = GetTreeDivisionIndices(rootfile,place,divider,part);

    //TH1F *hev = new TH1F("hev", "Energy spectrum", 4096, -0.5, 4095.5);

    TString rootfilename;

    Int_t size = place.Sizeof();
    if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
    if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

// rootfile to read from
    TFile *fr = new TFile(rootfilename, "READ");

    TTree *tr_t = (TTree*)fr->Get("tr");

    TBranch *adcEvent   = tr_t->GetBranch("adc");
    adcEvent->SetAddress(evt_r.padc);

// rootfile to write to
    TFile *fw = new TFile( Form("/home/andreas/vip2/data/root/LNGS/20160211_%d.root",part), "RECREATE");
    fw->cd();

    TTree *tr = new TTree("tr", "SiPM and SDD data");
    tr->Branch("adc",   &evt.padc,   "ADC Channel/I");
 
    cout << ind_array[0] << " " << ind_array[1] << endl;


    for( Int_t i=ind_array[0]; i<ind_array[1]; i++ ){
        adcEvent->GetEntry(i);
        if( evt_r.padc[adc_channel] > 200 ){
            evt.padc[adc_channel] = evt_r.padc[adc_channel];
	    tr->Fill();
            
        }
    }

    fw->Write();
    //TFile hist_file(Form("/home/andreas/vip2/reports/1608_VIPReportLNF/sdd1_part%dLNGS.root",part), "CREATE");
    //hev->Write(); 
    //hist_file.Close();

    return;
   
}


void MakeQDCHistos(TString rootfile, TString place, Int_t qdc_channel){

  TH1F *total;
  TH1F *sdd_background;
  TH1F *scinti_all;
  TH1F *coinc;
  TH1F *scinti_segment;
  TH1F *sdd_cut;
  TH1F *sdd_tdc_cut;
  TH1F *qdc_cutoff;

  Int_t *adc_back = new Int_t();
  Int_t *scinti_coinc = new Int_t();
  Int_t *scinti_single = new Int_t();

  *adc_back = 0;
  *scinti_coinc = 0;
  *scinti_single = 0;

  TLegend *leg;

  Int_t size = place.Sizeof();
  Int_t tdc_channel;

  if(size == 4){tdc_channel = Qdc2Tdc_smi[qdc_channel];}
  if(size == 5){tdc_channel = Qdc2Tdc_lngs[qdc_channel];}
 

//  Int_t adc_channel = SDDToPadc[sdd];


  total = FillHistoWithCondition(rootfile, place, "qdc", qdc_channel, 0, 5000, qdc_channel); // select all events


  scinti_segment = FillHistoWithCondition(rootfile, place, "tdc", tdc_channel, 5000, 20000, qdc_channel);
 // scinti_all = FillHistoWithCondition(rootfile, place, "tdc", 13, 3800, 20000, qdc_channel);
  //coinc = FillHistoWithCondition(rootfile, place, "tdc", 13, 4300, 15000, qdc_channel);

 // TH1F *missed_events = (TH1F*)total->Clone("missed_events");
 // missed_events->Add(scinti_all,-1);

  sdd_cut = FillQdcHistoAdcCut(rootfile, place, qdc_channel, adc_back, scinti_coinc, scinti_single,0); // only with 0 as last parameter the 3 numbers are calculated
//  sdd_tdc_cut = FillQdcHistoAdcCut(rootfile, place, qdc_channel, adc_back, scinti_coinc, scinti_single,1);
  if(size==4){qdc_cutoff = MakeQdcCutoffHisto("smi",qdc_channel,10000);}
  if(size==5){qdc_cutoff = MakeQdcCutoffHisto("lngs",qdc_channel,10000);}


/*
  TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
  cc->cd();

  gPad->SetLogy();

  total->SetLineColor(1);
  total->Draw();

  scinti_all->SetLineColor(2);
  scinti_all->Draw("same");

  scinti_segment->SetLineColor(4);
  scinti_segment->Draw("same");

  missed_events->SetLineColor(7);
  missed_events->Draw("same");

  gStyle->SetOptStat(0);

  leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetHeader(Form("QDC %d spectrum with TDC cuts",qdc_channel));
  leg->AddEntry(total,"whole spectrum","l");
  leg->AddEntry(scinti_all,"entry in tdc[13]","l");
  leg->AddEntry(scinti_segment,"entry in layer-tdc channel","l");
  leg->AddEntry(missed_events,"no entry in tdc[13]","l");
  leg->Draw();

*/
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();

  gPad->SetLogy();

  total->SetLineColor(1);
  total->Draw();

  sdd_cut->SetLineColor(2);
  sdd_cut->Draw("same");

  qdc_cutoff->SetLineColor(4);
  qdc_cutoff->Draw("same");

  scinti_segment->SetLineColor(7);
  scinti_segment->Draw("same");

  gStyle->SetOptStat(0);

  leg = new TLegend(0.4,0.7,0.9,0.9);
  leg->SetHeader(Form("QDC %d spectrum with ADC cuts",qdc_channel));
  leg->AddEntry(qdc_cutoff,"QDC Cutoff");
  leg->AddEntry(total,"whole spectrum","l");
  leg->AddEntry(sdd_cut,"Events in SDD Background (> 7000 eV)","l");
  leg->AddEntry(scinti_segment,"Events in the coincidence timing of the scintillator segment","l");

//  leg->AddEntry((TObject*)0,Form("Events in SDD Background and over QDC Cutoff: %d",*adc_back),"");
//  leg->AddEntry((TObject*)0,Form("Events in SDD Background and over QDC Cutoff and in tdc[13]: %d",*scinti_coinc),"");
//  leg->AddEntry((TObject*)0,Form("Events in SDD Background and over QDC Cutoff and in tdc of layer: %d",*scinti_single),"");

  leg->Draw();

  c2->Print(Form("/home/andreas/vip2/reports/1602_BackgroundAnalysis2/QDC_with_ADC_Cut_LNGS/QDC_channel_%d.pdf",qdc_channel));
 // c2->Print("test.pdf");
  //c2->Print("/home/andreas/vip2/test.pdf");
  c2->Delete();

}




void PrintQDCHistos(){

  for(Int_t i = 0; i < 32; i++){

    MakeQDCHistos("BackRej_20160217_2.root","lngs",i);

  }

}

void AnalyseBackground(TString rootfile, TString place, Int_t part, Int_t lbound_ev, Int_t ubound_ev){

// this function doesnot work for data at LNGS taken before the 20.1.2016
// also for LNGS the tdc[3] channel is not included

  EventStruct  evt_r;

  TH1F *adc_mul_hist = new TH1F("adc_mul_hist","SDD multiplicity in background",10,0,10);
  TH1F *qdc_mul_hist = new TH1F("qdc_mul_hist","Scinti multiplicity in background",31,0,31);

  Int_t qdc_over_thresh[NQDC_CH] = {0};
  Int_t adc_channel;

  Int_t tdc_layer_scinti[8] = {0};
  Int_t tdc_layer_coinc[8] = {0};
  Int_t tdc_layer_sdd[8] = {0};

  Int_t tdc_layer_sdd_coinc_flag = 0; // any scinti layer tdc has coinc timing

  Int_t tdc_outer_coinc_flag = 0;
  Int_t tdc_inner_coinc_flag = 0;
  Int_t tdc_single_flag = 0;

  Int_t tdc_outer_flag;
  Int_t tdc_inner_flag;

  Int_t tdc_outer_scinti_flag = 0;
  Int_t tdc_inner_scinti_flag = 0;
  Int_t tdc_layer_scinti_flag = 0;

  Int_t tdc_layer_coinc_flag = 0; // tdc inner and outer layer coincidence
  
  Int_t tdc_scinti_flag = 0;
  Int_t tdc_coinc_flag = 0;
  Int_t tdc_sdd_flag = 0;

  Int_t qdc_outer_flag = 0;
  Int_t qdc_inner_flag = 0;
  Int_t qdc_coinc_flag = 0;
  Int_t qdc_flag = 0;
  Int_t qdc_single_flag = 0;

  Int_t qdc_mul;
  Int_t adc_mul;

  Int_t bin_ev[7];
  Double_t slope[7];
  Double_t offset[7];
  Int_t test;

  Int_t Tdc_Smi2Lngs[8] = {0, 1, 2, 3, 4, 5, 8, 10};
  Int_t dummie;
  
  Int_t event_counter = 0;

  Int_t qdc_single_counter = 0;
  Int_t qdc_coinc_counter = 0;
  Int_t qdc_counter = 0;

  Int_t tdc_single_counter = 0;
  Int_t tdc_coinc_counter = 0;
  Int_t tdc_scinti_counter = 0;
  Int_t tdc_sdd_counter = 0;
  Int_t tdc_layer_sdd_coinc_counter = 0;
  Int_t tdc_layer_coinc_counter = 0;
  Int_t tdc_layer_scinti_counter = 0;

//                            vorraussetzung -> das was gezählt wird
  Int_t dummie_counter1 = 0; // any qdc hit -> any tdc layer has coinc timing
  Int_t dummie_counter2 = 0; // any qdc hit -> any tdc layer has scinti timing
  Int_t dummie_counter3 = 0; // qdc coinc inner AND outer -> inner AND outer tdc layer show coincidence timing
  Int_t dummie_counter4 = 0; // qdc coinc inner AND outer -> tdc[13] has coinc timing
  Int_t dummie_counter5 = 0; // any tdc layer has coinc timing -> any qdc > thresh
  Int_t dummie_counter6 = 0; // any tdc layer has coinc timing -> tdc[13] has coinc timing
  Int_t dummie_counter7 = 0; // any tdc layer has coinc timing -> tdc[13] has SDD timing
  Int_t dummie_counter8 = 0; // any tdc layer has coinc timing -> inner AND outer tdc layer show coincidence timing
  Int_t dummie_counter12 = 0; // any tdc layer has coinc timing -> inner OR outer tdc layer show coincidence timing
  Int_t dummie_counter9 = 0; // outer AND inner tdc layer have coinc timing -> outer AND inner qdc layer have coinc timing
  Int_t dummie_counter10 = 0; // outer AND inner tdc layer have coinc timing -> tdc[13] has coinc timing
  Int_t dummie_counter11 = 0; // tdc[13] has SDD timing-> any tdc layer has coinc timing

 
  

  TString rootfilename;

  
  Int_t size = place.Sizeof();
  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;
 
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");
 
  TBranch *qdcEvent   = tree->GetBranch("qdc");
  qdcEvent->SetAddress(evt_r.qdc);
  TBranch *adcEvent   = tree->GetBranch("adc");
  adcEvent->SetAddress(evt_r.padc);
  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);


  Int_t nevent = tree->GetEntries();

  for(Int_t i = 0; i < 6; i++){

    slope[i+1] = GetSlope(rootfile,place,i+1);
    offset[i+1] = GetOffset(rootfile,place,i+1);

  }

  Int_t divider = 4;

  Int_t* ind_array;

  ind_array = GetTreeDivisionIndices(rootfile,place,divider,part);
  Int_t ind_lbound = ind_array[0];
  Int_t ind_ubound = ind_array[1];

  //for(Int_t i = ind_lbound; i < ind_ubound; i++){
  for(Int_t i = 0; i < nevent; i++){

      test = 0;

      tdcEvent->GetEvent(i);
      qdcEvent->GetEvent(i);
      adcEvent->GetEvent(i);

      adc_mul = 0;
      for(Int_t j = 0; j < 6; j++){ 
 
         
         adc_channel = SDDToPadc[j+1];
         
         bin_ev[j+1] = offset[j+1] + evt_r.padc[adc_channel] * slope[j+1];
         
         if(bin_ev[j+1] > lbound_ev && bin_ev[j+1] < ubound_ev){test = 1; // here is the condition for filling

           adc_mul += 1;

         }  
        //  cout << adc_channel << " " << bin_ev[j+1] << " " << size << " " << test << " " << slope[j+1] << " " << " " << offset[j+1] << " " << evt_r.padc[adc_channel] << endl;      
      }

      

     
     if(test == 1 && size == 4){// start of the loop for SMI

       adc_mul_hist->Fill(adc_mul);
       event_counter += 1;
       qdc_mul = 0;

       for(Int_t j = 0; j < NQDC_CH; j++){ // checking the QDC > threshold


         if(evt_r.qdc[j]>QdcCuts[j]){
 

            qdc_over_thresh[j]=1;
            qdc_mul += 1;
            qdc_flag = 1;

            if(Qdc2Outer[j] == 1){ qdc_outer_flag = 1; }
            if(Qdc2Outer[j] == 0){ qdc_inner_flag = 1; }

         }

       }

       qdc_mul_hist->Fill(qdc_mul);

       if(qdc_outer_flag == 1 && qdc_inner_flag == 1){ qdc_coinc_flag = 1;}
       if((qdc_outer_flag + qdc_inner_flag) == 1){ qdc_single_flag = 1; }

// ********************* setting the layer tdc arrays

       for(Int_t j = 0; j < 8; j++){ 

         if(evt_r.tdc[j]<2000){ tdc_layer_sdd[j] = 1; }
         if(evt_r.tdc[j]>2000 && evt_r.tdc[j]<5000){ tdc_layer_scinti[j] = 1; }
         if(evt_r.tdc[j]>5000){ tdc_layer_coinc[j] = 1; tdc_layer_sdd_coinc_flag = 1;}

         if(evt_r.tdc[j]>2000 && Tdc2Outer_smi[j] == 1){ tdc_outer_flag = 1; }
         if(evt_r.tdc[j]>2000 && Tdc2Outer_smi[j] == 0){ tdc_inner_flag = 1; }  


       }

// ******************** setting the scintillator tdc flags 

       if(tdc_outer_flag == 1 && tdc_inner_flag == 1) { tdc_layer_coinc_flag = 1; }  
       if((tdc_outer_flag + tdc_inner_flag ) == 1) { tdc_single_flag = 1; }

       if(evt_r.tdc[13] < 2000) { tdc_sdd_flag = 1; }
       if(evt_r.tdc[13] > 2000 && evt_r.tdc[13] < 5000) { tdc_scinti_flag = 1; }
       if(evt_r.tdc[13] > 5000) { tdc_coinc_flag = 1; }

 //************************* setting the counters

       if(qdc_flag == 1) {  qdc_counter += 1;  if(tdc_layer_sdd_coinc_flag == 1) { dummie_counter1 += 1; }} 
       if(qdc_coinc_flag == 1) { qdc_coinc_counter += 1; if(tdc_coinc_flag == 1) { dummie_counter4 += 1; } }
       if(qdc_single_flag == 1) { qdc_single_counter += 1; }

       if(tdc_single_flag == 1) { tdc_single_counter += 1; if(qdc_flag == 1) { dummie_counter8 += 1; } }
       if(tdc_coinc_flag == 1) { tdc_coinc_counter += 1; }

       if(tdc_scinti_flag == 1) { tdc_scinti_counter += 1; }

       if(tdc_layer_sdd_coinc_flag == 1) { tdc_layer_sdd_coinc_counter += 1; if(qdc_flag == 1) { dummie_counter2 += 1; } if( tdc_coinc_flag == 1 ) { dummie_counter3 +=1; } 
            if(tdc_sdd_flag == 1) { dummie_counter7 += 1; } }

       if(tdc_sdd_flag == 1) { tdc_sdd_counter += 1; if(tdc_layer_sdd_coinc_flag == 1) {dummie_counter5 += 1;} if( qdc_flag == 1 ) { dummie_counter6 += 1; }}

       if(tdc_layer_coinc_flag == 1){ tdc_layer_coinc_counter += 1; }


    for(Int_t k = 0; k < NQDC_CH; k++){
      qdc_over_thresh[k] = 0;
    }

    for(Int_t k = 0; k < 8; k++){
      tdc_layer_scinti[k] = 0;
      tdc_layer_coinc[k] = 0;
      tdc_layer_sdd[k] = 0;
    }

    tdc_outer_flag = 0;
    tdc_inner_flag = 0;
    tdc_single_flag = 0;
    tdc_layer_coinc_flag = 0; 
    tdc_layer_sdd_coinc_flag = 0;
  
    tdc_scinti_flag = 0;
    tdc_coinc_flag = 0;
    tdc_sdd_flag = 0;

    qdc_outer_flag = 0;
    qdc_inner_flag = 0;
    qdc_coinc_flag = 0;
    qdc_flag = 0;
    qdc_single_flag = 0;

    tdc_outer_coinc_flag = 0;
    tdc_inner_coinc_flag = 0;

    tdc_outer_scinti_flag = 0;
    tdc_inner_scinti_flag = 0;
    tdc_layer_scinti_flag = 0;

  Int_t tdc_layer_coinc_flag = 0; // tdc inner and outer layer coincidence
  
   } // end of if loop with background test

   if(test == 1 && size == 5){  // start of the if loop with background test for LNGS

       adc_mul_hist->Fill(adc_mul);
       event_counter += 1;
       qdc_mul = 0;
       Int_t tdc_channel;

       for(Int_t j = 0; j < NQDC_CH; j++){ // checking the QDC > threshold

        
           if(evt_r.qdc[j]>QdcCuts_100nsGate[j]){

            qdc_over_thresh[j]=1;
            qdc_mul += 1;
            qdc_flag = 1;

            if(Qdc2Outer[j] == 1){ qdc_outer_flag = 1; }
            if(Qdc2Outer[j] == 0){ qdc_inner_flag = 1; }

          }

       }

       qdc_mul_hist->Fill(qdc_mul);

       if(qdc_outer_flag == 1 && qdc_inner_flag == 1){ qdc_coinc_flag = 1;}
       if((qdc_outer_flag + qdc_inner_flag) == 1){ qdc_single_flag = 1; }

// ********************* setting the layer tdc arrays 

       for(Int_t j = 0; j < 8; j++){ 

         if( j == 3 ){continue;} // this tdc channel is broken ------  change here

         tdc_channel = Tdc_Smi2Lngs[j]; 

         if(evt_r.tdc[tdc_channel]<2000){ tdc_layer_sdd[j] = 1; } // layer tdc has sdd timing
         if(evt_r.tdc[tdc_channel]>2000 && evt_r.tdc[tdc_channel]<7000){ tdc_layer_scinti[j] = 1; } // layer tdc has scinti only timing
         if(evt_r.tdc[tdc_channel]>7000){ tdc_layer_coinc[j] = 1; tdc_layer_sdd_coinc_flag = 1;} // (any) layer tdc has coinc timing
 
      //   if( j == 3 && evt_r.tdc[tdc_channel] > 5000) { cout << j << " " << tdc_channel << " " << evt_r.tdc[tdc_channel] << " " << tdc_layer_sdd_coinc_flag << endl;}

         if(evt_r.tdc[tdc_channel]>5000 && Tdc2Outer_smi[j] == 1){ tdc_outer_coinc_flag = 1; } // any outer/inner layer tdc shows coincidence timing 
         if(evt_r.tdc[tdc_channel]>5000 && Tdc2Outer_smi[j] == 0){ tdc_inner_coinc_flag = 1; }  

         if(evt_r.tdc[tdc_channel]>2000 && evt_r.tdc[tdc_channel]<5000 && Tdc2Outer_smi[j] == 1){ tdc_outer_scinti_flag = 1; } // any outer/inner layer tdc shows scintillator timing 
         if(evt_r.tdc[tdc_channel]>2000 && evt_r.tdc[tdc_channel]<5000 && Tdc2Outer_smi[j] == 0){ tdc_inner_scinti_flag = 1; }
         if(evt_r.tdc[tdc_channel]>2000 && evt_r.tdc[tdc_channel]<5000){ tdc_layer_scinti_flag = 1; }  // any tdc layer shows scinti timing 


       }

// ******************** setting the scintillator tdc flags 

       if(tdc_outer_coinc_flag == 1 && tdc_inner_coinc_flag == 1) { tdc_layer_coinc_flag = 1; }  // one or more outer AND one or more inner layers show coincidence timing
       if((tdc_outer_coinc_flag + tdc_inner_coinc_flag ) == 1) { tdc_single_flag = 1; } // one or more outer OR one or more inner layers show coincidence timing

       if(evt_r.tdc[13] < 2000) { tdc_sdd_flag = 1; } // tdc 13 has sdd timing
       if(evt_r.tdc[13] > 2000 && evt_r.tdc[13] < 5000) { tdc_scinti_flag = 1; } // tdc 13 has scintillator timing
       if(evt_r.tdc[13] > 5000) { tdc_coinc_flag = 1; } // tdc 13 has coincidence timing

 //************************* setting the counters

       if(qdc_flag == 1) {  qdc_counter += 1;  if(tdc_layer_sdd_coinc_flag == 1) { dummie_counter1 += 1; } if(tdc_layer_scinti_flag == 1) {dummie_counter2 += 1; }} 

       if(qdc_coinc_flag == 1) { qdc_coinc_counter += 1; if(tdc_coinc_flag == 1) { dummie_counter4 += 1; } if(tdc_layer_coinc_flag == 1) {dummie_counter3 += 1; }}

       if(qdc_single_flag == 1) { qdc_single_counter += 1; }

       if(tdc_layer_sdd_coinc_flag == 1) { tdc_layer_sdd_coinc_counter += 1; if(qdc_flag == 1) { dummie_counter5 += 1; } if( tdc_coinc_flag == 1 ) { dummie_counter6 +=1; } }
       if(tdc_sdd_flag == 1){ { dummie_counter7 += 1; } if(tdc_layer_coinc_flag == 1) { dummie_counter8 +=1; } if( tdc_single_flag == 1 ) { dummie_counter12 += 1; }}

       if( tdc_layer_coinc_flag == 1){ tdc_layer_coinc_counter += 1;  if(qdc_coinc_flag == 1) {dummie_counter9 += 1;} if(tdc_coinc_flag == 1) {dummie_counter10 +=1;} }



//       if(tdc_single_flag == 1) { tdc_single_counter += 1; }
       if(tdc_coinc_flag == 1) { tdc_coinc_counter += 1; }
       if(tdc_scinti_flag == 1) { tdc_scinti_counter += 1; }
       if(tdc_sdd_flag == 1) { tdc_sdd_counter += 1; if(tdc_layer_sdd_coinc_flag == 1) {dummie_counter11 += 1;} }
       if(tdc_layer_scinti_flag == 1) { tdc_layer_scinti_counter += 1; }
 



    for(Int_t k = 0; k < NQDC_CH; k++){
      qdc_over_thresh[k] = 0;
    }

    for(Int_t k = 0; k < 8; k++){
      tdc_layer_scinti[k] = 0;
      tdc_layer_coinc[k] = 0;
      tdc_layer_sdd[k] = 0;
    }

    tdc_outer_flag = 0;
    tdc_inner_flag = 0;
    tdc_single_flag = 0;
    tdc_layer_coinc_flag = 0; 
    tdc_layer_sdd_coinc_flag = 0;
  
    tdc_scinti_flag = 0;
    tdc_coinc_flag = 0;
    tdc_sdd_flag = 0;

    tdc_outer_coinc_flag = 0;
    tdc_inner_coinc_flag = 0;

    tdc_outer_scinti_flag = 0;
    tdc_inner_scinti_flag = 0;
    tdc_layer_scinti_flag = 0;


    qdc_outer_flag = 0;
    qdc_inner_flag = 0;
    qdc_coinc_flag = 0;
    qdc_flag = 0;
    qdc_single_flag = 0;



   }


  } // end of the loop over ALL events

  cout << endl;

  cout << "The number of background events between " << lbound_ev << " eV and " << ubound_ev << " eV is: " << event_counter << endl;
  cout << "Among these events, the numbers are as follows:" << endl << endl;

  cout << "The number of qdc events (any scint > thresh) is: " << qdc_counter << endl;
  cout << "From these, there are: " << dummie_counter1 << " events in any scintillator tdc layer with coincidence timing = " << 100 *((float)dummie_counter1/(float)qdc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter2 << " events in any scintillator tdc layer with scintillator timing = " << 100 *((float)dummie_counter2/(float)qdc_counter) << "%" << endl << endl;

  cout << "The number of qdc coincidence events (inner AND outer layer) is: " << qdc_coinc_counter << endl;
  cout << "From these, there are: " << dummie_counter3 << " events with coincidence timing in inner AND outer layer = " << 100*((float)dummie_counter3/(float)qdc_coinc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter4 << " events with tdc[13] coincidence timing = " << 100*((float)dummie_counter4/(float)qdc_coinc_counter) << "%" << endl << endl;


  cout << "The number of events with coincidence timing in any scintillator tdc layer is: " << tdc_layer_sdd_coinc_counter << endl;
  cout << "From these, there are: " << dummie_counter5 << " events  with any qdc > thresh = " << 100*((float)dummie_counter5/(float)tdc_layer_sdd_coinc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter6 << " events  with tdc[13] coinc timing = " << 100*((float)dummie_counter6/(float)tdc_layer_sdd_coinc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter7 << " events  with tdc[13] SDD timing = " << 100*((float)dummie_counter7/(float)tdc_layer_sdd_coinc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter12 << " events  with inner OR outer tdc layer coinc timing = " << 100*((float)dummie_counter12/(float)tdc_layer_sdd_coinc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter8 << " events with coinc timing in inner AND outer tdc layer = " << 100*((float)dummie_counter8/(float)tdc_layer_sdd_coinc_counter) << "%" << endl << endl;

  cout << "The number of qdc single layer (exactly one layer) events is: " << qdc_single_counter << endl << endl;

  cout << "The number of events with coincidence timing in inner AND outer layer is: " << tdc_layer_coinc_counter << endl;
  cout << "From these, there are: " << dummie_counter9 << " events  with inner AND outer qdc layer hit = " << 100*((float)dummie_counter9/(float)tdc_layer_coinc_counter) << "%" << endl;
  cout << "From these, there are: " << dummie_counter10 << " events  with tdc[13] coinc timing = " << 100*((float)dummie_counter10/(float)tdc_layer_coinc_counter) << "%" << endl << endl;
 

  cout << "The number of tdc coincidence events (in tdc[13]) is: " << tdc_coinc_counter << endl;
  cout << "The number of tdc scintillator events (in tdc[13]) is: " << tdc_scinti_counter << " ...this should be 0!" << endl << endl;

  cout << "The number of tdc SDD events (in tdc[13]) is: " << tdc_sdd_counter << endl;
  cout << "From these, there are: " << dummie_counter11 << " events with any tdc layer having coinc timing" << endl << endl;

  cout << "The number of events with any tdc layer having scintillator timing is : " << tdc_layer_scinti_counter << endl << endl;

  cout << "The rejection ratio is around " << 100*((float)tdc_layer_sdd_coinc_counter/(float)event_counter) << "%" << endl;
  cout << "... in part " << part << " of " << divider << endl << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  adc_mul_hist->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();
  qdc_mul_hist->Draw();

/*                            Vorraussetung -> das was gezählt wird
  Int_t dummie_counter1 = 0; // any qdc hit -> any tdc layer has coinc timing
  Int_t dummie_counter2 = 0; // any qdc hit -> any tdc layer has scinti timing
  Int_t dummie_counter3 = 0; // qdc coinc inner AND outer -> inner AND outer tdc layer show coincidence timing
  Int_t dummie_counter4 = 0; // qdc coinc inner AND outer -> tdc[13] has coinc timing
  Int_t dummie_counter5 = 0; // any tdc layer has coinc timing -> any qdc > thresh
  Int_t dummie_counter6 = 0; // any tdc layer has coinc timing -> tdc[13] has coinc timing
  Int_t dummie_counter7 = 0; // any tdc layer has coinc timing -> tdc[13] has SDD timing
  Int_t dummie_counter8 = 0; // any tdc layer has coinc timing -> inner AND outer tdc layer show coincidence timing
  Int_t dummie_counter12 = 0; // any tdc layer has coinc timing -> inner OR outer tdc layer show coincidence timing
  Int_t dummie_counter9 = 0; // outer AND inner tdc layer have coinc timing -> outer AND inner qdc layer have coinc timing
  Int_t dummie_counter10 = 0; // outer AND inner tdc layer have coinc timing -> tdc[13] has coinc timing
  Int_t dummie_counter11 = 0; // tdc[13] has SDD timing-> any tdc layer has coinc timing
*/


}



void AnalyseSpectrum(TString rootfile, TString place, Int_t sdd){

// this function prints the amount of counts in the ROI, the Cu,Mn, Ti line!
  TH1F *histo;
  //histo = FillScaledHistogram(rootfile, place, sdd); // change here!!

//  histo->Draw();

  Int_t lbound = 7629;
  Int_t ubound = 7829;
  

  TAxis *axis = histo->GetXaxis();

  Int_t bmin = axis->FindBin(lbound);
  Int_t bmax = axis->FindBin(ubound);
  

  Int_t counts_back = histo->Integral(bmin,bmax);

  counts_back -= histo->GetBinContent(bmin)*(lbound-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
  counts_back -= histo->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-ubound)/axis->GetBinWidth(bmax);
  

  cout << "The number of background counts in the ROI is: " << counts_back << endl;


  Int_t lbound_cu = 7947;
  Int_t ubound_cu = 8147;  

  Int_t bmin_cu = axis->FindBin(lbound_cu);
  Int_t bmax_cu = axis->FindBin(ubound_cu);
  

  Int_t counts_cu = histo->Integral(bmin_cu,bmax_cu);

  counts_cu -= histo->GetBinContent(bmin_cu)*(lbound_cu-axis->GetBinLowEdge(bmin_cu))/axis->GetBinWidth(bmin_cu);
  counts_cu -= histo->GetBinContent(bmax_cu)*(axis->GetBinUpEdge(bmax_cu)-ubound_cu)/axis->GetBinWidth(bmax_cu);

  counts_cu = counts_cu - counts_back;

  cout << "The number of Copper counts is: " << counts_cu << endl;


  Int_t lbound_mn = 5798;
  Int_t ubound_mn = 5998;
  
  Int_t bmin_mn = axis->FindBin(lbound_mn);
  Int_t bmax_mn = axis->FindBin(ubound_mn);
  

  Int_t counts_mn = histo->Integral(bmin_mn,bmax_mn);

  counts_mn -= histo->GetBinContent(bmin_mn)*(lbound_mn-axis->GetBinLowEdge(bmin_mn))/axis->GetBinWidth(bmin_mn);
  counts_mn -= histo->GetBinContent(bmax_mn)*(axis->GetBinUpEdge(bmax_mn)-ubound_mn)/axis->GetBinWidth(bmax_mn);

  cout << "The number of Mn counts is: " << counts_mn << endl;


  Int_t lbound_ti = 4410;
  Int_t ubound_ti = 4610;
  
  Int_t bmin_ti = axis->FindBin(lbound_ti);
  Int_t bmax_ti = axis->FindBin(ubound_ti);
  

  Int_t counts_ti = histo->Integral(bmin_ti,bmax_ti);

  counts_ti -= histo->GetBinContent(bmin_ti)*(lbound_ti-axis->GetBinLowEdge(bmin_ti))/axis->GetBinWidth(bmin_ti);
  counts_ti -= histo->GetBinContent(bmax_ti)*(axis->GetBinUpEdge(bmax_ti)-ubound_ti)/axis->GetBinWidth(bmax_ti);

  cout << "The number of Ti counts is: " << counts_ti << endl;



  TH1F *histo_unsalced;
  histo_unscaled = HistoFromTree(rootfile, place, sdd);
  TAxis *axis_unscaled = histo_unscaled->GetXaxis();

  Int_t bmin_over = axis_unscaled->FindBin(3850);
  Int_t bmax_over = axis_unscaled->FindBin(4100);


  Int_t counts_over = histo_unscaled->Integral(bmin_over,bmax_over);

  cout << "The number of Overflow counts is: " << counts_over << endl;

  // ------------ addition

  Int_t lbound1 = 500;
  Int_t ubound1 = 4000;
  
  Int_t bmin1 = axis->FindBin(lbound1);
  Int_t bmax1 = axis->FindBin(ubound1);
  

  Int_t counts1 = histo->Integral(bmin1,bmax1);

  cout << "The number of Hz in region 1 is: " << (float)counts1*1000/(float)18900 << endl;

  Int_t lbound2 = 4000;
  Int_t ubound2 = 7000;
  
  Int_t bmin2 = axis->FindBin(lbound2);
  Int_t bmax2 = axis->FindBin(ubound2);
  

  Int_t counts2 = histo->Integral(bmin2,bmax2);

  cout << "The number of Hz in region 2 is: " << (float)counts2*1000/(float)18900 << endl;

  Int_t lbound3 = 7000;
  Int_t ubound3 = 40000;
  
  Int_t bmin3 = axis->FindBin(lbound3);
  Int_t bmax3 = axis->FindBin(ubound3);
  

  Int_t counts3 = histo->Integral(bmin3,bmax3);

  cout << "The number of Hz in region 3 is: " << (float)counts3*1000/(float)18900 << endl;

  cout << "Combined Rate: " << (float)(counts1+counts2+counts3)*1000/18900 << endl;


}


void GetRateParts(TString rootfile, Int_t divider, TString place){

  Int_t size = place.Sizeof();
  TString rootfilename;

  EventStruct  evt_r;

  Double_t rate_final;
  Int_t* ind_array;
  Int_t ind_lbound;
  Int_t ind_ubound;

  Int_t part;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  for(Int_t i = 1; i < divider+1; i++){

    part = i;

    ind_array =  GetTreeDivisionIndices(rootfile, place, divider, part);

    ind_lbound = ind_array[0];
    ind_ubound = ind_array[1];
  
    rate_final = GetRate(rootfile, ind_lbound, ind_ubound, place);

    cout << rate_final << endl;
  }

}
/*
void MakeROIPlots(TString rootfile, TString place, Int_t sdd){

  TH1F *hist;
  hist=FillScaledHistogram(rootfile,place,sdd); // 20 eV binning originally

  hist->Rebin(3);
  hist->GetXaxis()->SetRangeUser(7000,10000);
  hist->GetYaxis()->SetTitle("Counts/60 eV");
  hist->GetXaxis()->SetTitle("Energy in eV");
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->SetTitle(Form("SDD %d energy spectrum",sdd));
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  

  hist->Draw();
  
  //c1->Print(Form("/home/andreas/vip2/reports/1608_VIPReportLNF/sdd%d-60eVbin-7-10.png",sdd));

  c1->Close();


}*/

TH1F* FillSummedEnergyHisto(TString rootfile, TString place, Int_t low_edge_ev, Int_t ch_number){

 Int_t ev_range = 3200;
 Int_t bin_nr;
 Double_t bin_cont;

 TH1F *histo1;
 TH1F *histo2;
 TH1F *histo3;
 TH1F *histo4;
 TH1F *histo5;
 TH1F *histo6;

 histo1=FillScaledHistogram(rootfile,place,1,low_edge_ev,ch_number);
 histo1->Scale(1,"width");
 histo2=FillScaledHistogram(rootfile,place,2,low_edge_ev,ch_number);
 histo2->Scale(1,"width");
 histo3=FillScaledHistogram(rootfile,place,3,low_edge_ev,ch_number);
 histo3->Scale(1,"width");
 histo4=FillScaledHistogram(rootfile,place,4,low_edge_ev,ch_number);
 histo4->Scale(1,"width");
 histo5=FillScaledHistogram(rootfile,place,5,low_edge_ev,ch_number);
 histo5->Scale(1,"width");
 histo6=FillScaledHistogram(rootfile,place,6,low_edge_ev,ch_number);
 histo6->Scale(1,"width");



 TH1F *histo_ev = new TH1F("histo_ev","Energy Spectrum",ev_range,low_edge_ev,low_edge_ev+ev_range);

 for(Int_t ev = low_edge_ev; ev < 9600; ev++){

	bin_cont = 0;

	bin_nr = histo1->FindBin(ev);
	bin_cont = histo1->GetBinContent(bin_nr);

	bin_nr = histo2->FindBin(ev);
	bin_cont = bin_cont + histo2->GetBinContent(bin_nr);

	bin_nr = histo3->FindBin(ev);
	bin_cont = bin_cont + histo3->GetBinContent(bin_nr);

	bin_nr = histo4->FindBin(ev);
	bin_cont = bin_cont + histo4->GetBinContent(bin_nr);

	bin_nr = histo5->FindBin(ev);
	bin_cont = bin_cont + histo5->GetBinContent(bin_nr);

	bin_nr = histo6->FindBin(ev);
	bin_cont = bin_cont + histo6->GetBinContent(bin_nr);
	//cout << ev << " " << bin_nr << " " << bin_cont << endl;

	histo_ev->SetBinContent(ev-low_edge_ev,bin_cont);
 }

 gStyle->SetOptStat(0);
 histo_ev->GetXaxis()->SetTitle("Energy [eV]");
 histo_ev->Draw();

 return histo_ev;

}

void MakeStabilityPlots(){

 

  Double_t x[8] ={1.,2.,3.,4.,5.,6.,7.,8.};
  Double_t x_err[8] = {0.};

  Double_t mnka[8] = {771.65,770.98,770.65,770.16,775.36,775.05,775.13,768.2};
  Double_t mnka_err[8] = {0.103,0.086,0.085,0.085,0.091,0.086,0.091,0.085};
  Double_t fwhm[8] = {152.48,153.24,150.49,151.47,154.79,154.08,153.31,156.58};
  Double_t fwhm_err[8] = {0.104,0.088,0.086,0.085,0.093,0.088,0.092,0.089};

  TGraphErrors *gmnkaerr = new TGraphErrors(8, x, mnka, x_err, mnka_err );
  TGraphErrors *gfwhmerr = new TGraphErrors(8, x, fwhm, x_err, fwhm_err );

    gmnkaerr->SetMarkerStyle(20);
    gmnkaerr->SetMarkerColor(4);
    gmnkaerr->SetMarkerSize(1.3);
    gmnkaerr->GetYaxis()->SetTitle("ADC channel");
    gmnkaerr->GetYaxis()->SetTitleOffset(1.4);
    gmnkaerr->SetTitle("Peak Stability of the Mn K-alpha line");

    gfwhmerr->SetMarkerStyle(20);
    gfwhmerr->SetMarkerColor(4);
    gfwhmerr->SetMarkerSize(1.3);
    gfwhmerr->GetYaxis()->SetTitle("FWHM in eV");
    gfwhmerr->GetYaxis()->SetTitleOffset(1.4);
    gfwhmerr->SetTitle("FWHM of Mn K-alpha in eV");


  TCanvas *cmn = new TCanvas("cmn", "cmn", 1, 1, 900, 1000);
  cmn->cd();
  
  gmnkaerr->Draw("ap");

  TCanvas *cfwhm = new TCanvas("cfwhm", "cfwhm", 1, 1, 900, 1000);
  cfwhm->cd();
  
  gfwhmerr->Draw("ap");
}

