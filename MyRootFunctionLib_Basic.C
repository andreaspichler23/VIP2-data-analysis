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


using namespace std;

/* a few lines to rule them all:
TFile *f = new TFile("BackRej_1023.root")
tr->Draw("adc[3]>>htemp(500,0,4200)","tdc[13]>5000 && adc[3]>80")
tr->Draw("adc[3]>>htemp1(500,0,4200)","tdc[13]>5000", "same")
*/



TH1F* HistoFromTree(TString rootfile, TString place, Int_t sdd){ // not working for SMI because of the adc -> padc TBranch  name and size change

  EventStruct  evt_r;
  Int_t adc_channel_list[6] = {0,2,3,5,6,7};
  Int_t adc_channel = adc_channel_list[sdd-1];

  Int_t size = place.Sizeof();
  TString rootfilename;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  //cout << rootfilename << endl;

  //TFile *f = new TFile(rootfilename, "READ");
  TFile *f = TFile::Open(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");
  
  TBranch *adcEvent;


  //cout << 1 << endl;
  //adcEvent = tree->GetBranch("padc");
  adcEvent = tree->GetBranch("adc");
  adcEvent->SetAddress(evt_r.padc);

  Int_t nevent = tree->GetEntries();
  //cout << nevent << endl;
  //cout << 2 << endl;
  TH1F *spect = new TH1F("spect","", 1024, 0, 4096);


  for( Int_t i = 0; i < nevent; i++){
    //cout << "getevent:" << endl;
    adcEvent->GetEvent(i);
    //cout << "fill:" << endl;
    spect->Fill(evt_r.padc[adc_channel]);
 
  }
  
  spect->SetDirectory(0); // to decouple it from the open file directory -> otherwise when you close the rootfile, the histogram is also gone
  
  f->Close();
  delete f;
  //spect->Draw();

  //cout << 3;
  
  
 return spect;

}


void SaveFinalHisto(TString inFile, TString outFile, TString place){ 
    // saving histograms from tree of final rootfile to another rootfile with scintillator rejection cut

  
  Int_t adcChannelList[6] = {0,2,3,5,6,7};
  Double_t energyList[16];
  Short_t trgid;
  

  Int_t size = place.Sizeof();
  TString inFileName, outFileName;

  if(size==4)inFileName = ROOT_PATH_SMI + "/" + inFile;
  if(size==5)inFileName = ROOT_PATH_LNGS + "/1-618files-final/" + inFile + ".root";
  
  if(size==4)outFileName = ROOT_PATH_SMI + "/" + outFile;
  if(size==5)outFileName = ROOT_PATH_LNGS + "/1-618files-final/" + outFile + ".root";

  //cout << rootfilename << endl;


  TFile *inF = TFile::Open(inFileName, "READ");
  TFile *outF = TFile::Open(outFileName, "UPDATE");
  TTree *tree = (TTree*)inF->Get("tr");
  
  TBranch *energyB;
  TBranch *trgidB;
  
  TString name1 = inFile + "sdd1";
  TString name2 = inFile + "sdd2";
  TString name3 = inFile + "sdd3";
  TString name4 = inFile + "sdd4";
  TString name5 = inFile + "sdd5";
  TString name6 = inFile + "sdd6";
 

  //adcEvent = tree->GetBranch("padc");
  energyB = tree->GetBranch("energy");
  trgidB = tree->GetBranch("trgid");
  
  energyB->SetAddress(energyList);
  trgidB->SetAddress(&trgid);

  Int_t nevent = tree->GetEntries();

  TH1F *sdd1H = new TH1F(name1,"", 9000, 1000, 10000);
  TH1F *sdd2H = new TH1F(name2,"", 9000, 1000, 10000);
  TH1F *sdd3H = new TH1F(name3,"", 9000, 1000, 10000);
  TH1F *sdd4H = new TH1F(name4,"", 9000, 1000, 10000);
  TH1F *sdd5H = new TH1F(name5,"", 9000, 1000, 10000);
  TH1F *sdd6H = new TH1F(name6,"", 9000, 1000, 10000);


  for( Int_t i = 0; i < nevent; i++){
    
    energyB->GetEvent(i);
    trgidB->GetEvent(i);
    
    if( trgid == 1 ){
        
        sdd1H->Fill(energyList[0]);
        sdd2H->Fill(energyList[2]);
        sdd3H->Fill(energyList[3]);
        sdd4H->Fill(energyList[5]);
        sdd5H->Fill(energyList[6]);
        sdd6H->Fill(energyList[7]);
        
    }
    
 
  }
  
  
  //sdd1H->Draw();
  outF->cd();
  //gDirectory->pwd();
  //gDirectory->ls();
  outF->Write();
  
  inF->Close();
  outF->Close();
  
  
  
 return;

}

Double_t GetOffset(TString rootfile, TString place, Int_t sdd){

  TH1F *histo;
  
  histo = HistoFromTree(rootfile, place, sdd);
  
  histo -> GetXaxis()->SetRangeUser(200,3000);
  

  TSpectrum *s = new TSpectrum();
  s->Search(histo, 2, "nodraw", 0.05 );  // 2. attribute is sigma                    

  Int_t peakN = s->GetNPeaks();

//  cout << peakN << endl;

  Double_t *pX = s->GetPositionX();
      
     
  Int_t idx[peakN]; // thats an index array
  TMath::Sort(peakN,pX,idx,0); // idx now contains indices of sorted array
//  cout << pX[idx[0]] << " " << pX[idx[1]] <<  " " << pX[idx[2]] << endl;
 
  Double_t slope = (MnKa1-TiKa1)/(pX[idx[2]]-pX[idx[0]]);
  Double_t offset = TiKa1 - slope * pX[idx[0]];

  return offset;

}
Double_t GetSlope(TString rootfile, TString place, Int_t sdd){

  TH1F *histo;
  
  histo = HistoFromTree(rootfile, place, sdd);

  histo -> GetXaxis()->SetRangeUser(200,3000);
  

  TSpectrum *s = new TSpectrum();

  s->Search(histo, 2, "nodraw", 0.05 );  // 2. attribute is sigma                    

  Int_t peakN = s->GetNPeaks();

//  cout << peakN << endl;

  Double_t *pX = s->GetPositionX();
      
     
  Int_t idx[peakN]; // thats an index array
  TMath::Sort(peakN,pX,idx,0); // idx now contains indices of sorted array
//  cout << pX[idx[0]] << " " << pX[idx[1]] <<  " " << pX[idx[2]] << endl;
 
  Double_t slope = (MnKa1-TiKa1)/(pX[idx[2]]-pX[idx[0]]);
  Double_t offset = TiKa1 - slope * pX[idx[0]];
  
  return slope;

}

int* GetTreeDivisionIndices(TString rootfile, TString place, Int_t divider, Int_t part){

// this function returns start and end indices for a part of the given rootfile on the base of the Mn k-alpha count rate of SDD 1
    
  

  Double_t slope;
  Double_t offset;
  Int_t bin_ev, bin_ch;
  Int_t mn_counter = 0;
 
  static int ind_array[2];

  EventStruct  evt_r;

  Int_t size = place.Sizeof();
  TString rootfilename;

  
  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  //TFile *f = new TFile(rootfilename, "READ");
  TFile *f = TFile::Open(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TBranch *adcEvent   = tree->GetBranch("adc");
  adcEvent->SetAddress(evt_r.padc);  

  slope = GetSlope( rootfile, place, 1);
  offset = GetOffset( rootfile, place, 1);

  Int_t nevent = tree->GetEntries();

  for( Int_t i = 0; i < nevent; i++){

      adcEvent->GetEvent(i);
      bin_ch = evt_r.padc[0];
      bin_ev = offset + bin_ch * slope;

      if(bin_ev > 5748 && bin_ev < 6048){ mn_counter += 1; } 
  
  }


 Int_t mn_event_lb = (part - 1) * mn_counter/divider;
 Int_t mn_event_ub = part * mn_counter/divider;
 
 if( part > divider ){ cout << " PART > DIVIDER; THIS SHOULD NOT BE !!!"  << endl; goto end;}
 
 mn_counter = 0;

 for( Int_t j = 0; j < nevent; j++ ){

    adcEvent->GetEvent(j);
    bin_ch = evt_r.padc[0];
    bin_ev = offset + bin_ch * slope;

    if(bin_ev > 5748 && bin_ev < 6048){ mn_counter += 1; }

    if( mn_counter == mn_event_lb ) { ind_array[0] = j; }
    if( mn_counter == mn_event_ub ) { ind_array[1] = j; break;}

 }

// cout<< nevent << " " << ind_array[0] << " " << ind_array[1] << endl;
 end:
 f->Close();
 delete f;
 
 
 return ind_array;

}



Int_t GetRunTime(TString rootfile, TString place, Int_t divider, Int_t part, Int_t complete){

// this function returns the runtime of the given rootfile on the base of the Labview millisecond clock - it goes through the events event by event and looks at the ms time tags
// these time tags should be steadily rising
// if the new event takes place x<10 seconds after the previous one,  x seconds are added to the total count
// it counts the seconds in which some events happens - if there is no event for 10 seconds or more, it is counted as a gap and this time is not taken into account
// also events which have some wrong timing information are discarded for this analysis
// the labview ms clock takes values from -2*10^9 - 2 * 10^9 (= 4*10^6 seconds ~ 46 days 7 h 6 min ) ... in that region consecutive events have consecutive timestamps 
// once the clock reaches the upper limit it resets back to the negative limit
// useful: tr->Draw("evid:clk>>htmp(1000,0,1000000000,1000,0,50000)") .... limits for ms clock to be changed

    
    // update: switch to milliseconds as the time unit 
    //update 2: switch to being able to get run time of parts of the rootfile

  /*Double_t slope;
  Double_t offset;
  Double_t bin_ev; 
  Int_t bin_ch;
  Int_t ti_counter = 0;
  Double_t ti_rate = 0.107125;*/
  Double_t sec;
  Int_t hour, day, min;

  Int_t sec_curr = 0, sec_buff = 0, sec_counter = 0;
  Int_t fail_counter = 0;
  Int_t t_gap;
 
  Int_t* ind_array;
  
  
  Int_t lbound = 0;
  Int_t ubound = 0;

  EventStruct  evt_r;

  Int_t size = place.Sizeof();
  TString rootfilename;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;
  
  //cout << rootfilename << endl;

  //TFile *f = new TFile(rootfilename, "READ");
  TFile *f = TFile::Open(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  //TBranch *adcEvent   = tree->GetBranch("adc");
  //adcEvent->SetAddress(evt_r.padc);  

  TBranch *clkEvent   = tree->GetBranch("clk");
  clkEvent->SetAddress(&evt_r.clk);

  
  //cout << "starting time: " << sec_buff << endl;

  Int_t nevent = tree->GetEntries();
  
  if( complete == 1 ){ lbound = 1; ubound = nevent; }
  if( complete == 0 ){ ind_array =  GetTreeDivisionIndices(rootfile, place, divider, part); 
  
    lbound = ind_array[0];
    ubound = ind_array[1];
  
  }
  
  clkEvent->GetEvent(lbound);
  sec_buff = (int)evt_r.clk/1000 - 1;
  

  //cout << "lbound: " << lbound << " ubound: " << ubound << " nevent: " << nevent << endl;

  for( Int_t i = lbound; i <= ubound; i++){

    clkEvent->GetEvent(i);
    sec_curr = (int)evt_r.clk/1000; // second of current event

     

    t_gap = sec_curr - sec_buff;
    //cout << i << endl;		

    if( t_gap < 0 ) { fail_counter += 1;// cout << "BIG jump backward at: sec previous: " << sec_buff << " sec current: " << sec_curr << endl;
    }
    if( t_gap > 0 && t_gap < 10 ) {sec_counter += t_gap; }
    if( t_gap > 10 && i != 0){ fail_counter += 1; 
	//cout << "BIG jump forward at: sec previous: " << sec_buff << " sec current: " << sec_curr << endl;
    }
    
    sec_buff = sec_curr;  

  }

    
/*
  slope = GetSlope( rootfile, place, 1);
  offset = GetOffset( rootfile, place, 1);
  //cout << slope << " " << offset << endl;
  Int_t nevent = tree->GetEntries();

  for( Int_t i = 0; i < nevent; i++){

      adcEvent->GetEvent(i);
      bin_ch = evt_r.padc[0];
      bin_ev = offset + bin_ch * slope;

      if(bin_ev > 4300 && bin_ev < 4720){ ti_counter += 1; }       
  
  }


 cout << ti_counter << endl;
 
 sec = ti_counter / ti_rate;
 */

 sec = (int)sec_counter/1;
 //cout << "The duration of this rootfile is approx " << sec << " seconds" << endl;
 day  = (Int_t)sec/(86400);
 hour = (Int_t)(sec - (day * 86400))/(3600);
 min  = (Int_t)(sec - (day * 86400 + hour * 3600))/(60);

 cout << "That is about " <<  day << " days, " << hour << " hours, " << min << " minutes" << endl;
 cout << "Also there were around " << fail_counter << " events with no ms clock timing information or big gaps between events" << endl;

 sec = (Int_t)sec;
 
 f->Close();
 delete f;
 return sec;

}

TH1F* FillHistoWithCondition(TString rootfile, TString place, TString cond_module, Int_t cond_channel, Int_t cond_low_limit, Int_t cond_up_limit, Int_t hist_channel){

// something should be added to also be able to plot not-qdc histograms
// this function returns a histogram of the (qdc) - channel "hist_channel", filled only with events which are between the low_limit and the up_limit in the channel "cond_channel"; the module (tdc, qdc, adc), from which the filling condition is taken, is given by "cond_module" 

  EventStruct  evt_r;
//  Int_t adc_channel_list[6] = {0,2,3,5,6,7};
//  Int_t adc_channel = adc_channel_list[sdd-1];

  Int_t size = place.Sizeof();
  TString rootfilename;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  const char *adc_str = "adc";
  const char *tdc_str = "tdc";
  const char *qdc_str = "qdc";

  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TH1F *hist = new TH1F("hist","", 1000, 0, 4000);

  TBranch *qdcEvent   = tree->GetBranch("qdc");
  qdcEvent->SetAddress(evt_r.qdc);
  TBranch *adcEvent   = tree->GetBranch("adc");
  adcEvent->SetAddress(evt_r.padc);
  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);

  Int_t nevent = tree->GetEntries();

  if(cond_module.EqualTo(tdc_str)){

    for( Int_t i = 0; i < nevent; i++){

      tdcEvent->GetEvent(i);
      qdcEvent->GetEvent(i);
      adcEvent->GetEvent(i);

        if(evt_r.tdc[cond_channel] > cond_low_limit && evt_r.tdc[cond_channel] < cond_up_limit){

          hist->Fill(evt_r.qdc[hist_channel]);
     
        }
 
    } 

  }

  if(cond_module.EqualTo(qdc_str)){

    for( Int_t i = 0; i < nevent; i++){

      tdcEvent->GetEvent(i);
      qdcEvent->GetEvent(i);
      adcEvent->GetEvent(i);

        if(evt_r.qdc[cond_channel] > cond_low_limit && evt_r.qdc[cond_channel] < cond_up_limit){

          hist->Fill(evt_r.qdc[hist_channel]);
     
        }
 
    } 

  }


  if(cond_module.EqualTo(adc_str)){

    Int_t sdd = PadcToSDD[cond_channel];

    Double_t slope = GetSlope(rootfile,place,sdd);
    Double_t offset = GetOffset(rootfile,place,sdd);
    Int_t bin_ev;

    for( Int_t i = 0; i < nevent; i++){

      tdcEvent->GetEvent(i);
      qdcEvent->GetEvent(i);
      adcEvent->GetEvent(i);

      bin_ev = offset + evt_r.padc[cond_channel] * slope;


        if(bin_ev > cond_low_limit && bin_ev < cond_up_limit){

          hist->Fill(evt_r.qdc[hist_channel]);
     
        }
 
    } 

  } 
  
//  hist->Draw();

  return hist;

}

TH1F* FillQdcHistoAdcCut(TString rootfile, TString place, Int_t qdc_channel, Int_t *adc_back, Int_t *scinti_coinc, Int_t *scinti_single, Int_t with_tdc_cut){

// this function fills a qdc histogram of "qdc_channel" with all events which are in the sdd background (any sdd > 7000 eV currently)
  EventStruct  evt_r;

  Int_t size = place.Sizeof();
  TString rootfilename;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  const char *adc_str = "adc";
  const char *tdc_str = "tdc";
  const char *qdc_str = "qdc";

  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TH1F *hist = new TH1F("hist","", 1000, 0, 4000);

  TBranch *qdcEvent   = tree->GetBranch("qdc");
  qdcEvent->SetAddress(evt_r.qdc);
  TBranch *adcEvent   = tree->GetBranch("adc");
  adcEvent->SetAddress(evt_r.padc);
  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);

  Int_t nevent = tree->GetEntries();
  Double_t slope[7];
  Double_t offset[7];

  Int_t adc_channel;
  Int_t tdc_channel = Qdc2Tdc_smi[qdc_channel];

  for(Int_t i = 0; i < 6; i++){

    slope[i+1] = GetSlope(rootfile,place,i+1);
    offset[i+1] = GetOffset(rootfile,place,i+1);

  }

  Int_t bin_ev[7];
  Int_t back_test;


  if(with_tdc_cut == 0){

    for( Int_t i = 0; i < nevent; i++){

      tdcEvent->GetEvent(i);
      qdcEvent->GetEvent(i);
      adcEvent->GetEvent(i);

      back_test = 0;


      for(Int_t j = 0; j < 6; j++){ 
 
         adc_channel = SDDToPadc[j+1];
         bin_ev[j+1] = offset[j+1] + evt_r.padc[adc_channel] * slope[j+1];
         if(bin_ev[j+1] > 7000){back_test = 1;}  // here is the condition for filling
         
      }

      if(back_test == 1){

          hist->Fill(evt_r.qdc[qdc_channel]);
          
          if(evt_r.qdc[qdc_channel] > QdcCuts[qdc_channel]) {*adc_back+=1;}
          if(evt_r.qdc[qdc_channel] > QdcCuts[qdc_channel] && evt_r.tdc[tdc_channel] > 5000) {*scinti_single+=1;}
          if(evt_r.qdc[qdc_channel] > QdcCuts[qdc_channel] && evt_r.tdc[13] > 5000) {*scinti_coinc+=1;}
     
      }

 
     }

  }

  if(with_tdc_cut == 1){

    for( Int_t i = 0; i < nevent; i++){

      tdcEvent->GetEvent(i);
      qdcEvent->GetEvent(i);
      adcEvent->GetEvent(i);

      back_test = 0;


      for(Int_t j = 0; j < 6; j++){ 
 
         adc_channel = SDDToPadc[j+1];
         bin_ev[j+1] = offset[j+1] + evt_r.padc[adc_channel] * slope[j+1];
         if(bin_ev[j+1] > 7000){back_test = 1;}  // here is the condition for filling
         
      }

      if(back_test == 1){

          if(evt_r.tdc[13]>5000){hist->Fill(evt_r.qdc[qdc_channel]);}
     
      }

 
     }

  }

//  hist->Draw();

  return hist;


}


TH1F* FillScaledHistogram(TString rootfile, TString place, Int_t sdd, Int_t low_edge_ev, Int_t ch_number){

// fills a histogram scaled in energy in a histogram starting from low_edge_ev... and from there ch_number adc channels upwards (1 ch ~ 9 eV)

  TH1F *histo;
  histo = HistoFromTree(rootfile, place, sdd); // only used for getting the slope and offset
  histo->Draw();
  //cout << 4 << endl;
  histo -> GetXaxis()->SetRangeUser(200,3000);
  

  TSpectrum *s = new TSpectrum();
  s->Search(histo, 2, "nodraw", 0.05 );  // 2. attribute is sigma                    

  Int_t peakN = s->GetNPeaks();

  //cout << peakN << endl;

  Double_t *pX = s->GetPositionX();
  //cout << px[0] << " " << pX[1] << endl;
     
  Int_t idx[peakN]; // thats an index array
  TMath::Sort(peakN,pX,idx,0); // idx now contains indices of sorted array
  cout << pX[idx[0]] << " " << pX[idx[1]] <<  " " << pX[idx[2]] << endl;
 
  Double_t slope = (MnKa1-TiKa1)/(pX[idx[2]]-pX[idx[0]]);
  cout << slope << endl;
  Double_t offset = TiKa1 - slope * pX[idx[0]];
  cout << offset << endl;

  EventStruct  evt_r;
  Int_t adc_channel_list[6] = {0,2,3,5,6,7};
  Int_t adc_channel = adc_channel_list[sdd-1];

  Int_t size = place.Sizeof();
  TString rootfilename;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;


  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TBranch *adcEvent;

  if(size==4)adcEvent = tree->GetBranch("adc");
  if(size==5)adcEvent = tree->GetBranch("adc");

  adcEvent->SetAddress(evt_r.padc);

  Int_t nevent = tree->GetEntries();

  //Double_t ch_number = 300; // amount of channels from the original adc histogram to include in the energy calibrated one
  //Double_t low_edge_ev = 7000;
  Double_t high_edge_ev = low_edge_ev + ch_number * slope;

  //Int_t low_edge_ev_int = (Int_t)low_edge_ev;
  Int_t high_edge_ev_int = (Int_t)high_edge_ev;
    

  TH1F *spect = new TH1F("spect","Energy spectrum", ch_number, low_edge_ev, high_edge_ev_int); 

  //TH1F *spect = new TH1F("spect","Energy spectrum", 100, 7000, 10200);
  Double_t bin_ev;

  for( Int_t i = 0; i < nevent; i++){

 
    adcEvent->GetEvent(i);  // here might be a problem if the address of the tbranch does not have the same size as the event in the branch

    bin_ev = offset + slope * evt_r.padc[adc_channel];
 
   // cout << bin_ev << " " << offset << " " << evt_r.padc[adc_channel] << " " << slope << endl;

    if(evt_r.padc[adc_channel]>0) {spect->Fill(bin_ev);}   
 
  }
    
  //spect->Draw();

  return spect;
}

TH1F* MakeQdcCutoffHisto(TString place, Int_t qdc_channel, Int_t fill_number){

   Int_t size = place.Sizeof();
   TH1F *hist = new TH1F("hist","",1000,0,4000);

   if(size == 4){

     for(Int_t i = 0; i < fill_number; i++){

        hist->Fill(QdcCuts[qdc_channel]); 

     } 

   }

   if(size == 5){

     for(Int_t i = 0; i < fill_number; i++){

        hist->Fill(QdcCuts_100nsGate[qdc_channel]); 

     }

  }  



  return hist;

} 



/*
void HistoFromTree(TString rootfile){

  EventStruct  evt_r; 

  TString rootfilename = ROOT_PATH_LNGS + "/" + rootfile;
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");
 
  TBranch *adcEvent   = tree->GetBranch("adc");
  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);
  adcEvent->SetAddress(evt_r.padc);

  Int_t nevent = tree->GetEntries();

  Int_t nevent_part = (Int_t)nevent/6;

  Int_t lbound = 0 * nevent_part;
  Int_t ubound = 1 * nevent_part;

  cout << nevent_part << endl;

  TH1F *spect = new TH1F("spect","", 1000, 0, 5000);
  TH1F *back = new TH1F("back","", 1000, 0, 5000);

  for( Int_t i = 0; i < nevent; i++){

    adcEvent->GetEvent(i);
    tdcEvent->GetEvent(i);

    spect->Fill(evt_r.padc[0]);

    if(evt_r.tdc[13]>6000){back->Fill(evt_r.padc[0]);}

 
  }

  TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
  cc->cd();

  gPad->SetLogy();
  back->SetLineColor(2);
 
//  htime->GetXaxis()->SetTitle("SDD Timing [ns * 10]");
//  htime->GetYaxis()->SetTitle("Scintillator Timing [ns * 10]");

  spect->Draw();
  back->Draw("same");

}
*/

void MakeOneQDCPlot(TString rootfile, TString place)
{
    Int_t size = place.Sizeof();

    TString rootfilename;

    if(size == 4){rootfilename = ROOT_PATH_SMI + "/" + rootfile;}
    if(size == 5){rootfilename = ROOT_PATH_LNGS + "/" + rootfile;}
    

    TFile *f = new TFile(rootfilename, "READ");
    f->cd();

    TH1F* hvq;
   // Int_t list [4] = {15,16,19,21};

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    cc->Divide(4,8);
    cc->cd();
    
    cc->SetGridx();
    cc->SetGridy();

    for( Int_t i = 0 ; i < NScinti ; i++){

       cc->cd(i+1);
       gPad->SetLogy();
       hvq = (TH1F*)f->Get( Form("hq[%d]", i ) );
       hvq->GetXaxis()->SetTitle( Form("QDC [%d]", i) );
       
       hvq->Draw();    
       
    }

}


void MakeOnePADCPlot(TString rootfile, TString place){

    Int_t size;
    size = place.Sizeof();
    TString rootfilename;
    if(size==4){ rootfilename = ROOT_PATH_SMI + "/" + rootfile;}
    if(size==5){ rootfilename = ROOT_PATH_LNGS + "/" + rootfile;}
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
       hsp -> GetXaxis()->SetRangeUser(0,4500);
       hsp->GetXaxis()->SetTitle( Form("PADC [%d]" , i));
       
       hsp->Draw();    
       
    }

}


void MakeOnePADCPlot2Files(TString rootfile1, TString rootfile2)
{ 
   
   Double_t scaling_factor_lngs = 343320;
   Double_t scaling_factor_smi = 6420; // 343320 sec for the ~ 3 day period end of october at smi; 342000 for the withoutcurrent files from july at smi (wild guess); 6420 for "highrate" with the other source
 

   Double_t scaling_factor_nogauge = 0.604; // 52200 sec for no gauge data 
   Double_t scaling_factor_gauge = 40.25;   // 40.25 d for BackRej_20160217_2.root at LNGS
   
   TH1F* hsp1;
   TH1F* hsp2;

   TString rootfilename1 = ROOT_PATH_LNGS + "/" + rootfile1;
   TString rootfilename2 = ROOT_PATH_LNGS + "/" + rootfile2;

   TFile *f1 = new TFile(rootfilename1, "READ");
   TFile *f2 = new TFile(rootfilename2, "READ");

   

   TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
   cc->cd();
   gPad->SetLogy();

   hsp1 = (TH1F*)f1->Get("rate"); // ------------------
   hsp2 = (TH1F*)f2->Get("rate");

   TH1F *gauge = new TH1F("with gauge","With Gauge",1000,0,10000);

   hsp2 -> GetXaxis()->SetRangeUser(-20,20);
   hsp1 -> GetXaxis()->SetRangeUser(-20,20);

  // Double_t counts_one = hsp1->GetEntries();
  // Double_t counts_two = hsp2->GetEntries();

  // scaling_factor = counts_one/counts_two;

   hsp2->SetLineColor(2);
   
   hsp1->Scale(1./scaling_factor_gauge);
   hsp2->Scale(1./scaling_factor_nogauge);

//   back_red->Divide(hsp1,hsp2,1,1);

  // hsp1->Fit("gaus","","",750,810);
  // hsp1->GetFunction("gaus")->SetLineColor(4);

  // hsp2->Fit("gaus","","",750,810);
  // hsp2->GetFunction("gaus")->SetLineColor(2);

   hsp1->Rebin(4);
   hsp2->Rebin(4);

   gStyle->SetOptStat(0);
   TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
   leg->SetHeader("rate with/without pressure gauge");
   leg->AddEntry(hsp1,"with pressure gauge","l");
   leg->AddEntry(hsp2,"no pressure gauge","l");
   

   hsp2->Draw();
   hsp1->Draw("same");
   leg->Draw("same");


   
}

/*
void MakeOneQDCPlot3Files(Int_t qdc_ch, TString rootfile1, TString rootfile2, TString rootfile3)
{ 
   
   //Double_t scaling_factor_lngs = 343320;
   //Double_t scaling_factor_smi = 6420; // 343320 sec for the ~ 3 day period end of october at smi; 342000 for the withoutcurrent files from july at smi (wild guess); 6420 for "highrate" with the other source
   
   TH1F* hsq1;
   TH1F* hsq2;
   TH1F* hsq3;

   TString rootfilename1 = ROOT_PATH_LNGS + "/" + rootfile1;
   TString rootfilename2 = ROOT_PATH_LNGS + "/" + rootfile2;
   TString rootfilename3 = ROOT_PATH_LNGS + "/" + rootfile3;

   TFile *f1 = new TFile(rootfilename1, "READ");
   TFile *f2 = new TFile(rootfilename2, "READ");
   TFile *f3 = new TFile(rootfilename3, "READ");
   

   TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
   cc->cd();
   gPad->SetLogy();

   hsq1 = (TH1F*)f1->Get(Form("hq[%d]",qdc_ch)); // ------------------
   hsq2 = (TH1F*)f2->Get(Form("hq[%d]",qdc_ch));
   hsq3 = (TH1F*)f3->Get(Form("hq[%d]",qdc_ch));


   hsq2 -> GetXaxis()->SetRangeUser(0,1000);
   hsq1 -> GetXaxis()->SetRangeUser(0,1000);
   hsq3 -> GetXaxis()->SetRangeUser(0,1000);

  // Double_t counts_one = hsp1->GetEntries();
  // Double_t counts_two = hsp2->GetEntries();

  // scaling_factor = counts_one/counts_two;

   hsq2->SetLineColor(2);
   hsq3->SetLineColor(3);

//   hsp1->Scale(1./scaling_factor_lngs);
//   hsp2->Scale(1./scaling_factor_smi);

   hsq2->Draw();
   hsq1->Draw("same");
   hsq3->Draw("same");
   gStyle->SetOptStat(0);

   TFitResultPtr *ptr;


   TFitResultPtr ptr_tmp = hsq3->Fit("gaus","S");
   Double_t mean_tmp = ptr_tmp->Parameter(1);
   Double_t sigma_tmp = ptr_tmp->Parameter(2);
 
 
   TFitResultPtr r = hsq3->Fit("gaus","S","",mean_tmp-4*sigma_tmp,mean_tmp+3*sigma_tmp);
   Double_t mean = r->Parameter(1);
   Double_t sigma = r->Parameter(2);

   Double_t cutoff_tmp = mean_tmp + 4* sigma_tmp;  
   Double_t cutoff =   mean + 4 * sigma;
   cout << "QDC number: " << qdc_ch << " old cutoff = " << cutoff_tmp << " new cutoff = " << cutoff << endl;

   leg = new TLegend(0.5,0.7,0.9,0.9);
   leg->SetHeader("Different QDC Gate width");
   leg->AddEntry(hsq1,"Gate width = 150 ns","l");
   leg->AddEntry(hsq2,"Gate width = 100 ns","l");
   leg->AddEntry(hsq3,"Gate width = 80 ns","l");
   leg->Draw();

   cc -> Print(Form("/home/andreas/vip2/reports/1601_ScintiAndQdcChanges/QdcGateWidth_Pictures/Scintillator_%d.pdf",qdc_ch));
   cc -> Delete();
  
}


void SaveQdcGatePictures(TString rootfile1, TString rootfile2, TString rootfile3){
 
    for( Int_t i=0; i<NScinti; i++ ){
      
      MakeOneQDCPlot3Files(i, rootfile1, rootfile2, rootfile3);

    }

}
*/

void MakeOneTDCPlot(TString rootfile, TString place)
{

    Int_t size = place.Sizeof();
    TString rootfilename;
    
    if(size==4){rootfilename = ROOT_PATH_SMI + "/" + rootfile;}
    if(size==5){rootfilename = ROOT_PATH_LNGS + "/" + rootfile;}
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


void MakeOneQDCPlot_DiffTrig(TString rootfile){

  EventStruct  evt_r; 

  TString rootfilename = ROOT_PATH_LNGS + "/" + rootfile;
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TH1F *h_self = new TH1F("h_self","",1000,0,4000);
  TH1F *h_other = new TH1F("h_other","",1000,0,4000);
  TH1F *h_all = new TH1F("h_all","",1000,0,4000);

  TBranch *qdcEvent   = tree->GetBranch("qdc");
  TBranch *tdcEvent   = tree->GetBranch("tdc");

  qdcEvent->SetAddress(evt_r.qdc);
  tdcEvent->SetAddress(evt_r.tdc);

  Int_t nevent = tree->GetEntries();

  for(Int_t i = 0; i < nevent; i++){

    qdcEvent->GetEvent(i);
    tdcEvent->GetEvent(i);

    if ( evt_r.tdc[6] > 2000){ h_self -> Fill(evt_r.qdc[4]);}
    if ( evt_r.tdc[6] < 2000){ h_other -> Fill(evt_r.qdc[4]);}
    h_all -> Fill(evt_r.qdc[4]);

  }

  TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
  cc->cd();
 
 // gPad->SetLogy();

  h_other->SetLineColor(2);
 // h_all->SetLineColor(1);

  //h_all->Draw();
  //h_self->GetYaxis()->SetRangeUser(0,18000);
  
  h_other->Draw();
  
 // h_self->Draw("same");
  
  
  
  

 
//  

  

}


void GetMnPeakPos(){

   gROOT->Reset();

   Int_t MnKa = 5898;
   Int_t TiKa = 4511;
   Double_t fwhm;
   Double_t fwhm_err;
   Double_t slope;

   TH1F* htmp;
  
   Double_t thr = 0.05;
  // int NumberOfFiles = 9;
   //Int_t day_list [9] = {15,16,17,18,23,24,25,29,30};
   const char *name_list[8];
   name_list[0] = "20151005_1.root";
   name_list[1] = "20151006_1.root";
   name_list[2] = "20151006_2.root";
   name_list[3] = "20151007_1.root";
   name_list[4] = "20151007_2.root";
   name_list[5] = "20151008_1.root";
   name_list[6] = "20151008_2.root";
   name_list[7] = "20151008_3.root";

   Int_t padc_channels[6] = {0,2,3,5,6,7}; 


for(int j = 1; j < 7; j++){ 

   ofstream outfile;
   outfile.open(Form("/home/andreas/vip2/reports/0915_SetupTest/MnPeakPos_sdd%d.txt",j));

   for(int i = 0; i < 8; i++){

 
     cout << "start with file: " << name_list[i] << endl << endl;
      
     // TString rootfile = Form("06%d_without.root",day_list[i]);

      TString rootfile = name_list[i];
      TString rootfilename = ROOT_PATH + "/" + rootfile;

      TFile *f1 = new TFile(rootfilename, "READ");
      
      htmp = (TH1F*)f1->Get(Form("hp[%d]",padc_channels[j-1]));

      htmp -> GetXaxis()->SetRangeUser(400,1000);

   //   htmp -> Draw();

      TSpectrum *s = new TSpectrum();
      s->Search(htmp, 2, "", thr );                

      

      Int_t peakN = s->GetNPeaks();

      //cout << peakN << endl;

      Double_t *pX = s->GetPositionX();
      
     
     Int_t idx[peakN];
     TMath::Sort(peakN,pX,idx,0);
  //   cout << pX[idx[0]] << " " << pX[idx[1]] <<  " " << pX[idx[2]] << endl;
 
     slope = (MnKa-TiKa)/(pX[idx[2]]-pX[idx[0]]);
   //  cout << slope << endl;
     
  // Double_t *pY = s->GetPositionY();

      TFitResultPtr r = htmp->Fit("gaus","S","",pX[idx[2]]-40,pX[idx[2]]+40);
      Double_t mean = r->Parameter(1);
      Double_t sigma = r->Parameter(2);
      fwhm_err = r->ParError(2);
      fwhm_err = fwhm_err * slope * 2.355;
      fwhm = sigma * 2.355 * slope;

      outfile << pX[idx[2]] << " " << mean << " " << fwhm << " +/- " << fwhm_err << endl;
   
      cout << "Number of peaks: " << peakN << endl;
      cout << "Mn peak at: " << pX[idx[2]] << " " << mean << endl << endl;
  
   }

   outfile.close();

}
   
}

void DrawQDCHitPattern(TString rootfile, Int_t Scinti){

  gROOT->Reset();

  EventStruct  evt_r; 

  TString rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TH2F *hpos = new TH2F("hpos",Form("Hit pattern for Scintillator Layer: %d, Column: %d",Qdc2Layer[Scinti],Qdc2Column[Scinti]), 8, -0.5, 7.5, 6, -0.5, 5.5);

  TBranch *qdcEvent   = tree->GetBranch("qdc");
  TBranch *adcEvent   = tree->GetBranch("adc");

  qdcEvent->SetAddress(evt_r.qdc);
  adcEvent->SetAddress(evt_r.padc);

  Int_t nevent = tree->GetEntries();

 // cout << nevent << endl;

  for( Int_t i = 0; i < nevent; i++){

    qdcEvent->GetEvent(i);
    adcEvent->GetEvent(i);

   // cout << evt_r.qdc[Scinti] << " " << QdcCuts[NScinti] << endl;

    if( evt_r.qdc[Scinti] > QdcCuts[Scinti]){ 

      for (Int_t j = 0; j < NScinti; j++){

        if( evt_r.qdc[j] > QdcCuts[j]){

          evt_r.layer = Qdc2Layer[j];
          evt_r.col   = Qdc2Column[j];
          hpos->Fill( evt_r.col, evt_r.layer );

          //cout << evt_r.layer << " " << evt_r.col << endl;
        }

      } 

    }

  }

  //hpos -> Draw("BOX,TEXT");
   TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
   cc->cd();
   hpos -> Draw("BOX,TEXT");

   cc -> Print(Form("/home/andreas/vip2/reports/1215_BackgroundAnalysis/QDCHitPattern%d.pdf",Scinti));
  

}




void DrawTDCHitPattern(TString rootfile, Int_t tdc_ch){

  gROOT->Reset();

  EventStruct  evt_r; 

  TString rootfilename = ROOT_PATH_LNGS + "/" + rootfile;
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

//  TH2F *hpos = new TH2F("hpos",Form("Hit pattern for Scintillator Layer: %d, Column: %d",Tdc2Layer_lngs[tdc_ch],Tdc2Column_lngs[tdc_ch]), 5, -0.5, 4.5, 5, -0.5, 4.5);
  TH2F *hpos = new TH2F("hpos","TDC hit pattern", 5, -0.5, 4.5, 5, -0.5, 4.5);
  TBranch *tdcEvent   = tree->GetBranch("tdc");

  tdcEvent->SetAddress(evt_r.tdc);

  Int_t nevent = tree->GetEntries();

  Int_t tdc_channel;

  for( Int_t i = 0; i < nevent; i++){

    tdcEvent->GetEvent(i);

/*     
    if( evt_r.tdc[tdc_ch] > 500){ 

      for (Int_t j = 0; j < 11; j++){

        if( evt_r.tdc[j] > 500){

          evt_r.layer = Tdc2Layer_lngs[j];
          evt_r.col   = Tdc2Column_lngs[j];
          hpos->Fill( evt_r.col, evt_r.layer );

          //cout << evt_r.layer << " " << evt_r.col << endl;
        }

      } 

    }

*/

    for(Int_t j = 0; j < 11; j++){// draw the complete hit distribution
  
      if(evt_r.tdc[j]>500){ 

        if(j == 3) {continue;}

        evt_r.layer = Tdc2Layer_lngs[j];
        evt_r.col = Tdc2Column_lngs[j];
        hpos->Fill( evt_r.col, evt_r.layer );
      }

    }
 }

  //hpos -> Draw("BOX,TEXT");
   TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
   cc->cd();
   hpos -> Draw("BOX,TEXT");

   cc -> Print(Form("/home/andreas/vip2/reports/1602_BackgroundAnalysis2/TDC_Hit_Pattern/tdc_20160421.pdf"));
  

}

void SaveHitPatternPictures(TString rootfile){
 
    for( Int_t i=0; i<NScinti; i++ ){
      
      DrawTDCHitPattern(rootfile,i);

    }

}

void DrawRejectionHisto(TString rootfile, TString place, Int_t sdd){

  Double_t slope = GetSlope(rootfile, place, sdd);
  Double_t offset = GetOffset(rootfile, place, sdd);
  Double_t bin_ev;
  Int_t tdc_channel;
  Int_t adc_channel = SDDToPadc[sdd];

  EventStruct  evt_r; 

  Int_t test = 0;

  Int_t size = place.Sizeof();
  TString rootfilename;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TBranch *adcEvent   = tree->GetBranch("adc");
  adcEvent->SetAddress(evt_r.padc);
  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);

  Int_t nevent = tree->GetEntries();

  TH1F *spect = new TH1F("spect","", 1000, 0, 50000);
  TH1F *spect_rej = new TH1F("spect_rej","", 1000, 0, 50000);
 

  for( Int_t i = 0; i < nevent; i++){

 
    adcEvent->GetEvent(i);
    tdcEvent->GetEvent(i);

    test = 0;

    for(Int_t j = 0; j < 8; j++){

      if(j==3){continue;}// change here

      tdc_channel = TdcChannels[j];
  
      if(evt_r.tdc[tdc_channel] > 500){ test = 1; break;}

    }

    bin_ev = offset + slope * evt_r.padc[adc_channel];

    if(bin_ev>1000){spect->Fill(bin_ev);}
  
   if(bin_ev>1000 && test == 1){spect_rej -> Fill(bin_ev);} 
 
  }

  TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
  cc->cd();

  gPad->SetLogy();
  gStyle->SetOptStat(0);

  spect->GetXaxis()->SetTitle("Energy [eV]");
  spect->GetYaxis()->SetTitle("Counts/50eV");
  spect->GetYaxis()->SetRangeUser(0.5,200000);




  spect->Draw();

  spect_rej->SetLineColor(2);

  spect_rej->Draw("same");
  TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
  leg->AddEntry(spect,"Complete SDD spectrum","l");
  leg->AddEntry(spect_rej,"Rejected by scintillator Veto","l");
  leg->Draw("same");
}


void AnalyseTDCHitPattern(TString rootfile, TString place){

  Int_t size = place.Sizeof();
  TString rootfilename;

  EventStruct  evt_r; 

  Int_t top_outer;
  Int_t top_inner;
  Int_t bottom_inner;
  Int_t bottom_outer;
  Int_t left_outer;
  Int_t left_inner;
  Int_t right_inner;
  Int_t right_outer;

  TH2F *det_eff = new TH2F("det_eff","detection efficiency from TDC", 5,-0.5,4.5,5,-0.5,4.5);


  if(size==4){rootfilename = ROOT_PATH_SMI + "/" + rootfile; top_outer = 4; top_inner = 5; bottom_inner = 1; bottom_outer = 7; left_outer = 0; left_inner = 3; right_inner = 2;  right_outer = 6;}
  if(size==5){rootfilename = ROOT_PATH_LNGS + "/" + rootfile; top_outer = 4; top_inner = 5; bottom_inner = 1; bottom_outer = 8; left_outer = 0; left_inner = 3; right_inner = 2; right_outer = 10;}

  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);

  Int_t tdc_layer_scinti[11]={0};
  Int_t tdc_layer[11]={0};
  Int_t tdc_layer_coinc[11]={0};

  Int_t main_counter[32] = {0};
  Int_t coinc_counter[32] = {0};


  Int_t tdc_channel;

  Int_t nevent = tree->GetEntries();

  for(Int_t i = 0; i < nevent; i++){

  tdcEvent->GetEvent(i);

    for(Int_t j = 0; j < 8; j++){ 

      if(size==4){tdc_channel = TdcConnections_smi[j];}
      if(size==5){tdc_channel = TdcConnections_lngs[j];}

      if(evt_r.tdc[tdc_channel]>500 && evt_r.tdc[tdc_channel]<7000){ tdc_layer_scinti[tdc_channel] = 1; } // layer tdc has scinti only timing
      if(evt_r.tdc[tdc_channel]>500){ tdc_layer[tdc_channel] = 1; } // layer tdc has a signal
      if(evt_r.tdc[tdc_channel]>7000){ tdc_layer_coinc[tdc_channel] = 1; } // layer tdc has coinc timing

       
    }


  if(tdc_layer_scinti[top_outer] == 1 && tdc_layer_scinti[bottom_inner] == 1){ main_counter[0] += 1; if(tdc_layer_scinti[top_inner] == 1) {coinc_counter[0] += 1;} }
  if(tdc_layer_scinti[top_inner] == 1 && tdc_layer_scinti[bottom_outer] == 1){ main_counter[1] += 1; if(tdc_layer_scinti[bottom_inner] == 1) {coinc_counter[1] += 1;} } 

  if(tdc_layer_scinti[left_inner] == 1 && tdc_layer_scinti[top_outer] == 1){ main_counter[2] += 1; if(tdc_layer_scinti[top_inner] == 1) {coinc_counter[2] += 1;} }
  if(tdc_layer_scinti[left_inner] == 1 && tdc_layer_scinti[bottom_outer] == 1){ main_counter[3] += 1; if(tdc_layer_scinti[bottom_inner] == 1) {coinc_counter[3] += 1;} } 

  if(tdc_layer_scinti[right_inner] == 1 && tdc_layer_scinti[top_outer] == 1){ main_counter[4] += 1; if(tdc_layer_scinti[top_inner] == 1) {coinc_counter[4] += 1;} }
  if(tdc_layer_scinti[right_inner] == 1 && tdc_layer_scinti[bottom_outer] == 1){ main_counter[5] += 1; if(tdc_layer_scinti[bottom_inner] == 1) {coinc_counter[5] += 1;} } 

  if(tdc_layer_scinti[right_inner] == 1 && tdc_layer_scinti[top_inner] == 1){ main_counter[6] += 1; if(tdc_layer_scinti[top_outer] == 1) {coinc_counter[6] += 1;} }
  if(tdc_layer_scinti[right_inner] == 1 && tdc_layer_scinti[bottom_inner] == 1){ main_counter[7] += 1; if(tdc_layer_scinti[bottom_outer] == 1) {coinc_counter[7]+= 1;} }

  if(tdc_layer_scinti[left_inner] == 1 && tdc_layer_scinti[top_inner] == 1){ main_counter[8] += 1; if(tdc_layer_scinti[top_outer] == 1) {coinc_counter[8] += 1;} }
  if(tdc_layer_scinti[left_inner] == 1 && tdc_layer_scinti[bottom_inner] == 1){ main_counter[9] += 1; if(tdc_layer_scinti[bottom_outer] == 1) {coinc_counter[9]+= 1;} }

  if(tdc_layer_scinti[bottom_outer] == 1 && tdc_layer_scinti[top_inner] == 1){ main_counter[10] += 1; if(tdc_layer_scinti[top_outer] == 1) {coinc_counter[10] += 1;} }
  if(tdc_layer_scinti[top_outer] == 1 && tdc_layer_scinti[bottom_inner] == 1){ main_counter[11] += 1; if(tdc_layer_scinti[bottom_outer] == 1) {coinc_counter[11]+= 1;} }

  if(tdc_layer_scinti[bottom_inner] == 1 && tdc_layer_scinti[right_outer] == 1){ main_counter[12] += 1; if(tdc_layer_scinti[right_inner] == 1) {coinc_counter[12]+= 1;} }
  if(tdc_layer_scinti[bottom_inner] == 1 && tdc_layer_scinti[left_outer] == 1){ main_counter[13] += 1; if(tdc_layer_scinti[left_inner] == 1) {coinc_counter[13]+= 1;} }



  for(Int_t k = 0; k < 11; k++){
      tdc_layer_scinti[k] = 0;
      tdc_layer_coinc[k] = 0;
      tdc_layer[k] = 0;
  }

 }

 det_eff -> SetBinContent(3,4,100*(float)coinc_counter[0]/(float)main_counter[0]);
 det_eff -> SetBinContent(3,2,100*(float)coinc_counter[1]/(float)main_counter[1]);
 det_eff -> SetBinContent(2,4,100*(float)coinc_counter[2]/(float)main_counter[2]);
 det_eff -> SetBinContent(2,2,100*(float)coinc_counter[3]/(float)main_counter[3]);
 det_eff -> SetBinContent(4,4,100*(float)coinc_counter[4]/(float)main_counter[4]);
 det_eff -> SetBinContent(4,2,100*(float)coinc_counter[5]/(float)main_counter[5]);
 det_eff -> SetBinContent(4,5,100*(float)coinc_counter[6]/(float)main_counter[6]);
 det_eff -> SetBinContent(4,1,100*(float)coinc_counter[7]/(float)main_counter[7]);
 det_eff -> SetBinContent(2,5,100*(float)coinc_counter[8]/(float)main_counter[8]);
 det_eff -> SetBinContent(2,1,100*(float)coinc_counter[9]/(float)main_counter[9]);
 det_eff -> SetBinContent(3,5,100*(float)coinc_counter[10]/(float)main_counter[10]);
 det_eff -> SetBinContent(3,1,100*(float)coinc_counter[11]/(float)main_counter[11]);
 det_eff -> SetBinContent(4,3,100*(float)coinc_counter[12]/(float)main_counter[12]);
// det_eff -> SetBinContent(2,3,100*(float)coinc_counter[13]/(float)main_counter[13]);

// cout << coinc_counter[11] << " " << main_counter[11] << endl;

 det_eff->Draw("BOX,TEXT");


/*
  cout << endl;
  cout << "Events with top outer and bottom inner signal is: " << counter1 << " and from these there are: " << counter2 << " events in top inner = " << 100*((float)counter2/(float)counter1) << " %" << endl;
  cout << "Events with top inner && bottom outer signal is: " << counter3 << " and from these there are: " << counter4 << " events in bottom inner = " << 100*((float)counter4/(float)counter3) << " %" << endl << endl;

  cout << "Events with left inner and top outer signal is: " << counter5 << " and from these there are: " << counter6 << " events in top inner = " << 100*((float)counter6/(float)counter5) << " %" << endl;
  cout << "Events with left inner && bottom outer signal is: " << counter7 << " and from these there are: " << counter8 << " events in bottom inner = " << 100*((float)counter8/(float)counter7) << " %" << endl << endl;

  cout << "Events with right inner && top outer signal is: " << counter9 << " and from these there are: " << counter10 << " events in top inner = " << 100*((float)counter10/(float)counter9) << " %" << endl;
  cout << "Events with right inner && bottom outer signal is: " << counter11 << " and from these there are: " << counter12 << " events in bottom inner = " << 100*((float)counter12/(float)counter11) << " %" << endl << endl;

   cout << "Events with right inner && top inner signal is: " << counter13 << " and from these there are: " << counter14 << " events in top inner = " << 100*((float)counter14/(float)counter13) << " %" << endl;
  cout << "Events with right inner && bottom inner signal is: " << counter15 << " and from these there are: " << counter16 << " events in bottom inner = " << 100*((float)counter16/(float)counter15) << " %" << endl << endl;

  cout << endl;
*/
}


Double_t GetRate(TString rootfile, Int_t lbound, Int_t ubound, TString place){

  // returns the rate averaged over all the rates of the events in the rootfile from event number lbound to event number ubound
  Int_t size = place.Sizeof();
  TString rootfilename;

  EventStruct  evt_r;

  Double_t rate_final;

  if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TBranch *rate_branch = tree->GetBranch("rate");
  rate_branch->SetAddress(&evt_r.rate);

  //TBranch *tdcEvent   = tree->GetBranch("tdc");
  //tdcEvent->SetAddress(evt_r.tdc);

  Double_t sum = 0;
  Int_t counter = 0;

  for(Int_t i = lbound; i < ubound; i++){

   
    rate_branch->GetEvent(i);
    if( evt_r.rate > 0 ) {sum +=  evt_r.rate; counter += 1;}

  }

  rate_final = sum/counter;

  return rate_final;
  
}
/*




void MakeOneTDC2dPlot(TString rootfile){

  EventStruct  evt_r; 

  TString rootfilename = ROOT_PATH_SMI + "/" + rootfile;
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");

  TH2F *htime = new TH2F("htime","Scintillator vs SDD timing",34,-1000,16000,34,-1000,16000);

  TBranch *tdcEvent   = tree->GetBranch("tdc");
  tdcEvent->SetAddress(evt_r.tdc);

  Int_t nevent = tree->GetEntries();

  for( Int_t i = 0; i < nevent; i++){

    tdcEvent->GetEvent(i);

    htime->Fill(evt_r.tdc[12],evt_r.tdc[13]);
 
  }

  TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
  cc->cd();

  gPad->SetLogz();
 
  htime->GetXaxis()->SetTitle("SDD Timing [ns * 10]");
  htime->GetYaxis()->SetTitle("Scintillator Timing [ns * 10]");

  htime->Draw("BOX");

}
*/

