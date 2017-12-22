/**********************
ReadSlowData.C

16.10.2015
A. Pichler

**********************/



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


typedef struct{


    Int_t  year; 
    Int_t  month; 
    Int_t  day; 
    Int_t  hour; 
    Int_t  min;
    Int_t  sec;

    TDatime time_tdat;
    Double_t time_d;


    char AmOrPm;

    char ticks;

    Double_t room;
    Double_t ar_target;
    Double_t pcb_a;
    Double_t pcb_b;
    Double_t sdd1;
    Double_t argon_up;
    Double_t argon_down;
    Double_t sdd2;
    Double_t copper_int;
    Double_t copper_ext1;
    Double_t cooling_pad;

    Double_t vacuum_pres;
    Double_t current;
    Double_t argon_pres;
    Double_t voltage;

    Double_t heater_percent;
    Double_t ar_gas;

} SlowStruct;




int MakeSlowTree( TString slow_file, TString rootfilename){

    SlowStruct slow_s; 
    TString rootfile; 
    char  fline[512];
   

 
    FILE  *slowfile;
    if( (slowfile = fopen(SLOW_TEXT_PATH + "/" +slow_file, "r")) == NULL ){
        cout << "Cannot open file: " << slow_file << endl;
    }
    cout << slow_file << endl;

    rootfile = SLOW_ROOT_PATH + "/" + rootfilename; 
    TFile *f = new TFile( rootfile, "RECREATE");
    f->cd();

    TTree *tr = new TTree("tr", "Slow Data");

    tr->Branch("room",   &slow_s.room,   "Room temperature/D");
    tr->Branch("ar_target",   &slow_s.ar_target,   "Argon Target temperature/D");
    tr->Branch("pcb_a",  &slow_s.pcb_a,  "PCB A temperature/D");
    tr->Branch("pcb_b",  &slow_s.pcb_b,  "PCB B temperature/D");
    tr->Branch("sdd1",  &slow_s.sdd1,  "SDD 1 temperature/D");
    tr->Branch("argon_up",   &slow_s.argon_up,     "Argon Up temperature/D");
    tr->Branch("argon_down",   &slow_s.argon_down,     "Argon Down temperature/D");
    tr->Branch("sdd2",    &slow_s.sdd2,    "SDD 2 temperature/D");
    tr->Branch("copper_int",   &slow_s.copper_int,     "Internal Copper temperature/D");
    tr->Branch("copper_ext1",   &slow_s.copper_ext1,     "External Copper 1 temperature/D");
    tr->Branch("cooling_pad",    &slow_s.cooling_pad,     "Cooling Pad temperature/D");
    tr->Branch("vacuum_pres",   &slow_s.vacuum_pres,    "Vacuum pressure/D");
    tr->Branch("current",  &slow_s.current,   "Current/D");
    tr->Branch("argon_pres",    &slow_s.argon_pres,     "Argon pressure/D");
    tr->Branch("voltage",     &slow_s.voltage,      "Voltage/D");
    tr->Branch("heater_percent",       &slow_s.heater_percent,     "Heater percent/D");
    tr->Branch("ar_gas",   &slow_s.ar_gas,   "Argon Gas temperature/D");

    tr->Branch("time_b", &slow_s.time_tdat);
    tr->Branch("time_d", &slow_s.time_d);


    

    fgets( fline, 512, slowfile );// read the header


    while(fgets(fline,512,slowfile)){


     sscanf( fline, "%2d-%2d-%2d %2d:%2d:%2d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",

     &slow_s.day, &slow_s.month, &slow_s.year, &slow_s.hour, &slow_s.min, &slow_s.sec, &slow_s.ticks, &slow_s.room, &slow_s.ar_target, &slow_s.pcb_a, &slow_s.pcb_b, &slow_s.sdd1, &slow_s.argon_up, &slow_s.argon_down, &slow_s.sdd2, &slow_s.copper_int, &slow_s.copper_ext1, &slow_s.cooling_pad, &slow_s.vacuum_pres, &slow_s.current, &slow_s.argon_pres, &slow_s.voltage, &slow_s.heater_percent, &slow_s.ar_gas);

     slow_s.year = slow_s.year + 2000;
     //cout << slow_s.year << " " << slow_s.hour << " " << slow_s.min << " " << slow_s.month << " " << slow_s.day << endl;
     
     slow_s.time_tdat.Set(slow_s.year, slow_s.month, slow_s.day, slow_s.hour, slow_s.min, slow_s.sec);
     
     
    

//    if(slow_s.AmOrPm=='A' && slow_s.hour == 12){slow_s.hour = 0; slow_s.time_tdat.Set(slow_s.year, slow_s.month, slow_s.day, slow_s.hour, slow_s.min, slow_s.sec);}
//    if(slow_s.AmOrPm=='A' && slow_s.hour != 12){slow_s.time_tdat.Set(slow_s.year, slow_s.month, slow_s.day, slow_s.hour, slow_s.min, slow_s.sec);}
//    if(slow_s.AmOrPm=='P' && slow_s.hour == 12){slow_s.time_tdat.Set(slow_s.year, slow_s.month, slow_s.day, slow_s.hour, slow_s.min, slow_s.sec);}
//    if(slow_s.AmOrPm=='P' && slow_s.hour != 12){slow_s.hour = slow_s.hour + 12; slow_s.time_tdat.Set(slow_s.year, slow_s.month, slow_s.day, slow_s.hour, slow_s.min, slow_s.sec);}

   slow_s.time_d=slow_s.time_tdat.Convert();

 //  cout << slow_s.AmOrPm << endl;
   
    tr->Fill();

    }

   fclose(slowfile);
  // tr->Print();
   

   f->Write();
   f->Close();
  
   return 0;



}


void DrawSlowLogPictures( TString rootfile){

  TString rootfilename = SLOW_ROOT_PATH + "/" + rootfile;
  TFile *f = new TFile(rootfilename, "READ");
  TTree *tree = (TTree*)f->Get("tr");


  Double_t *room_t = new Double_t();
  TBranch *room = tree->GetBranch("room");
  room->SetAddress(room_t);

  Double_t *pcba_d = new Double_t();
  TBranch *pcb_a = tree->GetBranch("pcb_a");
  pcb_a->SetAddress(pcba_d);

  Double_t *pcbb_d = new Double_t();
  TBranch *pcb_b = tree->GetBranch("pcb_b");
  pcb_b->SetAddress(pcbb_d);

  Double_t *sdd1_d = new Double_t();
  TBranch *sdd1 = tree->GetBranch("sdd1");
  sdd1->SetAddress(sdd1_d);

  Double_t *sdd2_d = new Double_t();
  TBranch *sdd2 = tree->GetBranch("sdd2");
  sdd2->SetAddress(sdd2_d);

  Double_t *arpres_d = new Double_t();
  TBranch *argon_pres = tree->GetBranch("argon_pres");
  argon_pres->SetAddress(arpres_d);

  Double_t *artar_d = new Double_t();
  TBranch *ar_target = tree->GetBranch("ar_target");
  ar_target->SetAddress(artar_d);

  Double_t *argas_d = new Double_t();
  TBranch *ar_gas = tree->GetBranch("ar_gas");
  ar_gas->SetAddress(argas_d);

  Double_t *current_d = new Double_t();
  TBranch *current = tree->GetBranch("current");
  current->SetAddress(current_d);

  Double_t *argonup_d = new Double_t();
  TBranch *argon_up = tree->GetBranch("argon_up");
  argon_up->SetAddress(argonup_d);

  Double_t *argondown_d = new Double_t();
  TBranch *argon_down = tree->GetBranch("argon_down");
  argon_down->SetAddress(argondown_d);

  Double_t *coolingpad_d = new Double_t();
  TBranch *cooling_pad = tree->GetBranch("cooling_pad");
  cooling_pad->SetAddress(coolingpad_d);

  Double_t *copperint_d = new Double_t();
  TBranch *copper_int = tree->GetBranch("copper_int");
  copper_int->SetAddress(copperint_d);

  Double_t *heater_d = new Double_t();
  TBranch *heater_percent = tree->GetBranch("heater_percent");
  heater_percent->SetAddress(heater_d);


  Double_t *time_dou = new Double_t();
  TBranch *time_d = tree->GetBranch("time_d");
  time_d->SetAddress(time_dou);



  Int_t nentries = tree->GetEntries();
  Double_t room_temp[nentries];
  Double_t time_d_arr[nentries];
  Double_t current_d_arr[nentries];
  Double_t sdd1_d_arr[nentries];
  Double_t sdd2_d_arr[nentries];
  Double_t arpres_d_arr[nentries];
  Double_t pcba_d_arr[nentries];
  Double_t pcbb_d_arr[nentries];
  Double_t artar_d_arr[nentries];
  Double_t argas_d_arr[nentries];
  Double_t argonup_d_arr[nentries];
  Double_t argondown_d_arr[nentries];
  Double_t coolingpad_d_arr[nentries];
  Double_t copperint_d_arr[nentries];
  Double_t heater_d_arr[nentries];

//  cout << 1 << endl;
  for(Int_t i=0; i < nentries; i++){

    room->GetEvent(i);
    room_temp[i] = room_t[0];

    pcb_a->GetEvent(i);
    pcba_d_arr[i] = pcba_d[0];

    pcb_b->GetEvent(i);
    pcbb_d_arr[i] = pcbb_d[0];

    current->GetEvent(i);
    current_d_arr[i] = current_d[0];

    sdd1->GetEvent(i);
    sdd1_d_arr[i] = sdd1_d[0];

    sdd2->GetEvent(i);
    sdd2_d_arr[i] = sdd2_d[0];

    argon_pres->GetEvent(i);
    arpres_d_arr[i] = arpres_d[0];

    ar_target->GetEvent(i);
    artar_d_arr[i] = artar_d[0];

    ar_gas->GetEvent(i);
    argas_d_arr[i] = argas_d[0];

    argon_up->GetEvent(i);
    argonup_d_arr[i] = argonup_d[0];

    argon_down->GetEvent(i);
    argondown_d_arr[i] = argondown_d[0];

    cooling_pad->GetEvent(i);
    coolingpad_d_arr[i] = coolingpad_d[0];

    copper_int->GetEvent(i);
    copperint_d_arr[i] = copperint_d[0];

    heater_percent->GetEvent(i);
    heater_d_arr[i] = heater_d[0];

    time_d->GetEvent(i);
    time_d_arr[i] = time_dou[0];   

  // cout << room_temp[i] << endl;

  }

 // cout << 2 << endl;

  TGraph* slow_g[14];
 // cout << 3 << endl;


//  cout << 4 << endl;



  TDatime *start_t = new TDatime(2016,6,20,10,20,0); // start and ending time for the plots
  Double_t start_d = start_t->Convert();

  TDatime *end_t = new TDatime(2016,6,20,24,0,0);
  Double_t end_d = end_t->Convert();

 
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Divide(1,4);

  slow_g[0] = new TGraph(nentries, time_d_arr, sdd1_d_arr);
  slow_g[0]->GetYaxis()->SetTitle("SDD1");
  
  slow_g[1] = new TGraph(nentries, time_d_arr, current_d_arr);
  slow_g[1]->GetYaxis()->SetTitle("Current");

  slow_g[2] = new TGraph(nentries, time_d_arr, sdd2_d_arr);
  slow_g[2]->GetYaxis()->SetTitle("SDD2");
  
  slow_g[3] = new TGraph(nentries, time_d_arr, room_temp);
  slow_g[3]->GetYaxis()->SetTitle("Room Temperature");
  

  for(Int_t i = 0; i < 4; i ++){
   
   c1->cd(i+1);
   slow_g[i]->GetXaxis()->SetTimeDisplay(1);
   slow_g[i]->GetXaxis()->SetTimeFormat("%H:%M");
   slow_g[i]->GetXaxis()->SetTimeOffset(0);
   slow_g[i]->GetXaxis()->SetLabelSize(0.1);
   slow_g[i]->GetXaxis()->SetRangeUser(start_d,end_d);
   slow_g[i]->GetYaxis()->SetLabelSize(0.1);
   slow_g[i]->GetYaxis()->SetTitleOffset(0.5);
   slow_g[i]->GetYaxis()->SetTitleSize(0.1);
   slow_g[i]->Draw();

  }

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->Divide(1,4);

  slow_g[4] = new TGraph(nentries, time_d_arr, arpres_d_arr);
  slow_g[4]->GetYaxis()->SetTitle("Argon Pressure");
  
  slow_g[5] = new TGraph(nentries, time_d_arr, pcba_d_arr);
  slow_g[5]->GetYaxis()->SetTitle("PCB_A Temp");

  slow_g[6] = new TGraph(nentries, time_d_arr, pcbb_d_arr);
  slow_g[6]->GetYaxis()->SetTitle("PCB_B Temp");
  
  slow_g[7] = new TGraph(nentries, time_d_arr, heater_d_arr);
  slow_g[7]->GetYaxis()->SetTitle("Heater Percent");

  for(Int_t i = 4; i < 8; i ++){
   
   c2->cd(i-3);
   slow_g[i]->GetXaxis()->SetTimeDisplay(1);
   slow_g[i]->GetXaxis()->SetTimeFormat("%H:%M");
   slow_g[i]->GetXaxis()->SetTimeOffset(0);
   slow_g[i]->GetXaxis()->SetLabelSize(0.1);
   slow_g[i]->GetXaxis()->SetRangeUser(start_d,end_d);
   slow_g[i]->GetYaxis()->SetLabelSize(0.1);
   slow_g[i]->GetYaxis()->SetTitleOffset(0.5);
   slow_g[i]->GetYaxis()->SetTitleSize(0.1);
   slow_g[i]->Draw();

  }
 


/*
  Double_t room_temp[nentries];
  Double_t time_d_arr[nentries];
  Double_t current_d_arr[nentries];
  Double_t sdd1_d_arr[nentries];
  Double_t sdd2_d_arr[nentries];
  Double_t arpres_d_arr[nentries];
  Double_t pcba_d_arr[nentries];
  Double_t artar_d_arr[nentries];
  Double_t argonup_d_arr[nentries];
  Double_t argondown_d_arr[nentries];
  Double_t coolingpad_d_arr[nentries];
  Double_t copperint_d_arr[nentries];
  Double_t heater_d_arr[nentries];
*/
  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  c3->Divide(1,4);

  slow_g[8] = new TGraph(nentries, time_d_arr, argondown_d_arr);
  slow_g[8]->GetYaxis()->SetTitle("Argon Down Temp");
  
  slow_g[9] = new TGraph(nentries, time_d_arr, argonup_d_arr);
  slow_g[9]->GetYaxis()->SetTitle("Argon Up Temp");

  slow_g[10] = new TGraph(nentries, time_d_arr, artar_d_arr);
  slow_g[10]->GetYaxis()->SetTitle("Argon Target Temp");
  
  slow_g[11] = new TGraph(nentries, time_d_arr, argas_d_arr);
  slow_g[11]->GetYaxis()->SetTitle("Argon Gas Temp");

  for(Int_t i = 8; i < 12; i ++){
   
   c3->cd(i-7);
   slow_g[i]->GetXaxis()->SetTimeDisplay(1);
   slow_g[i]->GetXaxis()->SetTimeFormat("%H:%M");
   slow_g[i]->GetXaxis()->SetTimeOffset(0);
   slow_g[i]->GetXaxis()->SetLabelSize(0.1);
   slow_g[i]->GetXaxis()->SetRangeUser(start_d,end_d);
   slow_g[i]->GetYaxis()->SetLabelSize(0.1);
   slow_g[i]->GetYaxis()->SetTitleOffset(0.5);
   slow_g[i]->GetYaxis()->SetTitleSize(0.1);
   slow_g[i]->Draw();

  }

  TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
  c4->Divide(1,2);

  slow_g[12] = new TGraph(nentries, time_d_arr, coolingpad_d_arr);
  slow_g[12]->GetYaxis()->SetTitle("Cooling Pad Temp");

  slow_g[13] = new TGraph(nentries, time_d_arr, copperint_d_arr);
  slow_g[13]->GetYaxis()->SetTitle("Copper Internal");

  for(Int_t i = 12; i < 14; i ++){
   
   c4->cd(i-11);
   slow_g[i]->GetXaxis()->SetTimeDisplay(1);
   slow_g[i]->GetXaxis()->SetTimeFormat("%H:%M");
   slow_g[i]->GetXaxis()->SetTimeOffset(0);
   slow_g[i]->GetXaxis()->SetLabelSize(0.1);
   slow_g[i]->GetXaxis()->SetRangeUser(start_d,end_d);
   slow_g[i]->GetYaxis()->SetLabelSize(0.1);
   slow_g[i]->GetYaxis()->SetTitleOffset(0.5);
   slow_g[i]->GetYaxis()->SetTitleSize(0.1);
   slow_g[i]->Draw();

  }

}




