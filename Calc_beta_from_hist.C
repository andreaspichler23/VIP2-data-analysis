{ 
gROOT->Reset();
//get histograms from root files:

    TString DATA_PATH = "/home/andreas/vip2/data/root/";

    f1 = new TFile (DATA_PATH+"withoutcurrent_sdd1.root","READ");
    TH1F *onewithout = (TH1F*)f1->Get("hev");

    f2 = new TFile (DATA_PATH+"withcurrent_sdd1.root","READ");
    TH1F *onewith = (TH1F*)f2->Get("hev");

// start commenting out here when only 1 sdd-----------------------------------

    f3 = new TFile (DATA_PATH+"withoutcurrent_sdd3.root","READ");
    TH1F *threewithout = (TH1F*)f3->Get("hev");

    f4 = new TFile (DATA_PATH+"withcurrent_sdd3.root","READ");
    TH1F *threewith = (TH1F*)f4->Get("hev");


    f5 = new TFile (DATA_PATH+"withoutcurrent_sdd4.root","READ");
    TH1F *fourwithout = (TH1F*)f5->Get("hev");

    f6 = new TFile (DATA_PATH+"withcurrent_sdd4.root","READ");
    TH1F *fourwith = (TH1F*)f6->Get("hev");


    f7 = new TFile (DATA_PATH+"withoutcurrent_sdd5.root","READ");
    TH1F *fivewithout = (TH1F*)f7->Get("hev");

    f8 = new TFile (DATA_PATH+"withcurrent_sdd5.root","READ");
    TH1F *fivewith = (TH1F*)f8->Get("hev");


    f9 = new TFile (DATA_PATH+"withoutcurrent_sdd6.root","READ");
    TH1F *sixwithout = (TH1F*)f9->Get("hev");

    f10 = new TFile (DATA_PATH+"withcurrent_sdd6.root","READ");
    TH1F *sixwith = (TH1F*)f10->Get("hev");

//    TH1F *sum_with = (TH1F*)onewith->Clone("sum_with");
//    TH1F *sum_without = (TH1F*)onewithout->Clone("sum_without");

    onewith->Add(threewith,1);
    onewith->Add(fourwith,1);
    onewith->Add(fivewith,1);
    onewith->Add(sixwith,1);

    onewithout->Add(threewithout,1);
    onewithout->Add(fourwithout,1);
    onewithout->Add(fivewithout,1);
    onewithout->Add(sixwithout,1);



// end comment out-----------------------------


    onewith->SetLineColor(2);

//normalise by their entries

    Double_t with_ent = onewith->GetEntries();
    Double_t without_ent = onewithout->GetEntries();
 
    Double_t scaler = with_ent/without_ent;

    onewithout->Scale(scaler);


// subtract without current from with current

   TH1F *substract = (TH1F*)onewith->Clone("substract");
   substract->SetLineColor(9);

   substract -> Add(onewithout,-1);
  //  onewithout->Draw();
  //  onewith->Draw("same");
   substract->Draw();

// get the integral of values in the region of interest

   Double_t lbound = 7580; // as the expected energy seems to be 7679 eV
   Double_t ubound = 7780;

   Double_t lbin = substract->FindBin(lbound);
   Double_t ubin = substract->FindBin(ubound);

   Double_t integral_sub = substract->Integral(lbin,ubin);
   Double_t integral_with = onewith->Integral(lbin,ubin);
   Double_t integral_without = onewithout->Integral(lbin,ubin);

   

// Calculate the errors

   Double_t error = 0;
//   Double_t error_with = 0;
//   Double_t error_without = 0;
/*
   for ( Int_t i = lbin ; i <= ubin ; i++){

      error_with = onewith->GetBinError(i);
      error_without = onewithout->GetBinError(i);

      error = error + TMath::Abs(error_with - error_without);

   }
*/ 

   error = TMath::Sqrt(integral_with+integral_without);
   cout << "Substracted counts: "<< integral_sub << " +/- " << error << endl;
   cout << "Counts with current: " << integral_with << " +/- " << TMath::Sqrt(integral_with) << endl;
   cout << "Counts without current: " << integral_without << " +/- " << TMath::Sqrt(integral_without) <<endl;
   cout << with_ent << endl;






}
