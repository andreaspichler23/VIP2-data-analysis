{

TFile *vip_results = new TFile("vip_results.root","RECREATE","VIP results");

//TH1F *with_current = new TH1F("with_current","Data with current",4096,0,24576);
//TH1F *without_current = new TH1F("without_current","Data without current",4096,0,24576);
TH1F *subtracted = new TH1F("subtracted","Subtracted Spectrum",4096,0,24576);



ifstream txt_file;
txt_file.open("histo_laura.txt");

Int_t with_current_bin;
Int_t without_current_bin;
Double_t subtracted_bin;
Int_t subtracted_bin_int;

for(Int_t i = 0; i < 4096; i++){


 txt_file >> with_current_bin;

 with_current -> SetBinContent(i,with_current_bin);

// cout << with_current_bin << endl;



 txt_file >> without_current_bin;

 without_current -> SetBinContent(i,without_current_bin);

// cout << without_current_bin << endl;


 txt_file >> subtracted_bin;

 subtracted_bin_int = (Int_t)subtracted_bin;

 subtracted -> SetBinContent(i,subtracted_bin_int);

// cout << subtracted_bin_int << endl;

 

}

with_current->GetXaxis()->SetTitle("Energy[eV]");
without_current->GetXaxis()->SetTitle("Energy[eV]");
subtracted->GetXaxis()->SetTitle("Energy[eV]");

with_current->GetYaxis()->SetTitle("Counts/6eV");
without_current->GetYaxis()->SetTitle("Counts/6eV");
subtracted->GetYaxis()->SetTitle("Counts/6eV");

vip_results->Write();

}
