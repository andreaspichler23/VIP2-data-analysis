{

TFile *f = new TFile("20160217_2.root");
TTree *tree = (TTree*)f->Get("tr");
TH1F *htmp = (TH1F*)f->Get("hp[0]");

ofstream txt_file;
txt_file.open("VIP2_energy_histogram.txt");


//Int_t ent = tr->GetEntries();
Int_t y;

for(Int_t i = 0; i < 2047; i++){

 y = htmp->GetBinContent(i);

 txt_file << i << " " << y << endl;

}

txt_file.close();

}
