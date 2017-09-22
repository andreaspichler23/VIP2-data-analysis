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

#include  "CalibFunction.C"
#include  "MyRootFunctionLib_Basic.C"

void fitGaussSubHist(TString histName, Int_t rebin){
    
    //fitGaussSubHist("SubHistSum",25)
    
    TFile *fin = new TFile("/home/andreas/vip2/data/root/LNGS/1-618files-final/energyHistograms.root");
    
    TH1F *subH = (TH1F*)fin->Get(histName);
    subH->Rebin(rebin);
    
    Int_t ll = 7500;
    Int_t ul = 7900;
    
    TF1 *fitFunc = new TF1("fitFunc", gaussFunc, ll, ul, 3 );
    
    fitFunc->FixParameter(1,7729);
    fitFunc->FixParameter(2,85);
    
    fitFunc->SetParameter(0,1000);;
    subH->GetXaxis()->SetRangeUser(ll,ul);
    subH->Fit("fitFunc","RI");
    
    Double_t mpv = fitFunc->GetParameter(0);   
    mpv = mpv / rebin;
    
    Double_t sigma = fitFunc->GetParError(0);
    sigma = TMath::Sqrt(sigma / rebin);
    
    cout << " Pauli violating counts: " << mpv << " +- " << sigma << endl;
    
    subH->Draw();
    
}

void FitHistogram(  TString histString, Int_t savePlot, Int_t rebin ){
    
    
    Double_t counts;
    Int_t lowerL = 7000;
    Int_t upperL = 9500;
    Int_t nPar = 11;
    
    Double_t fitResult[nPar+4];
    
    ofstream fitFile;
    fitFile.open("/home/andreas/vip2/reports/Analysis/CuFitFile.txt",ios::app);
    
    TFile *inF = TFile::Open("/home/andreas/vip2/data/root/LNGS/1-618files-final/energyHistograms.root");
    TH1F *hist = (TH1F*)inF->Get(histString);
    
    hist->Rebin(rebin);
    hist->GetXaxis()->SetRangeUser(lowerL,upperL);
    counts = hist->Integral();
    
    Double_t cuKa1GainStart = (counts / 6) * rebin;
    Double_t cuKbGainStart = (counts / 30) * rebin;
    Double_t niKa1GainStart = (counts / 100) * rebin;
    
//    cout << "counts: " << counts << endl;
//    cout << "cu ka1 gain start: " << cuKa1GainStart << endl;
//    cout << "cu kb gain start: " << cuKbGainStart << endl;
//    cout << "ni gain start: " << niKa1GainStart << endl;
    
    Double_t backStart = (counts / 4000) * rebin;
    Double_t backSlopeStart = -(1./2000) * rebin;
    
//    cout << "back slope start: " << backSlopeStart << endl;
    
    Double_t sigmaStart = 75;
    
    
    TF1 *fitFunc = new TF1("fitFunc", RoiCuFunc, lowerL, upperL, nPar );
    
    fitFunc->SetParName(0,"Background constant");
    fitFunc->SetParName(10,"Background Slope");
    
    fitFunc->SetParName(1,"Cu Ka1 Gain");
    fitFunc->SetParName(4,"Cu Kb Gain");
    fitFunc->SetParName(7,"Ni Ka1 Gain");
    
    fitFunc->SetParName(2,"Cu Ka1 Mean");
    fitFunc->SetParName(5,"Cu Kb Mean");
    fitFunc->SetParName(8,"Ni Ka1 Mean");
    
    fitFunc->SetParName(3,"Cu Ka1 Sigma");
    fitFunc->SetParName(6,"Cu Kb Sigma");
    fitFunc->SetParName(9,"Ni Ka1 Sigma");
    
    
    fitFunc->SetParameter(0,backStart);
    fitFunc->SetParLimits(0,0,backStart*2);
    fitFunc->SetParameter(10,backSlopeStart);
    
    fitFunc->SetParameter(1,cuKa1GainStart);
    fitFunc->SetParameter(4,cuKbGainStart);
    fitFunc->SetParameter(7,niKa1GainStart);
    //fitFunc->SetParLimits(1,0,cuKa1GainStart*2);
    
    fitFunc->SetParameter(2,CuKa1);
    fitFunc->SetParLimits(2,CuKa1-10,CuKa1+30);
    fitFunc->SetParameter(5,CuKb1);
    fitFunc->SetParLimits(5,CuKb1-10,CuKb1+30);
    fitFunc->SetParameter(8,NiKa1);
    fitFunc->SetParLimits(8,NiKa1-10,NiKa1+10);
    
    fitFunc->SetParameter(3,sigmaStart);
    fitFunc->SetParameter(6,sigmaStart);
    fitFunc->SetParameter(9,sigmaStart);
    
    fitFunc->SetParLimits(9,sigmaStart-20,sigmaStart+20);
    
    hist->Fit("fitFunc", "RI");
    
    fitFunc->GetParameters(fitResult);
    fitResult[nPar]     = fitFunc->GetChisquare();
    fitResult[nPar+1]   = fitFunc->GetNDF();
    fitResult[nPar+2]   = fitResult[nPar]/fitResult[nPar+1];
    fitResult[nPar+3]   = fitResult[3] * 2.355;
    
    
    
    //fitFile << histString << " ";
    
    cout << "fwhm is: " << fitResult[nPar+3] << " eV!" << endl;
    
    for( Int_t i = 0; i <= nPar+3; i++ ){
        
        //fitFile << fitResult[i] << " ";
        
        
    }
    fitFile << endl;
    //hist->Draw();
    
    
    
    // calculate counts in the roi 
    
//    Double_t roiLL = 7729 - fitResult[nPar+3] / 2;
//    Double_t roiUL = 7729 + fitResult[nPar+3] / 2;
//    
//    hist->Delete();
//    //delete hist;
//    
//    TH1F *hist2 = (TH1F*)inF->Get(histString); // again with 1 eV binning
//    
//    
//    TAxis *axis = hist2->GetXaxis();
//
//    Int_t bmin = axis->FindBin(roiLL);
//    Int_t bmax = axis->FindBin(roiUL);
//
//    Int_t countsROI = hist2->Integral(bmin,bmax);
//
//    fitFile << roiLL << " " << roiUL << " " << countsROI << " " << rebin << endl;
//
//    hist2->Delete();
    inF->Close();
    delete inF;
    fitFile.close();
 
    
}

void FitHistogramLoop(Int_t rebin){
    
    TString histString;
    
    for( Int_t i = 1; i <= 6; i++){
        
        
        histString = Form("withCurrentsdd%d",i);
        cout << endl << histString << " START" << endl;
        FitHistogram(histString,0,rebin);
        
        
    }
    
    histString = "withCurrentSum";
    cout << endl << histString << " START" << endl;
    FitHistogram(histString,0,rebin);
    
    for( Int_t i = 1; i <= 6; i++){
        
        
        histString = Form("noCurrentSmallsdd%d",i);
        cout << endl << histString << " START" << endl;
        FitHistogram(histString,0,rebin);
        
        
    }
    
    histString = "noCurrentSmallSum";
    cout << endl << histString << " START" << endl;
    FitHistogram(histString,0,rebin);

    
}

Double_t calcRoiContribution(Double_t mean, Double_t gain, Double_t sigma, Int_t roiLL, Int_t roiUL, Int_t rebin){
    
    //TF1 *gaussian = new TF1("gaussian","gaus(0)",7000,9000);
    

    TF1 *gaussian = new TF1("gaussian",gaussFunc,7000,9000,3);
    
    gaussian->SetParameter(0,gain);
    gaussian->SetParameter(1,mean);
    gaussian->SetParameter(2,sigma);
    
    Double_t cont = gaussian->Integral(roiLL,roiUL);
    
    cont = cont / rebin;
    
    //gaussian->Draw();
    
    return cont;
    
}

void calcROICounts(){
    
    // adjust runTime und histString
    
    TFile *f = new TFile("/home/andreas/vip2/data/root/LNGS/1-618files-final/energyHistograms.root");
    
    ifstream cuFitFile;
    cuFitFile.open("/home/andreas/vip2/reports/Analysis/CuFitFile.txt");
    Double_t calibList[210];
    
    TString withString = "with";
    TString noString = "no";
    TString histString;
    TString saveString;
    
    
    Double_t roiCountsWith[6], roiCountsNo[6];
    
    Int_t maxBinNumber = 9002;
    Double_t binCenter;
    Int_t binContent;
    Int_t roiLL = 7629;
    Int_t roiUL = 7829;
    
    Int_t sourceLL = 4000;
    Int_t sourceUL = 6700;
    
    Int_t counter = 0;
    Double_t totalRoiWith = 0;
    Double_t totalRoiNo = 0;
    
    // read the calibration data
    
    for(Int_t i = 0; i <= 210; i++ ){
        
        cuFitFile >> calibList[i];
  
    }

    
    for( Int_t sdd = 1; sdd <= 6; sdd++){

        histString = Form(withString + "Currentsdd" + "%d", sdd );
        //cout << histString << endl;
        TH1F *withH = (TH1F*)f->Get(histString);
        
        roiCountsWith[sdd-1] = 0;

        for( Int_t binNumber = 1; binNumber <= maxBinNumber; binNumber++ ){
            
            binCenter = withH->GetBinCenter(binNumber);

            if( binCenter > roiLL && binCenter < roiUL ){
                
                //cout << binCenter << endl;
                binContent = withH -> GetBinContent(binNumber);
               // counter += 1;
                roiCountsWith[sdd-1] += binContent;
                
                //cout << binContent << " " << counter << endl;
                
            }

        }

        withH->Delete();
        cout << "sdd: " << sdd << " ROI counts with Current: " <<  " " << roiCountsWith[sdd-1] << endl;
        totalRoiWith += roiCountsWith[sdd-1];
    
    }
    
    for( Int_t sdd = 1; sdd <= 6; sdd++){

        histString = Form(noString + "CurrentSmallsdd" + "%d", sdd );
        //cout << histString << endl;
        TH1F *noH = (TH1F*)f->Get(histString);
        
        roiCountsNo[sdd-1] = 0;

        for( Int_t binNumber = 1; binNumber <= maxBinNumber; binNumber++ ){
            
            binCenter = noH->GetBinCenter(binNumber);

            if( binCenter > roiLL && binCenter < roiUL ){
                
                binContent = noH -> GetBinContent(binNumber);
               // counter += 1;
                roiCountsNo[sdd-1] += binContent;
                
                //cout << binContent << " " << counter << endl;
                
            }

        }

        noH->Delete();
        cout << "sdd: " << sdd << " ROI counts no Current: " <<  " " << roiCountsNo[sdd-1] << endl;
        totalRoiNo += roiCountsNo[sdd-1];
    
    }
    
    TH1F *noHSum = (TH1F*)f->Get("noCurrentSmallSum");
    TH1F *withHSum = (TH1F*)f->Get("withCurrentSum");
    
    
    Double_t lowBin, highBin;
    
    lowBin = noHSum->FindBin(roiLL);
    highBin = noHSum->FindBin(roiUL);
    
    Double_t binCenterL = withHSum->GetBinCenter(lowBin);
    Double_t binCenterH = withHSum->GetBinCenter(highBin-1);
    
    Double_t binContentL = noHSum->GetBinContent(lowBin);
    Double_t binContentH = noHSum->GetBinContent(highBin-1);
    
    Double_t noSum = noHSum->Integral(lowBin,highBin-1);
    Double_t withSum = withHSum->Integral(lowBin,highBin-1);
    
    // start removing the tails of Cu and Ni from the ROI counts ---------------------------
    
    Double_t cuGainWith = calibList[(6*15+2) - 1];
    Double_t cuGainNo   = calibList[(13*15+2) - 1];
    
    Double_t niGainWith = calibList[(6*15+8) - 1];
    Double_t niGainNo   = calibList[(13*15+8) - 1];
    
    Double_t cuSigmaWith = calibList[(6*15+4) - 1];
    Double_t cuSigmaNo   = calibList[(13*15+4) - 1];
    
    Double_t niSigmaWith = calibList[(6*15+10) - 1];
    Double_t niSigmaNo   = calibList[(13*15+10) - 1];
    
    Double_t cuMeanWith = calibList[(6*15+3) - 1];
    Double_t cuMeanNo   = calibList[(13*15+3) - 1];
    
    Double_t niMeanWith = calibList[(6*15+9) - 1];
    Double_t niMeanNo   = calibList[(13*15+9) - 1];
    
    
    //cout << cuSigmaWith << " " << cuSigmaNo << " " << niSigmaWith << " " << niSigmaNo << " " <<  cuMeanWith << " " << cuMeanNo << " " << niMeanWith << " " << niMeanNo << endl;
    
    Double_t cuka1WithCont = calcRoiContribution(cuMeanWith, cuGainWith, cuSigmaWith, 7629, 7829, 25);
    Double_t cuka2WithCont = calcRoiContribution(cuMeanWith-19.95, cuGainWith * 0.51, cuSigmaWith, 7629, 7829, 25);
    
    Double_t nika1WithCont = calcRoiContribution(niMeanWith, niGainWith, niSigmaWith, 7629, 7829, 25);
    Double_t nika2WithCont = calcRoiContribution(niMeanWith-17.26, niGainWith * 0.51, niSigmaWith, 7629, 7829, 25);
    
    Double_t totWithCont = cuka1WithCont + cuka2WithCont + nika1WithCont + nika2WithCont;
    
    // ---------------------- no current contribution
    
    Double_t cuka1NoCont = calcRoiContribution(cuMeanNo, cuGainNo, cuSigmaNo, 7629, 7829, 25);
    Double_t cuka2NoCont = calcRoiContribution(cuMeanNo-19.95, cuGainNo * 0.51, cuSigmaNo, 7629, 7829, 25);
    
    Double_t nika1NoCont = calcRoiContribution(niMeanNo, niGainNo, niSigmaNo, 7629, 7829, 25);
    Double_t nika2NoCont = calcRoiContribution(niMeanNo-17.26, niGainNo * 0.51, niSigmaNo, 7629, 7829, 25);
    
    Double_t totNoCont = cuka1NoCont + cuka2NoCont + nika1NoCont + nika2NoCont;
    
    Double_t subtrCountsWithout = noSum - totNoCont;
    Double_t subtrCountsWith    = withSum - totWithCont;

    
    //cout << "total counts ROI with current: " << totalRoiWith <<endl;
    //cout << "total counts ROI no current: " << totalRoiNo <<endl;
    
    cout << "cu ka1 contribution in roi with current: " << cuka1WithCont << endl;
    cout << "cu ka2 contribution in roi with current: " << cuka2WithCont << endl;
    
    cout << "ni ka1 contribution in roi with current: " << nika1WithCont << endl;
    cout << "ni ka2 contribution in roi with current: " << nika2WithCont << endl;
    
    cout << "cu ka1 contribution in roi without current: " << cuka1NoCont << endl;
    cout << "cu ka2 contribution in roi without current: " << cuka2NoCont << endl;
    
    cout << "ni ka1 contribution in roi without current: " << nika1NoCont << endl;
    cout << "ni ka2 contribution in roi without current: " << nika2NoCont << endl;
    
    cout << "tot contribution with current: " << totWithCont<< endl;
    cout << "tot contribution without current: " << totNoCont<< endl << endl;
    
    cout << "total counts ROI no current: " << noSum <<endl;
    cout << "total counts ROI with current: " << withSum <<endl;
    
    cout << "total counts ROI no current minus gaussian contributions: " << subtrCountsWithout <<endl;
    cout << "total counts ROI with current minus gaussian contributions: " << subtrCountsWith <<endl;
    
    
    cuFitFile.close();
    f->Close();
    
    

}

Double_t calcEnergy(Short_t channel, Double_t slope, Double_t offset, Double_t Ord2Coeff){
    
    
    Double_t energy = -( slope - 2*Ord2Coeff*MnKa1-TMath::Sqrt(4*Ord2Coeff*channel-4*Ord2Coeff*offset+slope*slope-4*Ord2Coeff*slope*MnKa1))/( 2*Ord2Coeff );
    
    return energy;
    
}

Double_t calcChannelWidth(Double_t energy, Double_t slope, Double_t Ord2Coeff){
    
    Double_t channelWidth = 1 / (slope + 2 * Ord2Coeff * (energy - MnKa1) );
    
    //cout << slope << " " << Ord2Coeff << " " << energy << " " << channelWidth << endfcalcroil;
    
    return channelWidth;
}

void AddCalibrationBranches(Int_t minFileNumber, Int_t maxFileNumber){
    
    // takes the information of many smaller rootfiles and writes the ttree to other rootfiles
    // adds the information of the fit results and the energy of every event 
    // Use $ROOTSYS/bin/hadd

//hadd result.root file1.root file2.root .. fileN.root   to add rootfiles with same structure

    Double_t calibList[150];
    Short_t adcChannel[16];
    Double_t energyEv[16];
    Double_t binLowE, binHighE;
    TString fileName, rootFileName, writeFileName;
    TRandom3 randG;
    Double_t randN;
    Int_t current;
    Double_t energyTmp, channelWidthTmp;
    
    Double_t Ord2CoeffWith[6] = {5.34213e-07,1.3487e-07,2.92066e-07,3.57834e-07,5.28058e-07,3.6126e-07};
    Double_t Ord2CoeffWithout[6] = {3.81858e-07,0.,4.88842e-07,3.11132e-07,4.51601e-07,3.49158e-07};
    
    ifstream calibFile;
    ifstream fileList;
    ofstream logFile;
    fileList.open(WORK_PATH + "/filelist/1-618Files-LNGS.txt");
    calibFile.open(WORK_PATH + "/calibrationlist/1-618Files-LNGS.txt");
    logFile.open("/home/andreas/vip2/reports/Analysis/logfile.txt");
    
    for( Int_t fileCount = 1; fileCount <= maxFileNumber; fileCount++){
        
        fileList >> fileName;
        
        current = 0;
        if (fileName.Contains("with")){ current = 1;}
        
        for( Int_t sdd = 1; sdd <= 6; sdd++ ){

            for( Int_t i = 0; i <= 24; i++ ){

                calibFile >> calibList[(sdd-1) * 25 + i];

            }
        }
        
        if( fileCount < minFileNumber ){ continue; } 
        
        rootFileName =  ROOT_PATH_LNGS + "/1-618files-original/" + fileName + ".root"; 
        //rootFileName =  ROOT_PATH_LNGS + "/1-618files-original/544withCurrent.root";
        //writeFileName = ROOT_PATH_LNGS + "/1-618files-calibrated/544withCurrent.root";
        writeFileName =  ROOT_PATH_LNGS + "/1-618files-calibrated/" + fileName + ".root";

        //rootFileName = ROOT_PATH_LNGS + "/test/1noCurrent.root";
        //rootFileName = ROOT_PATH_LNGS + "/1noCurrent.root";
        cout << rootFileName << endl;
        TFile *froot = TFile::Open(rootFileName, "READ");
        TTree *tree = (TTree*)froot->Get("tr");

        tree->SetBranchAddress("adc",adcChannel);


//        TObjArray *branchArray = tree->GetListOfBranches();
//        Int_t nbranches = branchArray->GetEntries();
//        if( nbranches > 17 ){ 
//
//            TBranch *calOld = (TBranch*)tree->GetBranch("cal");
//            TBranch *energyOld = (TBranch*)tree->GetBranch("energy");
//            branchArray->Remove(calOld);
//            branchArray->Remove(energyOld);
//            
//            //tree->Write("",TObject::kOverwrite);
//            froot->Write("", TObject::kOverwrite);
//            tree->Print();

//        }


        TBranch *cal = tree->Branch("cal",calibList,"Calibration Parameter[150]/D");
        TBranch *energy = tree->Branch("energy",energyEv,"Energy in eV[16]/D");
        Double_t nentries = tree->GetEntries();

        

        for( Int_t i = 0; i < nentries; i++ ){
            
            tree->GetEntry(i);
            
            binLowE = ((adcChannel[0] - calibList[4]) / calibList[2]) - 1/(2 * calibList[2]);
            binHighE = ((adcChannel[0] - calibList[4]) / calibList[2]) + 1/(2 * calibList[2]);
            randN = randG.Rndm();
            energyEv[0] = binLowE + randN * ( binHighE - binLowE );
            
            if( energyEv[0] > MnKa1 ){
                
                if(current == 0){
                    
                    energyTmp = calcEnergy(adcChannel[0], calibList[2], calibList[4], Ord2CoeffWithout[0]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[2], Ord2CoeffWithout[0]);

                }else{
                
                    energyTmp = calcEnergy(adcChannel[0], calibList[2], calibList[4], Ord2CoeffWith[0]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[2], Ord2CoeffWith[0]);

                }

                binLowE = energyTmp - channelWidthTmp/2;
                binHighE = energyTmp + channelWidthTmp/2;
                
                energyEv[0] = binLowE + randN * ( binHighE - binLowE );

            }
            
            
            
            binLowE = ((adcChannel[2] - calibList[29]) / calibList[27]) - 1/(2 * calibList[27]);
            binHighE = ((adcChannel[2] - calibList[29]) / calibList[27]) + 1/(2 * calibList[27]);
            randN = randG.Rndm();
            energyEv[2] = binLowE + randN * ( binHighE - binLowE );
            
            if( (energyEv[2] > MnKa1) && (current == 1) ){ // change here ?? ... for sdd 2 no current the 2nd order parameter == 0 -> linear fit

                if(current == 1){
                    
                    energyTmp = calcEnergy(adcChannel[2], calibList[27], calibList[29], Ord2CoeffWithout[1]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[27], Ord2CoeffWithout[1]);

                }

                binLowE = energyTmp - channelWidthTmp/2;
                binHighE = energyTmp + channelWidthTmp/2;
                
                energyEv[2] = binLowE + randN * ( binHighE - binLowE );

            }
            
            binLowE = ((adcChannel[3] - calibList[54]) / calibList[52]) - 1/(2 * calibList[52]);
            binHighE = ((adcChannel[3] - calibList[54]) / calibList[52]) + 1/(2 * calibList[52]);
            randN = randG.Rndm();
            energyEv[3] = binLowE + randN * ( binHighE - binLowE );
            
            if( energyEv[3] > MnKa1 ){
                
                if(current == 0){
                    
                    energyTmp = calcEnergy(adcChannel[3], calibList[52], calibList[54], Ord2CoeffWithout[2]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[52], Ord2CoeffWithout[2]);

                }else{
                
                    energyTmp = calcEnergy(adcChannel[3], calibList[52], calibList[54], Ord2CoeffWith[2]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[52], Ord2CoeffWith[2]);

                }

                binLowE = energyTmp - channelWidthTmp/2;
                binHighE = energyTmp + channelWidthTmp/2;
                
                energyEv[3] = binLowE + randN * ( binHighE - binLowE );

            }
            
            binLowE = ((adcChannel[5] - calibList[79]) / calibList[77]) - 1/(2 * calibList[77]);
            binHighE = ((adcChannel[5] - calibList[79]) / calibList[77]) + 1/(2 * calibList[77]);
            randN = randG.Rndm();
            energyEv[5] = binLowE + randN * ( binHighE - binLowE );
            
            if( energyEv[5] > MnKa1 ){
                
                if(current == 0){
                    
                    energyTmp = calcEnergy(adcChannel[5], calibList[77], calibList[79], Ord2CoeffWithout[3]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[77], Ord2CoeffWithout[3]);

                }else{
                
                    energyTmp = calcEnergy(adcChannel[5], calibList[77], calibList[79], Ord2CoeffWith[3]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[77], Ord2CoeffWith[3]);

                }

                binLowE = energyTmp - channelWidthTmp/2;
                binHighE = energyTmp + channelWidthTmp/2;
                
                energyEv[5] = binLowE + randN * ( binHighE - binLowE );

            }
            
            binLowE = ((adcChannel[6] - calibList[104]) / calibList[102]) - 1/(2 * calibList[102]);
            binHighE = ((adcChannel[6] - calibList[104]) / calibList[102]) + 1/(2 * calibList[102]);
            randN = randG.Rndm();
            energyEv[6] = binLowE + randN * ( binHighE - binLowE );
            
            if( energyEv[6] > MnKa1 ){
                
                if(current == 0){
                    
                    energyTmp = calcEnergy(adcChannel[6], calibList[102], calibList[104], Ord2CoeffWithout[4]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[102], Ord2CoeffWithout[4]);

                }else{
                
                    energyTmp = calcEnergy(adcChannel[6], calibList[102], calibList[104], Ord2CoeffWith[4]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[102], Ord2CoeffWith[4]);

                }

                binLowE = energyTmp - channelWidthTmp/2;
                binHighE = energyTmp + channelWidthTmp/2;
                
                energyEv[6] = binLowE + randN * ( binHighE - binLowE );

            }
            
            binLowE = ((adcChannel[7] - calibList[129]) / calibList[127]) - 1/(2 * calibList[127]);
            binHighE = ((adcChannel[7] - calibList[129]) / calibList[127]) + 1/(2 * calibList[127]);
            randN = randG.Rndm();
            energyEv[7] = binLowE + randN * ( binHighE - binLowE );
            
            if( energyEv[7] > MnKa1 ){
                
                if(current == 0){
                    
                    energyTmp = calcEnergy(adcChannel[7], calibList[127], calibList[129], Ord2CoeffWithout[5]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[127], Ord2CoeffWithout[5]);

                }else{
                
                    energyTmp = calcEnergy(adcChannel[7], calibList[127], calibList[129], Ord2CoeffWith[5]);
                    channelWidthTmp = calcChannelWidth(energyTmp, calibList[127], Ord2CoeffWith[5]);

                }

                binLowE = energyTmp - channelWidthTmp/2;
                binHighE = energyTmp + channelWidthTmp/2;
                
                energyEv[7] = binLowE + randN * ( binHighE - binLowE );

            }
            
            

            energy->Fill();
            cal->Fill();

        }
        
        if( calibList[6] > 200 || calibList[31] > 200 || calibList[56] > 200 || calibList[71] > 200 || calibList[96] > 200 || calibList[121] > 200 ){
                
            logFile << "fwhm too high in file: " << rootFileName << endl;
                
        }
            
        if( calibList[2] < 0.08 || calibList[27] < 0.08 || calibList[52] < 0.08 || calibList[77] < 0.08 || calibList[102] < 0.08 || calibList[127] < 0.08 ){
                
            logFile << "slope too low in file: " << rootFileName << endl;
                
        }
            
        if( calibList[12] < 600 || calibList[37] < 0.08 || calibList[62] < 600 || calibList[87] < 600 || calibList[112] < 600 || calibList[137] < 600 ){
                
            logFile << "mnka1 position too low in file: " << rootFileName << endl;
                
        }

        TFile *fwrite = TFile::Open(writeFileName, "NEW");
        TTree *tr = tree->CloneTree();
        
        //tree->Print();
        //gDirectory->pwd();
        //fwrite->ls();
        //tree->Write("",TObject::kOverwrite);
        //fwrite->Write("", TObject::kOverwrite);
        
        tr->Write();
        
        //tree->Print();
        froot->Close();
        fwrite->Close();
        delete fwrite;
        delete froot;
    }
    
    calibFile.close();
    fileList.close();
    logFile.close();
    
    
    
    return;
}

void PlotBackground(){
    
    ifstream bgFile, bgFileWith;
    bgFile.open("/home/andreas/vip2/reports/Analysis/backgroundNoCurrent.txt");
    bgFileWith.open("/home/andreas/vip2/reports/Analysis/backgroundWithCurrent.txt");
    
    Int_t sddN;
    Int_t splitsWith = 81;
    Int_t splits = 116;
    
    Double_t utList[splits];
    Double_t fileNList[splits];
    Double_t sourceR1[splits], backgR1[splits], highER1[splits], rejR1[splits], cuR1[splits];
    Double_t sourceR2[splits], backgR2[splits], highER2[splits], rejR2[splits], cuR2[splits];
    Double_t sourceR3[splits], backgR3[splits], highER3[splits], rejR3[splits], cuR3[splits];
    Double_t sourceR4[splits], backgR4[splits], highER4[splits], rejR4[splits], cuR4[splits];
    Double_t sourceR5[splits], backgR5[splits], highER5[splits], rejR5[splits], cuR5[splits];
    Double_t sourceR6[splits], backgR6[splits], highER6[splits], rejR6[splits], cuR6[splits];
    
    Double_t utListWith[splitsWith];
    Double_t fileNListWith[splitsWith];
    Double_t sourceR1With[splitsWith], backgR1With[splitsWith], highER1With[splitsWith], rejR1With[splitsWith], cuR1With[splitsWith];
    Double_t sourceR2With[splitsWith], backgR2With[splitsWith], highER2With[splitsWith], rejR2With[splitsWith], cuR2With[splitsWith];
    Double_t sourceR3With[splitsWith], backgR3With[splitsWith], highER3With[splitsWith], rejR3With[splitsWith], cuR3With[splitsWith];
    Double_t sourceR4With[splitsWith], backgR4With[splitsWith], highER4With[splitsWith], rejR4With[splitsWith], cuR4With[splitsWith];
    Double_t sourceR5With[splitsWith], backgR5With[splitsWith], highER5With[splitsWith], rejR5With[splitsWith], cuR5With[splitsWith];
    Double_t sourceR6With[splitsWith], backgR6With[splitsWith], highER6With[splitsWith], rejR6With[splitsWith], cuR6With[splitsWith];
    
    for( Int_t fileN = 1; fileN <= splits; fileN++ ){
        
        bgFile >> sddN >> utList[fileN-1] >> sourceR1[fileN-1] >> backgR1[fileN-1] >> highER1[fileN-1] >> cuR1[fileN-1] >> rejR1[fileN-1];
        bgFile >> sddN >> utList[fileN-1] >> sourceR2[fileN-1] >> backgR2[fileN-1] >> highER2[fileN-1] >> cuR2[fileN-1] >> rejR2[fileN-1];
        bgFile >> sddN >> utList[fileN-1] >> sourceR3[fileN-1] >> backgR3[fileN-1] >> highER3[fileN-1] >> cuR3[fileN-1] >> rejR3[fileN-1];
        bgFile >> sddN >> utList[fileN-1] >> sourceR4[fileN-1] >> backgR4[fileN-1] >> highER4[fileN-1] >> cuR4[fileN-1] >> rejR4[fileN-1];
        bgFile >> sddN >> utList[fileN-1] >> sourceR5[fileN-1] >> backgR5[fileN-1] >> highER5[fileN-1] >> cuR5[fileN-1] >> rejR5[fileN-1];
        bgFile >> sddN >> utList[fileN-1] >> sourceR6[fileN-1] >> backgR6[fileN-1] >> highER6[fileN-1] >> cuR6[fileN-1] >> rejR6[fileN-1];
        
        //fileNList[fileN-1] = fileN;
        
//        cout << backgR1[fileN-1] << " " << utList[fileN-1] << endl;
        
        
    }
    
    for( Int_t fileN = 1; fileN <= splitsWith; fileN++ ){
        
        bgFileWith >> sddN >> utListWith[fileN-1] >> sourceR1With[fileN-1] >> backgR1With[fileN-1] >> highER1With[fileN-1] >> cuR1With[fileN-1] >> rejR1With[fileN-1];
        bgFileWith >> sddN >> utListWith[fileN-1] >> sourceR2With[fileN-1] >> backgR2With[fileN-1] >> highER2With[fileN-1] >> cuR2With[fileN-1] >> rejR2With[fileN-1];
        bgFileWith >> sddN >> utListWith[fileN-1] >> sourceR3With[fileN-1] >> backgR3With[fileN-1] >> highER3With[fileN-1] >> cuR3With[fileN-1] >> rejR3With[fileN-1];
        bgFileWith >> sddN >> utListWith[fileN-1] >> sourceR4With[fileN-1] >> backgR4With[fileN-1] >> highER4With[fileN-1] >> cuR4With[fileN-1] >> rejR4With[fileN-1];
        bgFileWith >> sddN >> utListWith[fileN-1] >> sourceR5With[fileN-1] >> backgR5With[fileN-1] >> highER5With[fileN-1] >> cuR5With[fileN-1] >> rejR5With[fileN-1];
        bgFileWith >> sddN >> utListWith[fileN-1] >> sourceR6With[fileN-1] >> backgR6With[fileN-1] >> highER6With[fileN-1] >> cuR6With[fileN-1] >> rejR6With[fileN-1];
        
        //fileNList[fileN-1] = fileN;
        
        cout << backgR1With[fileN-1] << " " << utListWith[fileN-1] << endl;
        
        
    }
    
    bgFile.close();
    
    
    TCanvas *c1 = new TCanvas("SDD1", "SDD1", 800, 600);
    c1->Divide(2,2);
    
    c1->cd(1);
    TGraph *grb1 = new TGraph(splits,utList,backgR1);
    grb1->SetMarkerStyle(20);
    grb1->GetYaxis()->SetRangeUser(0.0075,0.02);
    grb1->Draw();
    
    TGraph *grb1With = new TGraph(splitsWith,utListWith,backgR1With);
    grb1With->SetMarkerStyle(20);
    grb1With->SetMarkerColor(4);
    grb1With->Draw("same");
    
//    c1->cd(2);
//    TGraph *grr1 = new TGraph(splits,utList,rejR1);
//    grr1->SetMarkerStyle(20);
//    grr1->Draw();
    
    c1->cd(2);
    TGraph *grcu1 = new TGraph(splits,utList,cuR1);
    grcu1->SetMarkerStyle(20);
    grcu1->Draw();
    
    c1->cd(3);
    TGraph *grs1 = new TGraph(splits,utList,sourceR1);
    grs1->SetMarkerStyle(20);
    grs1->Draw();
    
    c1->cd(4);
    TGraph *grhe1 = new TGraph(splits,utList,highER1);
    grhe1->SetMarkerStyle(20);
    grhe1->Draw();
    
    
    // SDD 2 ------------------------------------------
    
//    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
//    c2->Divide(2,2);
//    
//    c2->cd(1);
//    TGraph *grb2 = new TGraph(splits,utList,backgR2);
//    grb2->SetMarkerStyle(20);
//    grb2->Draw();
//    
////    c2->cd(2);
////    TGraph *grr2 = new TGraph(splits,utList,rejR2);
////    grr2->SetMarkerStyle(20);
////    grr2->Draw();
//    
//    c2->cd(2);
//    TGraph *grcu2 = new TGraph(splits,utList,cuR2);
//    grcu2->SetMarkerStyle(20);
//    grcu2->Draw();
//    
//    c2->cd(3);
//    TGraph *grs2 = new TGraph(splits,utList,sourceR2);
//    grs2->SetMarkerStyle(20);
//    grs2->Draw();
//    
//    c2->cd(4);
//    TGraph *grhe2 = new TGraph(splits,utList,highER2);
//    grhe2->SetMarkerStyle(20);
//    grhe2->Draw();
//    
//    // SDD 3 ------------------------------------------
//    
//    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
//    c3->Divide(2,2);
//    
//    c3->cd(1);
//    TGraph *grb3 = new TGraph(splits,utList,backgR3);
//    grb3->SetMarkerStyle(20);
//    grb3->Draw();
//    
////    c3->cd(2);
////    TGraph *grr3 = new TGraph(splits,utList,rejR3);
////    grr3->SetMarkerStyle(20);
////    grr3->Draw();
//    
//    c3->cd(2);
//    TGraph *grcu3 = new TGraph(splits,utList,cuR3);
//    grcu3->SetMarkerStyle(20);
//    grcu3->Draw();
//    
//    c3->cd(3);
//    TGraph *grs3 = new TGraph(splits,utList,sourceR3);
//    grs3->SetMarkerStyle(20);
//    grs3->Draw();
//    
//    c3->cd(4);
//    TGraph *grhe3 = new TGraph(splits,utList,highER3);
//    grhe3->SetMarkerStyle(20);
//    grhe3->Draw();
//    
//    // SDD 4 ------------------------------------------
//    
//    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
//    c4->Divide(2,2);
//    
//    c4->cd(1);
//    TGraph *grb4 = new TGraph(splits,utList,backgR4);
//    grb4->SetMarkerStyle(20);
//    grb4->Draw();
//    
////    c4->cd(2);
////    TGraph *grr4 = new TGraph(splits,utList,rejR4);
////    grr4->SetMarkerStyle(20);
////    grr4->Draw();
//    
//    c4->cd(2);
//    TGraph *grcu4 = new TGraph(splits,utList,cuR4);
//    grcu4->SetMarkerStyle(20);
//    grcu4->Draw();
//    
//    c4->cd(3);
//    TGraph *grs4 = new TGraph(splits,utList,sourceR4);
//    grs4->SetMarkerStyle(20);
//    grs4->Draw();
//    
//    c4->cd(4);
//    TGraph *grhe4 = new TGraph(splits,utList,highER4);
//    grhe4->SetMarkerStyle(20);
//    grhe4->Draw();
//    
//    // SDD 5 ------------------------------------------
//    
//    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
//    c5->Divide(2,2);
//    
//    c5->cd(1);
//    TGraph *grb5 = new TGraph(splits,utList,backgR5);
//    grb5->SetMarkerStyle(20);
//    grb5->Draw();
//    
////    c5->cd(2);
////    TGraph *grr5 = new TGraph(splits,utList,rejR5);
////    grr5->SetMarkerStyle(20);
////    grr5->Draw();
//    
//    c5->cd(2);
//    TGraph *grcu5 = new TGraph(splits,utList,cuR5);
//    grcu5->SetMarkerStyle(20);
//    grcu5->Draw();
//    
//    c5->cd(3);
//    TGraph *grs5 = new TGraph(splits,utList,sourceR5);
//    grs5->SetMarkerStyle(20);
//    grs5->Draw();
//    
//    c5->cd(4);
//    TGraph *grhe5 = new TGraph(splits,utList,highER5);
//    grhe5->SetMarkerStyle(20);
//    grhe5->Draw();
//    
//    // SDD 6 ------------------------------------------
//    
//    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
//    c6->Divide(2,2);
//    
//    c6->cd(1);
//    TGraph *grb6 = new TGraph(splits,utList,backgR6);
//    grb6->SetMarkerStyle(20);
//    grb6->Draw();
//    
////    c6->cd(2);
////    TGraph *grr6 = new TGraph(splits,utList,rejR6);
////    grr6->SetMarkerStyle(20);
////    grr6->Draw();
//    
//    c6->cd(2);
//    TGraph *grcu6 = new TGraph(splits,utList,cuR6);
//    grcu6->SetMarkerStyle(20);
//    grcu6->Draw();
//    
//    c6->cd(3);
//    TGraph *grs6 = new TGraph(splits,utList,sourceR6);
//    grs6->SetMarkerStyle(20);
//    grs6->Draw();
//    
//    c6->cd(4);
//    TGraph *grhe6 = new TGraph(splits,utList,highER6);
//    grhe6->SetMarkerStyle(20);
//    grhe6->Draw();
    
    
    
}


void WriteParts2Rootfile(TString rootfile, TString place, Int_t divider, Int_t part){

   // writes parts of a given rootfile (like the 1st tenth of it) to another rootfile - currently only the tree and the ADC Branch filled with data only from sdd 1
 
    EventStruct  evt_r;
    EventStruct  evt_w;

    Int_t* ind_array;
    ind_array = GetTreeDivisionIndices(rootfile,place,divider,part);

    //TH1F *hev = new TH1F("hev", "Energy spectrum", 4096, -0.5, 4095.5);

    TString rootfilename;

    Int_t size = place.Sizeof();
    if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
    if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;

// rootfile to read from
    TFile *fr = new TFile(rootfilename, "READ");

    TTree *tr_r = (TTree*)fr->Get("tr"); // tr_r is now a pointer to the tree I want to read

    TBranch *adcEvent   = tr_r->GetBranch("adc"); // adc Event is a pointer to the Branch I want to read
    adcEvent->SetAddress(evt_r.padc); // evt_r padc is from the tree I want to read to

// rootfile to write to
    TFile *fw = new TFile( Form("/home/andreas/vip2/data/root/LNGS/01day-parts/2noCurrent_%d.root",part),"RECREATE");
    fw->cd();

    TTree *tr_w = new TTree("tr", "SDD data");
    tr_w->Branch("adc",   evt_w.padc,   "ADC Channel Data[16]/S"); // tr_w is the tree I want to write to; with a Branch 
 
    cout << ind_array[0] << " " << ind_array[1] << endl;

   
        
        

    for( Int_t i=ind_array[0]; i<=ind_array[1]; i++ ){
    //for( Int_t i=1; i<10; i++ ){
        
        adcEvent->GetEvent(i);
        
        //cout << "i " << i << endl;
        
        for( Int_t adcChannel=0; adcChannel<16; adcChannel++ ){    
                
            evt_w.padc[adcChannel] = evt_r.padc[adcChannel];
            //cout << "adc Entry " << evt_w.padc[adcChannel] << endl;

        }
        
        tr_w->Fill();
    }
     
    

    fw->Write();
    //TFile hist_file(Form("/home/andreas/vip2/reports/1608_VIPReportLNF/sdd1_part%dLNGS.root",part), "CREATE");
    //hev->Write(); 
    //hist_file.Close();
    fr->Close();
    delete fr;
    fw->Close();
    delete fw;

    return;
   
}

void WriteParts2RootfileLoop( TString rootfile, TString place, Int_t divider ){
    
    for(Int_t part = 1; part<=80; part++){
    
        WriteParts2Rootfile( rootfile, place, divider, part );
    
    }
    
    return;
}

void CheckBackground( TString rootfile,  TString place, Int_t divider, Int_t part, Int_t complete, Int_t constFracTime, Int_t sec_rootfile){
    
    // if complete == 1 .... the whole rootfile is looked at; 
    // if constFracTime == 1 the file is divided into bits with about 6 hours or so; the part needs to be taken care of manually that it does not happen that part > divider
    // if sec_rootfile == -1 ... the duration of the complete rootfile is calculated with getruntime(); else sec_rootfile should be the duration of the complete rootfile in seconds
    // rootfile_loop_count == 1 ... rootfile is openend; number gives the count of loops over the parts of this rootfile
    // slope and offset are calculated from the complete root file, as the 5 hour pieces might not be enough to make a calibration
    //complete must be 0 or 1
    
    ofstream logfile;   
//    TString logfileName = ANALYSIS_PATH + "/logfile.txt";
    TString logfileName = ANALYSIS_PATH + "/backgroundDatafile.txt";
    
    //cout << logfileName << endl;
    
    logfile.open(logfileName,std::ofstream::app);
    
    Int_t sec_duration, sec_duration_total, hour_duration;
    Double_t slope[6], offset[6];
    Int_t sdd;
    Int_t adc_channel;
    Int_t* ind_array;
    Int_t lbound_ind = 0;
    Int_t ubound_ind = 0;
    Double_t bin_ev[6];
    Int_t coinc_counter = 0;
    Double_t coinc_frac;
    Double_t energyList[16];
    Int_t sddEventCount = 0;
    
    Double_t avg_rate_bg[6] = {0.0082, 0.011, 0.008, 0.0062, 0.0083, 0.0074};
    
    Int_t hours_target, hours_rootfile;
    
    hours_target = 24;
    
    Int_t lbound_source = 4000;
    Int_t ubound_source = 6700;
    Int_t counter_source[6] = {0};
    Double_t rate_source[6] = {0};
    
    Int_t lbound_bg = 7000;
    Int_t ubound_bg = 30000;
    Int_t counter_bg[6] = {0};
    Double_t rate_bg[6] = {0};
    
    Int_t lbound_he = 30000;
    Int_t ubound_he = 70000;
    Int_t counter_he[6] = {0};
    Double_t rate_he[6] = {0};
    
    Int_t lbound_cu = 7000;
    Int_t ubound_cu = 10000;
    Int_t counter_cu[6] = {0};
    Double_t rate_cu[6] = {0};
    
    EventStruct  evt_r;

    Int_t size = place.Sizeof();
    TString rootfilename;

    if(size==4)rootfilename = ROOT_PATH_SMI + "/" + rootfile;
    if(size==5)rootfilename = ROOT_PATH_LNGS + "/" + rootfile;
    
    if( sec_rootfile == -1 ){sec_duration_total = GetRunTime(rootfile, place, 1, 1, 1);}
    else( sec_duration_total = sec_rootfile );
    
    
    hours_rootfile = (int) sec_duration_total / 3600;
    
    if( constFracTime == 1 ) { divider = (int) hours_rootfile / hours_target; }
    
    sec_duration = GetRunTime(rootfile, place, divider, part, complete);

    
//    for( sdd = 1; sdd < 7; sdd++  ){
//    
//        slope[sdd-1] = GetSlope(rootfile, place, sdd);
//        offset[sdd-1] = GetOffset(rootfile, place, sdd);
//    
//    }

    
    //TFile *f;
    TFile *f = TFile::Open(rootfilename, "READ");
    //if( rootfile_loop_count == 1 ){
    //TFile  *f = new TFile(rootfilename,"read");
        //cout << " OPENING ROOT FILE" << endl;
    //}
    
    TTree *tree = (TTree*)f->Get("tr");
    
    TBranch *adcEvent   = tree->GetBranch("adc");
    TBranch *timeEvent  = tree->GetBranch("time"); 
    TBranch *utEvent  = tree->GetBranch("ut"); 
    TBranch *trgidEvent = tree->GetBranch("trgid");
    TBranch *energyEvent = tree->GetBranch("energy");
    
    adcEvent->SetAddress(evt_r.padc);
    timeEvent->SetAddress(evt_r.time);
    trgidEvent->SetAddress(&evt_r.trgid);
    utEvent->SetAddress(&evt_r.ut);
    energyEvent->SetAddress(energyList);
    
    Int_t nEvent = tree->GetEntries();
    
    if( part > divider ){ cout << " PART > DIVIDER; THIS SHOULD NOT BE !!!"  << endl; goto end;}
    
    if( complete == 1 || divider == 1){ lbound_ind = 1; ubound_ind = nEvent; }
    else{ ind_array =  GetTreeDivisionIndices(rootfile, place, divider, part); 
  
        lbound_ind = ind_array[0];
        ubound_ind = ind_array[1];
  
    }
    
    cout << " lower event index: " << lbound_ind << " upper event index:  " << ubound_ind << endl;
    
    for(Int_t i = lbound_ind; i <= ubound_ind; i++){
        
        adcEvent->GetEntry(i);
        trgidEvent->GetEntry(i);
        utEvent->GetEntry(i);
        energyEvent->GetEvent(i);
        
        if( evt_r.trgid > 2 ){ coinc_counter += 1; }
        if( evt_r.trgid == 1 ){ sddEventCount += 1; }
        
        for( sdd = 1; sdd < 7; sdd++  ){
            
            adc_channel = SDDToPadc[sdd];
            
            //bin_ev[sdd-1] = evt_r.padc[adc_channel] * slope[sdd-1] + offset[sdd-1];
        
//            if( bin_ev[sdd-1] > lbound_source &&  bin_ev[sdd-1] < ubound_source){ counter_source[sdd-1] += 1; }
//            if( bin_ev[sdd-1] > lbound_bg &&  bin_ev[sdd-1] < ubound_bg)        { counter_bg[sdd-1] += 1;     }
//            if( bin_ev[sdd-1] > lbound_he &&  bin_ev[sdd-1] < ubound_he)        { counter_he[sdd-1] += 1;     }
            
            if( energyList[adc_channel] > lbound_source &&  energyList[adc_channel] < ubound_source){ counter_source[sdd-1] += 1; }
            if( energyList[adc_channel] > lbound_bg &&  energyList[adc_channel] < ubound_bg)        { counter_bg[sdd-1] += 1;     }
            if( energyList[adc_channel] > lbound_he &&  energyList[adc_channel] < ubound_he)        { counter_he[sdd-1] += 1;     }
            if( energyList[adc_channel] > lbound_cu &&  energyList[adc_channel] < ubound_cu)        { counter_cu[sdd-1] += 1;     }
        
        }
        
    }
    
    cout << "coinc counter: " << coinc_counter << " sdd Event count: " << sddEventCount << endl;
    coinc_frac = (float)coinc_counter / (float)(sddEventCount);
    
    
//    if( coinc_frac > 0.025 ){ 
//        
//        cout << " HIGH rejection rate: " << coinc_frac << endl;
//        
//        timeEvent->GetEntry(lbound_ind);              
//        cout << " START: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//            
//        timeEvent->GetEntry(ubound_ind); 
//        cout << " END: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//        
//        logfile << " HIGH rejection rate: " << coinc_frac << endl;
//        
//        timeEvent->GetEntry(lbound_ind);              
//        logfile << " START: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//            
//        timeEvent->GetEntry(ubound_ind); 
//        logfile << " END: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//    
//    }
//    
//    for( sdd = 1; sdd < 7; sdd++ ){
//        
//        rate_source[sdd-1] = (float)counter_source[sdd-1] / (float)sec_duration;
//        rate_bg[sdd-1] = (float)counter_bg[sdd-1] / (float)sec_duration;
//        rate_he[sdd-1] = (float)counter_he[sdd-1] / (float)sec_duration;
//
//        
//        if( rate_bg[sdd-1] > (avg_rate_bg[sdd-1] * 1.5) ){
//            
//            cout << "SDD: " << sdd << " source rate: " << rate_source[sdd-1] << " bg rate: " << rate_bg[sdd-1] << " high energy rate: " << rate_he[sdd-1] << endl;
//            
//            timeEvent->GetEntry(lbound_ind);              
//            cout << " START: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//            
//            timeEvent->GetEntry(ubound_ind); 
//            cout << " END: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//            
//            logfile << "SDD: " << sdd << " source rate: " << rate_source[sdd-1] << " bg rate: " << rate_bg[sdd-1] << " high energy rate: " << rate_he[sdd-1] << endl;
//            
//            timeEvent->GetEntry(lbound_ind);              
//            logfile << " START: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//            
//            timeEvent->GetEntry(ubound_ind); 
//            logfile << " END: year: " << evt_r.time[0] << " month: " <<  evt_r.time[1] << " day: " << evt_r.time[2] << " hour: " << evt_r.time[3] << " minute: " << evt_r.time[4] << endl;
//        }
//    
//    }
    
    cout << " duration of this part in seconds: " << sec_duration << endl;
    for( Int_t sdd = 1; sdd <= 6; sdd++ ){
        
        rate_source[sdd-1] = (float)counter_source[sdd-1] / (float)sec_duration;
        rate_bg[sdd-1] = (float)counter_bg[sdd-1] / (float)sec_duration;
        rate_he[sdd-1] = (float)counter_he[sdd-1] / (float)sec_duration;
        rate_cu[sdd-1] = (float)counter_cu[sdd-1] / (float)sec_duration;
        
        
        
        logfile << sdd << " " << evt_r.ut << " " << rate_source[sdd-1] << " " << rate_bg[sdd-1] << " " << rate_he[sdd-1] << " " << rate_cu[sdd-1] << " " << coinc_frac << endl;
        
        
    }
    
    cout << endl;
    end:
    logfile.close();
    f->Close();
    delete f;
    
    return;
    
    
}

void CheckBackgroundFile(TString rootfile,  TString place, Int_t divider, Int_t constFracTime){
    
    // constFracTime = 1 ... divider is newly calculated from the length of the rootfile and the target length of 1 part of "hours-target"; hours also needs to be reset in CheckBackground function!!
    // else: divider parameter is taken
    
    Int_t sec_duration_total;
    Int_t hours_target = 24;
    Int_t hours_rootfile;
    
    
    
    if( constFracTime == 1 ){
        
        sec_duration_total = GetRunTime(rootfile, place, 1, 1, 1);
    
        hours_rootfile = (int) sec_duration_total / 3600;
        
        divider = (int) hours_rootfile / hours_target;
        
    }
    
    //cout << "divider: " << divider << endl;
    
    for( Int_t part = 1; part <= divider; part ++ ){
        
        
        CheckBackground( rootfile,  place, divider, part, 0, constFracTime, sec_duration_total);
        
    }
    
    return;
    
}

void CheckBackgroundFilelist(TString filelistFile){
    
    
    FILE  *flist;
    //Int_t size = place.Sizeof();

    if( (flist = fopen(LIST_PATH + "/" +filelistFile, "r")) == NULL ){
        cout << "Cannot open file: " << filelistFile << endl;
    }

    cout << filelistFile << endl; // filelistFile is the name of the file which contains a list of strings, which correspond with a ".list" ending to the filenames of the files which include a list of 
    //binary files (already existing), and with a ".root ending to the root files

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
        
        
        //cout << rootfilename << " " << filelist << endl;
        
        CheckBackgroundFile(rootfilename,  "lngs", 1, 1);

        //cout << "Return value : " << read_status << endl << endl;
        cout << endl << endl << endl;
    }
    
    fclose(flist);
    return;
    
    
}

//void CheckBackgroundLoop(TString rootfile, )



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



  TH1F *histo_unscaled;
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

    // prints average rate in all int divider parts of the rootfile 
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

void MakeROIPlots(TString withName, TString noName){

  TFile *f = new TFile("/home/andreas/vip2/data/root/LNGS/1-618files-final/energyHistograms.root");
  
  TH1F *withH =  (TH1F*)f->Get(withName);
  TH1F *noH   =  (TH1F*)f->Get(noName);
  
  noH->Rebin(25);
  noH->GetXaxis()->SetRangeUser(7300,8400);
  noH->GetYaxis()->SetTitle("Counts/25 eV");
  noH->GetXaxis()->SetTitle("Energy in eV");
  noH->GetYaxis()->SetRangeUser(0,3100);
  noH->GetYaxis()->SetTitleOffset(1.4);
  noH->SetLineColor(2);
  noH->SetLineWidth(2);
  

  withH->Rebin(25);
  withH->GetXaxis()->SetRangeUser(7300,8400);
  withH->SetLineWidth(2);
  //withH->GetYaxis()->SetTitle("Counts/25 eV");
  //withH->GetXaxis()->SetTitle("Energy in eV");
 
  
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  

  noH->Draw();
  withH->Draw("same");
  
  //c1->Print(Form("/home/andreas/vip2/reports/1608_VIPReportLNF/sdd%d-60eVbin-7-10.png",sdd));

  //c1->Close();


}

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


void checkHistograms(){
    
    TFile *f = new TFile("/home/andreas/vip2/data/root/LNGS/1-618files-final/energyHistograms.root");
    
    TH1F *no1 = (TH1F*)f->Get("noCurrentSmallsdd1");
    TH1F *no2 = (TH1F*)f->Get("noCurrentSmallsdd2");
    TH1F *no3 = (TH1F*)f->Get("noCurrentSmallsdd3");
    TH1F *no4 = (TH1F*)f->Get("noCurrentSmallsdd4");
    TH1F *no5 = (TH1F*)f->Get("noCurrentSmallsdd5");
    TH1F *no6 = (TH1F*)f->Get("noCurrentSmallsdd6");
    
    TH1F *noSum = (TH1F*)f->Get("noCurrentSmallSum");
    
    no1->Draw();
    
    //f->Close();
    
}