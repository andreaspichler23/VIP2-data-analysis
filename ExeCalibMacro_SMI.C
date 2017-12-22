#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
//#include <chrono>
//#include <thread>
#include <unistd.h>    

// for single files better use PhaseOneCalibLoop in SDDcalib.c or maybe not after all

//g++ -Wall  ExeCalibMacro_SMI.C  -o ExeCalibMacro_SMI
// ./ExeCalibMacro_SMI

using namespace std;

int main(){
    
    string fileName;
    //string rootFilePath = "/home/andreas/vip2/data/root/LNGS/test/";
    string rootFileName;
    string finalString;
    string smallQuote = "'";
    string bigQuote = "\"";
    string firstPart = "root -q -b 'CalibMacro_SMI.C(";
    string Comma = ",";
    string sddString;
    string adcString;
    string lastPart = ",\"TiMnCu\",\"smi\")'";
    
    stringstream adcSS, sddSS;
    ifstream fileList;
    int adcChannelList[6] = {0,2,3,5,6,7};
    int adcChannel;
    int startFileNumber = 1;
    int endFileNumber = 15;
    //int fileNumber;
    
    
    //fileList.open("/home/andreas/vip2/filelist/1-618Files-LNGS.txt");
    fileList.open("/home/andreas/vip2/filelist/ListofFilelists_SMI.txt");
    
  for( int fileNumber = 1; fileNumber <= endFileNumber; fileNumber++ ){
    
    fileList >> fileName;
    //cout << fileName << endl;
    rootFileName = "15files-original/" + fileName + ".root";
    //fileNumber = atoi(rootFileName.c_str());
    cout << fileNumber << endl;
    
    if ( fileNumber < startFileNumber ){ continue; }
    //cout << fileNumber << endl;
    
    
    for( int sdd = 1; sdd < 7; sdd++ ){
        
        adcChannel = adcChannelList[sdd-1];
        adcSS << adcChannel;
        adcString = adcSS.str();
        
        sddSS << sdd;
        sddString = sddSS.str();
        
        finalString = firstPart + sddString + Comma + adcString + Comma + bigQuote + rootFileName + bigQuote + lastPart;
        system(finalString.c_str());
        
        //cout << finalString << endl;
        
        adcSS.str("");
        sddSS.str("");
        
    }
    
   

  }
    
    fileList.close();
    return 0;
    
}