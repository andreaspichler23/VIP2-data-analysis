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

//g++ -Wall  ExeCalibMacro.C  -o ExeCalibMacro
// ./ExeCalibMacro

using namespace std;

int main(){
    
    string fileName;
    //string rootFilePath = "/home/andreas/vip2/data/root/LNGS/test/";
    string rootFileName;
    string finalString;
    string smallQuote = "'";
    string bigQuote = "\"";
    string firstPart = "root -q -b 'CalibMacro.C(";
    string Comma = ",";
    string sddString;
    string adcString;
    string lastPart = ",\"TiMnCu\",\"lngs\")'";
    
    stringstream adcSS, sddSS;
    ifstream fileList;
    int adcChannelList[6] = {0,2,3,5,6,7};
    int adcChannel;
    int startFileNumber = 239;
    int endFileNumber = 239;
    //int fileNumber;
    
    
    fileList.open("/home/andreas/vip2/filelist/1-618Files-LNGS.txt");
    
  for( int fileNumber = 1; fileNumber <= endFileNumber; fileNumber++ ){
    
    fileList >> fileName;
    rootFileName = "1-618files-original/" + fileName + ".root";
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
       
        
        adcSS.str("");
        sddSS.str("");
        
    }
    
   

  }
    
    fileList.close();
    return 0;
    
}