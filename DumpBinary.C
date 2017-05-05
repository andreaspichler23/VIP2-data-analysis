/************************************
  DumpBinary.C
  2014/02/14 
    H.Shi
  Dump the binary data file to a 
    temporary text file. 
************************************/
#include  <stdlib.h>
#include  <stdio.h>
#include  <iostream>

#include  "TStyle.h"

#include  "ReadFullDaqBinData.h" 

using namespace std;


Int_t ReadData( TString datafilename )
{

    FILE *datafile;
    FILE *txtfile;

    Short_t word = -1;

    Short_t buff[4] = { -10 };

    // binary data for Phaes One setup in the lab: 
    TString binary_file = BIN_PATH + "/PhaseOne/" + datafilename; 
    TString txt_file    = WORK_PATH + "/work/temp.txt"; 

    if((datafile = fopen(binary_file, "r"))==NULL){
        cout << "Cannot open file: " << binary_file << endl;
    }

    if((txtfile = fopen(txt_file, "w"))==NULL){
        cout << "Cannot open file: " << txt_file << endl;
    }

    while( 1 )
    {
        if( fread( &word, sizeof(word), 1, datafile) != 1 )
        {
            if(  feof ( datafile ) ) {
                cout << " End of the file " << endl;
                break;
            }else {
                cout << " Unexpected termination ! " << endl;
                break;
            }
        }else{
            //cout << word << ", " ;

            buff[0] = buff[1]; 
            buff[1] = buff[2]; 
            buff[2] = buff[3]; 
            buff[3] = word; 

            fprintf( txtfile, "%04d, ", word );

            //if( buff[0] > 0 && buff[1] == 0 && buff[2] == 0 && buff[3] == -1 ){
            //if( buff[0] > 0 && buff[1] == -1 && buff[2] == -1 && buff[3] > 0 ){
            // After deleting the FADC module
            if( buff[0] >= 0 && buff[1] > 0 && buff[2] == -1 && buff[3] > 0 ){
                fprintf( txtfile, "\n" );
            }
        }
    }

    fclose(txtfile);
    fclose(datafile);

    return 0; 
}
