/*********************************************
            SDDStability.h 
              July 2009
              H. Shi
**********************************************/
#ifndef SDDStability_h_included_
#define SDDStability_h_included_
extern void    SDDProfile( TString filelist,  Bool_t  saveflag,  Bool_t saveplot );
extern void    SingleRunOneSDDCalib( ID nBus,  ID nSDD,   TString input_list, 
                                     Int_t nf, TString  ith_run, 
                                     Bool_t saveflag,  Bool_t plotflag, TString parfile );
extern void    MakeParRootFile( TString parfile,  TString filelist, Bool_t saveflag );
extern Int_t   GetNumberOfRuns( TString filelist ); 
extern TString GetRunNames( TString filelist,  Int_t nf );
#endif
