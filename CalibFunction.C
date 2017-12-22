/********************************************************** 
                TiCuCalibFunction.C 
                   25 Jun 2009
                   H. Shi
               Pile up not considered.
          Oct 2009 Pile up and tail added.
          May 2010 added TiMnCuFullFitFunction.

                 CalibFunction.C
                    16 May 2014 
                    Modified for VIP-2 use. 
***********************************************************/
#include  <stdio.h>
#include  <stdlib.h>
#include  <iostream>

#include  "TMath.h"
#include  "TString.h"

#include  "XrayLines.h"
#include  "CalibFunction.h"
//#include  "PhaseOneSDDConnectionTable.h"
//#include  "SDDclass.h"
#include  "common.h"

using namespace std;

Double_t gaussFunc(Double_t *x, Double_t *par){
    
    Double_t xx = x[0];
    Double_t gauss = par[0]/(sqrt(2*TMath::Pi())*par[2])*TMath::Exp(-((xx-par[1])*(xx-par[1]))/(2*par[2]*par[2]));
    
    return gauss;
    
}

Double_t tailFuncEnergy(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : main gauss sigma in ev
// par[1] : Main gaus mean 
// par[2] : Tail mean shift in eV
// par[3] : Fano ... unused
// par[4] : Constant noise ... unused
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    //Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    //Double_t sig   = sqrt( slope * ( par[1] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] ); // signa in channels
    Double_t sig = par[0];
    arg   = (x[0] - par[1] + par[2] ) / sig / SQRT2;
    Double_t argB1 = (x[0] - par[1] + par[2]) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t RoiCuFunc(Double_t *x, Double_t *par){
    
    // fit of the region of roi, nickel, and cu ka kb
    // this is for fitting an already scaled histogram
    
    
    Double_t xx = x[0];
    
    //par[0] = background constant
    //par[10] = background slope
    
    Double_t back = par[0] + (xx - 7000) * par[10];
    
    //par[1] = cu ka1 gain
    //par[2] = cu ka1 mean
    //par[3] = cu ka1 sigma
    
    Double_t cuKa1 = par[1]/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-par[2])*(xx-par[2]))/(2*par[3]*par[3]));
    
    Double_t cuKa2Gain = par[1] * 0.51;
    Double_t cuKa2Mean = par[2] - 19.95;
    
    Double_t cuKa2 = cuKa2Gain/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-cuKa2Mean)*(xx-cuKa2Mean))/(2*par[3]*par[3]));
    
    //par[4] = cu kb gain
    //par[5] = cu kb mean
    //par[6] = cu kb sigma
    
    Double_t cuKb = par[4]/(sqrt(2*TMath::Pi())*par[6])*TMath::Exp(-((xx-par[5])*(xx-par[5]))/(2*par[6]*par[6]));
    
    //par[7] = ni ka1 gain
    //par[8] = ni ka1 mean
    //par[9] = ni ka1 sigma
  
    
    Double_t niKa1 = par[7]/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-par[8])*(xx-par[8]))/(2*par[9]*par[9]));
    
    Double_t niKa2Gain = par[7] * 0.51;
    Double_t niKa2Mean = par[8] - 17.3;
    
    Double_t niKa2 = niKa2Gain/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-niKa2Mean)*(xx-niKa2Mean))/(2*par[9]*par[9]));
    
    // par[11] = cu ka1 tail mean shift in ev
    // par[12] = exponential beta ka
    // par[13] = tail gain ratio
    
    Double_t arg = 0.;
    Double_t sig = par[3];
    arg   = (x[0] - par[2] + par[11] ) / sig / SQRT2;
    Double_t argB1 = (x[0] - par[2] + par[11]) / sig / par[12];
    Double_t argB2 = 1./ par[12] / SQRT2;
    Double_t argB3 = 1./2. / par[12] / par[12];
    Double_t norm  = 1./2. / sig / par[12] * TMath::Exp(argB3);
    Double_t cuKaTailFunc = par[1] * par[13] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    
    // par[14] = exponential beta kb
    // par[15] = tail gain ratio cu kb
 
    Double_t BetaArg   = (x[0] - par[5] + par[11] ) / sig / SQRT2;
    Double_t BetaArgB1 = (x[0] - par[5] + par[11]) / sig / par[14];
    Double_t BetaArgB2 = 1./ par[14] / SQRT2;
    Double_t BetaArgB3 = 1./2. / par[14] / par[14];
    Double_t BetaNorm  = 1./2. / sig / par[14] * TMath::Exp(BetaArgB3);
    Double_t cuKbTailFunc = par[1] * par[15] * BetaNorm * TMath::Exp(BetaArgB1) * TMath::Erfc( BetaArg + BetaArgB2 );
    

    
    Double_t roiCuFunc = back + cuKa1 + cuKa2 + cuKb + niKa1 + niKa2 + cuKaTailFunc + cuKbTailFunc;
    
   
    return roiCuFunc; 
    
}


Double_t RoiCuFunc_NoTails(Double_t *x, Double_t *par){
    
    // fit of the region of roi, nickel, and cu ka kb
    // this is for fitting an already scaled histogram
    
    
    Double_t xx = x[0];
    
    //par[0] = background constant
    //par[10] = background slope
    
    Double_t back = par[0] + (xx - 7000) * par[10];
    
    //par[1] = cu ka1 gain
    //par[2] = cu ka1 mean
    //par[3] = cu ka1 sigma
    
    Double_t cuKa1 = par[1]/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-par[2])*(xx-par[2]))/(2*par[3]*par[3]));
    
    Double_t cuKa2Gain = par[1] * 0.51;
    Double_t cuKa2Mean = par[2] - 19.95;
    
    Double_t cuKa2 = cuKa2Gain/(sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-((xx-cuKa2Mean)*(xx-cuKa2Mean))/(2*par[3]*par[3]));
    
    //par[4] = cu kb gain
    //par[5] = cu kb mean
    //par[6] = cu kb sigma
    
    Double_t cuKb = par[4]/(sqrt(2*TMath::Pi())*par[6])*TMath::Exp(-((xx-par[5])*(xx-par[5]))/(2*par[6]*par[6]));
    
    //par[7] = ni ka1 gain
    //par[8] = ni ka1 mean
    //par[9] = ni ka1 sigma
  
    
    Double_t niKa1 = par[7]/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-par[8])*(xx-par[8]))/(2*par[9]*par[9]));
    
    Double_t niKa2Gain = par[7] * 0.51;
    Double_t niKa2Mean = par[8] - 17.3;
    
    Double_t niKa2 = niKa2Gain/(sqrt(2*TMath::Pi())*par[9])*TMath::Exp(-((xx-niKa2Mean)*(xx-niKa2Mean))/(2*par[9]*par[9]));

    
    Double_t roiCuFunc = back + cuKa1 + cuKa2 + cuKb + niKa1 + niKa2;
    
   
    return roiCuFunc; 
    
}

////////* background *//////////
//  Exponential over constant 
//  line background.   
Double_t backFunc(Double_t *x, Double_t *par){
    Double_t back = par[0] + par[1] * TMath::Exp( -par[2] * x[0] );//par[1] * x[0] + par[2] * x[0] * x[0];
    //Double_t back  = par[0] * ( x[0] - par[1] ) * ( x[0] - par[1] ) + par[2];
    return back;
}


////////* Bremsstrahlung *//////////
//  For the Ti Cu Zr source case - Same functional dependence as tail function 
//  Just the sigma of the error function is a fit parameter 
Double_t bremsFunc(Double_t *x, Double_t *par){
// par[0] : Main gaus mean zr 
// par[1] : Main gaus mean ti
// par[2] : Brems Ch - point where Erf is exactly 1/2 - kind of the bremsstrahlungs mean
// par[3] : Bremsstrahlungs sigma 
// par[4] : Bremsstrahlung Beta
// par[5] : Bremsstrahlung gain
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] )/( ZrKa1 - TiKa1 );
    arg   = (x[0] - par[2] ) / par[3] / SQRT2;
    Double_t argB1 = (x[0] - par[2] ) / par[3] / par[4];
    Double_t argB2 = 1./ par[4] / SQRT2;
    Double_t argB3 = 1./2. / par[4] / par[4];
    Double_t norm  = 1./2. / par[3] / par[4] * TMath::Exp(argB3);
    Double_t tail_func = par[5] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

//////* shelf function ////////
// Constant below peak mean channel
Double_t shelfFunc(Double_t *x, Double_t *par){

  // par[0] = Ti ka 1 channel
  // par[1] = mn or zr ka 1 channel
  // par[2] = fano
  // par[3] = const noise
  // par[4] = shelf gain
  // par[5] = mean channel of corresponding gaussian
  // par[6] = gain of the corresponding gaussian
  // par[7] = source parameter; 0 for ti mn, 1 for zr ti calibration
    Double_t shelf = 0; 
  /*  if( x[0] < par[1] ){
        shelf = par[0]; 
    }
    return shelf; */
    Double_t slope;
    if(par[7] == 0) { slope = ( par[1] - par[0] )/( MnKa1 - TiKa1 );}
    if(par[7] == 1) { slope = ( par[1] - par[0] )/( ZrKa1 - TiKa1 );}
   Double_t sig = sqrt(slope * par[2] * SiW * par[5] + par[3] * par[3]);
   
   

   shelf = par[4] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[5])/(sqrt(2)*sig) ); // 1/2 because erfc would be 2 normally on the left side

   return shelf;   

}

////////* TiKa1 Gaussian *////////
Double_t tika1Func(Double_t *x, Double_t *par){
/*
par[0] = Ti ka1 ch
1 = mn ka 1 ch OR zr ka 1 ch
2 = ti ka1 gain
3 = fano noise
4 = constant noise
5 = source parameter.. 0 for mn ti calibration, 1 for zr ti calibration
*/

    Double_t arg = 0;
    Double_t slope;
    if(par[5] == 0) { slope = ( par[1] - par[0] )/( MnKa1 - TiKa1 );}
    if(par[5] == 1) { slope = ( par[1] - par[0] )/( ZrKa1 - TiKa1 );}
    Double_t sig = sqrt (slope * par[3] * SiW * par[0] + par[4] * par[4]);
    if ( sig != 0 ){
        arg = ( x[0] - par[0] ) / sig;
    }
    Double_t tika1gauss = par[2] / sig * TMath::Exp(-0.5 * arg * arg);
    return tika1gauss;
}



////////* MnKa1 Gaussian *////////
////////  For Ti Mn source case /////
Double_t mnka1Func(Double_t *x, Double_t *par){
    Double_t arg = 0;
    Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    Double_t sig = sqrt (slope * par[3] * SiW * par[0] + par[4] * par[4]);
    if ( sig != 0 ){
        arg = ( x[0] - par[0] ) / sig;
    }
    Double_t mnka1gauss = par[2] / sig * TMath::Exp(-0.5 * arg * arg);
    return mnka1gauss;
}

////////* CuKa1 Gaussian *////////
Double_t cuka1Func(Double_t *x, Double_t *par){

    //par[0] = cu mean
    //par[1] = mn mean or zr ka1 mean
    //par[2] = ti mean
    //par[3] = cu gain
  // 4 = fano, 5 = const
// 6 = source parameter ... 0 for mn ti calibration; 1 for zr ti calibration
    Double_t arg  = 0;
    Double_t mean = par[0];
    Double_t slope;
    if(par[6] == 0) { slope = ( par[1] - par[2] )/( MnKa1 - TiKa1 );}
    if(par[6] == 1) { slope = ( par[1] - par[2] )/( ZrKa1 - TiKa1 );}
    Double_t sig = sqrt (slope * par[4] * SiW * mean + par[5] * par[5]);
    if ( sig != 0 ){
        arg = ( x[0] - mean ) / sig;
    }
    Double_t cuka1gauss = par[3] / sig * TMath::Exp(-0.5 * arg * arg);
    return cuka1gauss;
}

////////* ZrKa1 Gaussian *////////
Double_t zrka1Func(Double_t *x, Double_t *par){
/*
par[0] = zr ka1 channel
par[1] = tika1 channel
par[2] = zr gain
par[3] = fano
par[4] = const noise
*/
    Double_t arg = 0;
    Double_t slope = ( par[0] - par[1] )/( ZrKa1 - TiKa1 );
    Double_t mean  =  par[0];
    Double_t sig = sqrt (slope * par[3] * SiW * mean + par[4] * par[4]);
    if ( sig != 0 ){
        arg = ( x[0] - mean ) / sig;
    }
    Double_t zrka1gauss = par[2] / sig * TMath::Exp(-0.5 * arg * arg);
    return zrka1gauss;
}

////////* CaKa1 Gaussian *////////
Double_t caka1Func(Double_t *x, Double_t *par){
    Double_t arg = 0;
    Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    Double_t mean  =  par[0] - slope * ( MnKa1 - CaKa1 );
    Double_t sig = sqrt (slope * par[3] * SiW * mean + par[4] * par[4]);
    if ( sig != 0 ){
        arg = ( x[0] - mean ) / sig;
    }
    Double_t caka1gauss = par[2] / sig * TMath::Exp(-0.5 * arg * arg);
    return caka1gauss;
}

//////*  Tail Function - gaus, exp and erfc part *//////
Double_t catailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : Main gaus mean mn 
// par[1] : Main gaus mean ti
// par[2] : Tail mean shift
// par[3] : Fano 
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    Double_t sig   = sqrt( slope * ( par[1] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    arg   = (x[0] - par[1] + par[2] * slope ) / sig / SQRT2;
    Double_t argB1 = (x[0] - par[1] + par[2] * slope ) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t titailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : Main gaus mean mn or Zr 
// par[1] : Main gaus mean ti
// par[2] : Tail mean shift
// par[3] : Fano 
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
// par[8] : source parameter... 0 for mn ti, 1 for ti zr
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope;
    Double_t mean = par[1] - par[2] * slope;
    if(par[8] == 0) { slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );}
    if(par[8] == 1) { slope = ( par[0] - par[1] )/( ZrKa1 - TiKa1 );}
    Double_t sig   = sqrt( slope * mean * SiW * par[3] + par[4] * par[4] );
    arg   = (x[0] - mean ) / sig / SQRT2;
    Double_t argB1 = (x[0] - mean ) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t tibetatailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : mean ti kb
// par[1] : Main gaus mean ti
// par[2] : Tail mean shift
// par[3] : Fano 
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
// par[8] : ti kb to ka ratio
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] )/( MnKb1 - TiKa1 );
    Double_t sig   = sqrt( slope * ( par[1] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    arg   = (x[0] - par[1] + par[2] * slope ) / sig / SQRT2;
    Double_t argB1 = (x[0] - par[1] + par[2] * slope ) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[8] * par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t mntailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : Main gaus mean mn
// par[1] : Main gaus mean ti
// par[2] : Tail mean shift
// par[3] : Fano 
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    Double_t sig   = sqrt( slope * ( par[0] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    arg   = (x[0] - par[0] + par[2] * slope ) / sig / SQRT2;
    Double_t argB1 = (x[0] - par[0] + par[2] * slope ) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t mnbetatailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : mean mn kb
// par[1] : Main gaus mean ti
// par[2] : Tail mean shift
// par[3] : Fano 
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
// par[8] : mn kb to ka ratio
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] )/( MnKb1 - TiKa1 );
    Double_t sig   = sqrt( slope * ( par[0] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    arg   = (x[0] - par[0] + par[2] * slope ) / sig / SQRT2;
    Double_t argB1 = (x[0] - par[0] + par[2] * slope ) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[8] * par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t cutailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : Main gaus mean mn or Zr
// par[1] : Main gaus mean ti
// par[2] : Cu Ka1
// par[3] : Tail mean shift
// par[4] : Fano 
// par[5] : Constant noise
// par[6] : Cu Ka1 Gaus gain
// par[7] : Exponential beta 
// par[8] : Tail gain ratio
// par[9] : source parameter ... 0 for mn ti, 1 for zr ti calibration
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope;
    Double_t mean = par[2] - par[3] * slope;
    if(par[9] == 0) { slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );}
    if(par[9] == 1) { slope = ( par[0] - par[1] )/( ZrKa1 - TiKa1 );}
    Double_t sig   = sqrt( slope * mean * SiW * par[4] + par[5] * par[5] );
    arg   = (x[0] - mean ) / sig / SQRT2;
    Double_t argB1 = (x[0] - mean ) / sig / par[7];
    Double_t argB2 = 1./ par[7] / SQRT2;
    Double_t argB3 = 1./2. / par[7] / par[7];
    Double_t norm  = 1./2. / sig / par[7] * TMath::Exp(argB3);
    Double_t tail_func = par[6] * par[8] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}

Double_t zrtailFunc(Double_t *x, Double_t *par){
// --- tail description ---
// par[0] : Main gaus mean zr
// par[1] : Main gaus mean ti
// par[2] : Tail mean shift
// par[3] : Fano 
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Exponential beta 
// par[7] : Tail gain ratio
// par[8] : source parameter... 0 for ti mn, 1 for ti zr
//  sigma using same fano and constant noise as main peak
    Double_t arg = 0.;
    Double_t slope;
    Double_t mean = par[0] - par[2] * slope;
    if(par[9] == 0) { slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );}
    if(par[9] == 1) { slope = ( par[0] - par[1] )/( ZrKa1 - TiKa1 );}
    Double_t sig   = sqrt( slope * mean * SiW * par[3] + par[4] * par[4] );
    arg   = (x[0] - mean ) / sig / SQRT2;
    Double_t argB1 = (x[0] - mean) / sig / par[6];
    Double_t argB2 = 1./ par[6] / SQRT2;
    Double_t argB3 = 1./2. / par[6] / par[6];
    Double_t norm  = 1./2. / sig / par[6] * TMath::Exp(argB3);
    Double_t tail_func = par[5] * par[7] * norm * TMath::Exp(argB1) * TMath::Erfc( arg + argB2 );
    return tail_func;
}


////////* Pile-up function, with broadening of sigma *////////////
Double_t capileFunc(Double_t *x, Double_t *par){
// --- pile-up description ---
// par[0] : main gauss mean mn 
// par[1] : main gauss mean ti
// par[2] : pile mean shift
// par[3] : Fano
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Pile sigma factor
// par[7] : Pile up gain ratio
//  sigma broadens with respect to main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] ) / ( MnKa1 - TiKa1 );
    Double_t sig   = par[6] * sqrt( slope * ( par[1] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    if ( sig != 0 ){
        arg = ( x[0] - par[1] - par[2] * slope ) / sig;
    }
    Double_t capilegauss = par[7] * par[5] / sig * TMath::Exp(-0.5 * arg * arg);
    return   capilegauss;
}

Double_t tipileFunc(Double_t *x, Double_t *par){
// --- pile-up description ---
// par[0] : main gauss mean cu
// par[1] : main gauss mean ti
// par[2] : pile mean shift
// par[3] : Fano
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Pile sigma factor
// par[7] : Pile up gain ratio
//  sigma broadens with respect to main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] ) / ( MnKa1 - TiKa1 );
    Double_t sig   = par[6] * sqrt( slope * ( par[1] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    if ( sig != 0 ){
        arg = ( x[0] - par[1] - par[2] * slope ) / sig;
    }
    Double_t tipilegauss = par[7] * par[5] / sig * TMath::Exp(-0.5 * arg * arg);
    return   tipilegauss;
}

Double_t mnpileFunc(Double_t *x, Double_t *par){
// --- pile-up description ---
// par[0] : Main gauss mean mn
// par[1] : Main gauss mean ti
// par[2] : Pile mean shift
// par[3] : Fano
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Pile sigma factor
// par[7] : Pile up gain ratio
//  sigma broadens with respect to main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] ) / ( MnKa1 - TiKa1 );
    Double_t sig   = par[6] * sqrt( slope * ( par[0] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    if ( sig != 0 ){
        arg = ( x[0] - par[0] - par[2] * slope ) / sig;
    }
    Double_t mnpilegauss = par[7] * par[5] / sig * TMath::Exp(-0.5 * arg * arg);
    return   mnpilegauss;
}


Double_t cupileFunc(Double_t *x, Double_t *par){
// --- pile-up description ---
// par[0] : Main gauss mean cu
// par[1] : Main gauss mean ti
// par[2] : Pile mean shift
// par[3] : Fano
// par[4] : Constant noise
// par[5] : Gaus gain
// par[6] : Pile sigma factor
// par[7] : Pile up gain ratio
//  sigma broadens with respect to main peak
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] ) / ( MnKa1 - TiKa1 );
    Double_t sig   = par[6] * sqrt( slope * ( par[0] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    if ( sig != 0 ){
        arg = ( x[0] - par[0] - par[2] * slope ) / sig;
    }
    Double_t cupilegauss = par[7] * par[5] / sig * TMath::Exp(-0.5 * arg * arg);
    return   cupilegauss;
}

Double_t zrpileFunc(Double_t *x, Double_t *par){
    Double_t arg = 0.;
    Double_t slope = ( par[0] - par[1] ) / ( MnKa1 - TiKa1 );
    Double_t sig   = par[6] * sqrt( slope * ( par[1] - par[2] * slope ) * SiW * par[3] + par[4] * par[4] );
    if ( sig != 0 ){
        arg = ( x[0] - par[1] - par[2] * slope ) / sig;
    }
    Double_t zrpilegauss = par[7] * par[5] / sig * TMath::Exp(-0.5 * arg * arg);
    return   zrpilegauss;
}


//////// Escape peaks for Mn Ka1 and Ti Ka1 ////////////
Double_t tiescFunc(Double_t *x, Double_t *par){
//  --- escape peak description --- 
//  par[0] : main gauss mean mn ka1
//  par[1] : main gauss mean ti ka1
//  par[2] : ti ka1 gain 
//  par[3] : FANO 
//  par[4] : CSTN
//  par[5] : ti escape to ka1 ratio
    Double_t arg   = 0;
    Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    Double_t mean  = par[1] - SiKa * slope; 
    Double_t sig   = sqrt (slope * par[3] * SiW * mean + par[4] * par[4]);
    if ( sig != 0 ){
        arg = ( x[0] - mean ) / sig;
    }
    Double_t tiescgauss = par[5] * par[2] / sig * TMath::Exp(-0.5 * arg * arg);
    return tiescgauss;
}

Double_t mnescFunc(Double_t *x, Double_t *par){
//  --- escape peak description --- 
//  par[0] : main gauss mean mn ka1
//  par[1] : main gauss mean ti ka1
//  par[2] : ti ka1 gain 
//  par[3] : FANO 
//  par[4] : CSTN
//  par[5] : mn escape to ka1 ratio
    Double_t arg   = 0;
    Double_t slope = ( par[0] - par[1] )/( MnKa1 - TiKa1 );
    Double_t mean  = par[0] - SiKa * slope;    // for Ti Mn source case
    Double_t sig   = sqrt (slope * par[3] * SiW * mean + par[4] * par[4]);
    if ( sig != 0 ){
        arg = ( x[0] - mean ) / sig;
    }
    Double_t mnescgauss = par[5] * par[2] / sig * TMath::Exp(-0.5 * arg * arg);
    return mnescgauss;
}

//////////////////////////////////////////////////////////////////////
/*  Note : Two sources calibration function: 
Global:
    par[0] : pol1;            par[1] : pol2;           par[2]  :  pol3;
    par[3] : FANO;            par[4] : CONST_N;        par[20] :  TailBeta;
    par[21]: PileSig;         par[15]: PileShift;      par[17] :  TailShift;
    par[16]: PileGainR;
Ti:
    par[5]:  Ka1 mean         par[6]: Ka1 gain       par[7]: Kb mean
    par[8]:  KbGain/Ka1Gain   
Cu: 
    par[9]:  Ka1 mean         par[10]: Ka1 gain      par[11]: Kb mean
    par[12]: KbGain/Ka1Gain 
////////////////////////////////
    # Further for the Mn and Fe in between, apply two Gaussians:
Mn: par[13]: gainMn/Cu          
Fe: par[14]: gainFe/Cu
///////////////////////////////
*/////////////////////////////////////////////// /////////////////////
Double_t TiCuSourceFitFunc(Double_t *x, Double_t *par){
    Double_t ticu_function = 0.;
    Double_t arg[6] = {0};    
    Double_t slope        =  (par[9] - par[5])/(CuKa1 - TiKa1);    
    Double_t cuka2m        = par[9] - (CuKa1 - CuKa2) * slope;
    Double_t sig_CuKa1     = sqrt( slope * par[9] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKa2     = sqrt( slope * cuka2m * SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKb      = sqrt( slope * par[11]* SiW * par[3] + par[4] * par[4] );

    Double_t tika2m        = par[5] - (TiKa1 - TiKa2) * slope;           // TiKa2 mean channel    
    Double_t sig_TiKa1     = sqrt( slope * par[5] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa2     = sqrt( slope * tika2m * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKb      = sqrt( slope * par[7] * SiW * par[3] + par[4] * par[4] );

    if( sig_TiKa1 * sig_TiKa2 * sig_TiKb * sig_CuKa1 * sig_CuKa2 * sig_CuKb != 0){
        arg[0]  = (x[0] - par[5])/sig_TiKa1;
        arg[1]  = (x[0] - tika2m)/sig_TiKa2;
        arg[2]  = (x[0] - par[7])/sig_TiKb;
        arg[3]  = (x[0] - par[9])/sig_CuKa1;
        arg[4]  = (x[0] - cuka2m)/sig_CuKa2;
        arg[5]  = (x[0] - par[11])/sig_CuKb;
        ticu_function = backFunc(x,par)
            + par[6]   / sig_TiKa1 * TMath::Exp(-0.5 * arg[0] * arg[0])
            + TiKa2_RI / TiKa1_RI * par[6] / sig_TiKa2 * TMath::Exp(-0.5 * arg[1] * arg[1])
            + par[8]  * par[6] / sig_TiKb * TMath::Exp(-0.5 * arg[2] * arg[2])
            + par[10]  / sig_CuKa1 * TMath::Exp(-0.5 * arg[3] * arg[3])
            + CuKa2_RI / CuKa1_RI * par[10] / sig_CuKa2 * TMath::Exp(-0.5 * arg[4] * arg[4])
            + par[12]  * par[10] / sig_CuKb * TMath::Exp(-0.5 * arg[5] * arg[5]);
    }
    return ticu_function;
}


Double_t TiMnFullFitFunc(Double_t *x, Double_t *par)
{ 
    Double_t timn_function = 0.;
    Double_t arg[33] = {0.};    

   // Double_t mn_shelf = 0.;

    Double_t slope         = (par[9] - par[5])/(MnKa1 - TiKa1); // par[9] = init_mean_mn; par[5]=init_mean_ti; slope in ch / eV  

    Double_t mnka2m        = par[9] - (MnKa1 - MnKa2) * slope; // mn ka2 start value channel
    Double_t mnka1e        = par[9] - SiKa * slope; // channel with a mnka1 photon absorbed by Si -> lost energy not recorded
    Double_t mnka2e        = mnka2m - SiKa * slope; // same for mn-ka2
    Double_t mnkbe         = par[11]- SiKa * slope; // par[11]=init_mean_mnkb, same for mn-kb
   // Double_t pmnka1m       = par[9] + par[15] * slope; // par[15]=pileup mean shift in eV (60 eV), peak shift due to pulse pile up - at our rate prob = 0 
   // Double_t pmnka2m       = mnka2m + par[15] * slope; // same for mn-ka2
   // Double_t pmnkbm        = par[11]+ par[15] * slope; // same for mn-kb
    Double_t tmnka1m       = par[9] - par[14] * slope; // par[14]=Tail mean shift [eV]
    Double_t tmnka2m       = mnka2m - par[14] * slope; // same for mn-ka2
    Double_t tmnkbm        = par[11]- par[14] * slope; // same for mn-kb

    Double_t sig_MnKa1     = sqrt( slope * par[9] * SiW * par[3] + par[4] * par[4] ); // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_MnKa2     = sqrt( slope * mnka2m * SiW * par[3] + par[4] * par[4] ); 
    Double_t sig_MnKb      = sqrt( slope * par[11]* SiW * par[3] + par[4] * par[4] );
    Double_t sig_MnKa1e    = sqrt( slope * mnka1e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_MnKa2e    = sqrt( slope * mnka2e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_MnKbe     = sqrt( slope * mnkbe  * SiW * par[3] + par[4] * par[4] );
  //  Double_t sig_pMnKa1    = par[18] * sqrt( slope * pmnka1m* SiW * par[3] + par[4] * par[4] ); // par[18] is sigma broadening of pileup
  //  Double_t sig_pMnKa2    = par[18] * sqrt( slope * pmnka2m* SiW * par[3] + par[4] * par[4] );
  //  Double_t sig_pMnKb     = par[18] * sqrt( slope * pmnkbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKa1    = sqrt( slope * tmnka1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKa2    = sqrt( slope * tmnka2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKb     = sqrt( slope * tmnkbm * SiW * par[3] + par[4] * par[4] );

    Double_t tika2m        = par[5] - (TiKa1 - TiKa2) * slope;           // TiKa2 mean channel    
    Double_t tika1e        = par[5] - SiKa * slope;
    Double_t tika2e        = tika2m - SiKa * slope;
    Double_t tikbe         = par[7] - SiKa * slope; // par[7]=Ti-Kb channel
  //  Double_t ptika1m       = par[5] + par[15] * slope * 0.8;             // added 0.8 factor 23/02/2010
  //  Double_t ptika2m       = tika2m + par[15] * slope * 0.8;
  //  Double_t ptikbm        = par[7] + par[15] * slope * 0.8;
    Double_t ttika1m       = par[5] - par[14] * slope;
    Double_t ttika2m       = tika2m - par[14] * slope;
    Double_t ttikbm        = par[7] - par[14] * slope;

    Double_t sig_TiKa1     = sqrt( slope * par[5] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa2     = sqrt( slope * tika2m * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKb      = sqrt( slope * par[7] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa1e    = sqrt( slope * tika1e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa2e    = sqrt( slope * tika2e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKbe     = sqrt( slope * tikbe  * SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKa1    = par[18] * sqrt( slope * ptika1m* SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKa2    = par[18] * sqrt( slope * ptika2m* SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKb     = par[18] * sqrt( slope * ptikbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa1    = sqrt( slope * ttika1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa2    = sqrt( slope * ttika2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKb     = sqrt( slope * ttikbm * SiW * par[3] + par[4] * par[4] );

    if( sig_TiKa1 * sig_TiKa2 * sig_TiKb * sig_MnKa1 * sig_MnKa2 * sig_MnKb 
          //  * sig_CaKa1 * sig_CaKa2 * sig_CaKb 
          //  * sig_pTiKa1 * sig_pTiKa2 * sig_pTiKb  * sig_pMnKa1 * sig_pMnKa2 * sig_pMnKb  
            * sig_tTiKa1 * sig_tTiKa2 * sig_tTiKb * sig_tMnKa1 * sig_tMnKa2 * sig_tMnKb  
            * sig_MnKa1e * sig_MnKa2e * sig_MnKbe * sig_TiKa1e * sig_TiKa2e * sig_TiKbe 
          //  * sig_pCaKa1 * sig_pCaKa2 * sig_pCaKb * sig_tCaKa1 * sig_tCaKa2 * sig_tCaKb 
            != 0){
       
        arg[0]  = (x[0] - par[5])/sig_TiKa1;
        arg[1]  = (x[0] - tika2m)/sig_TiKa2;
        arg[2]  = (x[0] - par[7])/sig_TiKb;
        arg[21] = (x[0] - tika1e)/sig_TiKa1e;   // Escape peak
        arg[22] = (x[0] - tika2e)/sig_TiKa2e;
        arg[23] = (x[0] - tikbe)/sig_TiKbe;
        
        arg[3]  = (x[0] - par[9])/sig_MnKa1;
        arg[4]  = (x[0] - mnka2m)/sig_MnKa2;
        arg[5]  = (x[0] - par[11])/sig_MnKb;
        arg[6]  = (x[0] - mnka1e)/sig_MnKa1e;  // Escape peak
        arg[7]  = (x[0] - mnka2e)/sig_MnKa2e;
        arg[20] = (x[0] - mnkbe) /sig_MnKbe;
/*
        arg[24]  = (x[0] - par[22])/sig_CaKa1;
        arg[25]  = (x[0] - caka2m)/sig_CaKa2;
        arg[26]  = (x[0] - par[24])/sig_CaKb;
*//*
        arg[8]  = (x[0] - ptika1m)/sig_pTiKa1;
        arg[9]  = (x[0] - ptika2m)/sig_pTiKa2;
        arg[10] = (x[0] - ptikbm) /sig_pTiKb;
*/
        arg[11] = (x[0] - ttika1m)/sig_tTiKa1/SQRT2;
        Double_t argtika1B1 = (x[0] - ttika1m)/sig_tTiKa1/par[19]; //par[19]=tail beta slope
        Double_t argtika1B2 = 1./par[19]/SQRT2;
        Double_t argtika1B3 = 1./2./par[19]/par[19];
        Double_t normtika1  = 1./2./sig_tTiKa1/par[19]*TMath::Exp(argtika1B3);
        arg[12] = (x[0] - ttika2m)/sig_tTiKa2/SQRT2;
        Double_t argtika2B1 = (x[0] - ttika2m)/sig_tTiKa2/par[19];
        Double_t argtika2B2 = 1./par[19]/SQRT2;
        Double_t argtika2B3 = 1./2./par[19]/par[19];
        Double_t normtika2  = 1./2./sig_tTiKa2/par[19]*TMath::Exp(argtika2B3);
        arg[13] = (x[0] - ttikbm) /sig_tTiKb/SQRT2;
        Double_t argtikbB1 = (x[0] - ttikbm)/sig_tTiKb/par[19];
        Double_t argtikbB2 = 1./par[19]/SQRT2;
        Double_t argtikbB3 = 1./2./par[19]/par[19];
        Double_t normtikb  = 1./2./sig_tTiKb/par[19]*TMath::Exp(argtikbB3);
/*
        arg[14] = (x[0] - pmnka1m)/sig_pMnKa1;
        arg[15] = (x[0] - pmnka2m)/sig_pMnKa2;
        arg[16] = (x[0] - pmnkbm) /sig_pMnKb;
*/
        arg[17] = (x[0] - tmnka1m)/sig_tMnKa1/SQRT2;
        Double_t argmnka1B1 = (x[0] - tmnka1m)/sig_tMnKa1/par[20];
        Double_t argmnka1B2 = 1./par[20]/SQRT2;
        Double_t argmnka1B3 = 1./2./par[20]/par[20];
        Double_t normmnka1  = 1./2./sig_tMnKa1/par[20]*TMath::Exp(argmnka1B3);
        arg[18] = (x[0] - tmnka2m)/sig_tMnKa2/SQRT2;
        Double_t argmnka2B1 = (x[0] - tmnka2m)/sig_tMnKa2/par[20];
        Double_t argmnka2B2 = 1./par[20]/SQRT2;
        Double_t argmnka2B3 = 1./2./par[20]/par[20];
        Double_t normmnka2  = 1./2./sig_tMnKa2/par[20]*TMath::Exp(argmnka2B3);
        arg[19] = (x[0] - tmnkbm) /sig_tMnKb/SQRT2;
        Double_t argmnkbB1 = (x[0] - tmnkbm)/sig_tMnKb/par[20];
        Double_t argmnkbB2 = 1./par[20]/SQRT2;
        Double_t argmnkbB3 = 1./2./par[20]/par[20];
        Double_t normmnkb  = 1./2./sig_tMnKb/par[20]*TMath::Exp(argmnkbB3);
/*
        arg[27] = (x[0] - pcaka1m)/sig_pCaKa1;
        arg[28] = (x[0] - pcaka2m)/sig_pCaKa2;
        arg[29] = (x[0] - pcakbm) /sig_pCaKb;
        arg[30] = (x[0] - tcaka1m)/sig_tCaKa1/SQRT2;
        Double_t argcaka1B1 = (x[0] - tcaka1m)/sig_tCaKa1/par[17];
        Double_t argcaka1B2 = 1./par[17]/SQRT2;
        Double_t argcaka1B3 = 1./2./par[17]/par[17];
        Double_t normcaka1  = 1./2./sig_tCaKa1/par[17]*TMath::Exp(argcaka1B3);
        arg[31] = (x[0] - tcaka2m)/sig_tCaKa2/SQRT2;
        Double_t argcaka2B1 = (x[0] - tcaka2m)/sig_tCaKa2/par[17];
        Double_t argcaka2B2 = 1./par[17]/SQRT2;
        Double_t argcaka2B3 = 1./2./par[17]/par[17];
        Double_t normcaka2  = 1./2./sig_tCaKa2/par[17]*TMath::Exp(argcaka2B3);
        arg[32] = (x[0] - tcakbm) /sig_tCaKb/SQRT2;
        Double_t argcakbB1 = (x[0] - tcakbm)/sig_tCaKb/par[17];
        Double_t argcakbB2 = 1./par[17]/SQRT2;
        Double_t argcakbB3 = 1./2./par[17]/par[17];
        Double_t normcakb  = 1./2./sig_tCaKb/par[17]*TMath::Exp(argcakbB3);
*/

        
        timn_function = backFunc(x,par)
            + par[6]  / sig_TiKa1* TMath::Exp(-0.5 * arg[0] * arg[0]) //par[6] = ti ka1 gain   Ti ka1 gauss
            + TiKa2_RI/ TiKa1_RI * par[6] / sig_TiKa2 * TMath::Exp(-0.5 * arg[1] * arg[1]) // 100/50 for ti
            + par[8]  * par[6]   / sig_TiKb * TMath::Exp(-0.5 * arg[2] * arg[2]) // par[8]=ti kb to ka gain ratio
           // + par[16] * par[6]   / sig_pTiKa1 * TMath::Exp(-0.5 * arg[8] * arg[8]) //par[16] = pileup gain ratio
           // + par[16] * par[6]   * TiKa2_RI / TiKa1_RI / sig_pTiKa2 * TMath::Exp(-0.5 * arg[9] * arg[9]) 
           // + par[16] * par[6]   * par[8] / sig_pTiKb * TMath::Exp(-0.5 * arg[10] * arg[10])
            + par[6]  * par[15]  * normtika1 * TMath::Exp(argtika1B1) * TMath::Erfc(arg[11] + argtika1B2)// par[15]=ti tail gain ratio ka1
            + par[6]  * par[15]  * TiKa2_RI / TiKa1_RI * normtika2 * TMath::Exp(argtika2B1) * TMath::Erfc(arg[12] + argtika2B2) // ka2
            + par[6]  * par[16]  * par[8] * normtikb * TMath::Exp(argtikbB1) * TMath::Erfc(arg[13] + argtikbB2) // kb

            + par[10] / sig_MnKa1* TMath::Exp(-0.5 * arg[3] * arg[3]) // par[10] = mn ka1 gain
            + MnKa2_RI/ MnKa1_RI * par[10] / sig_MnKa2 * TMath::Exp(-0.5 * arg[4] * arg[4])
            + par[12] * par[10]  / sig_MnKb * TMath::Exp(-0.5 * arg[5] * arg[5]) // par[12]=mn kb to ka gain ratio (0.7 in one fit for lngs)
            //+ par[16] * par[10]  / sig_pMnKa1 * TMath::Exp(-0.5 * arg[14] * arg[14])
           // + par[16] * par[10]  * MnKa2_RI / MnKa1_RI / sig_pMnKa2 * TMath::Exp(-0.5 * arg[15] * arg[15])
           // + par[16] * par[10]  * par[12] / sig_pMnKb * TMath::Exp(-0.5 * arg[16] * arg[16])
            + par[10] * par[17]  * normmnka1 * TMath::Exp(argmnka1B1) * TMath::Erfc(arg[17] + argmnka1B2)
            + par[10] * par[17]  * MnKa2_RI / MnKa1_RI * normmnka2 * TMath::Exp(argmnka2B1) * TMath::Erfc(arg[18] + argmnka2B2)
            + par[10] * par[18]  * par[12] * normmnkb * TMath::Exp(argmnkbB1) * TMath::Erfc(arg[19] + argmnkbB2) 


            + par[13] * par[10]  / sig_MnKa1e * TMath::Exp(-0.5 * arg[6] * arg[6]) // Escape peak ... mn ka1
            + par[13] * par[10]  / sig_MnKa2e * TMath::Exp(-0.5 * arg[7] * arg[7]) // Escape peak ... mn k2
            + par[13] * par[10]  / sig_MnKbe  * TMath::Exp(-0.5 * arg[20]* arg[20])  // Escape peak ... mn kb

            + par[13] * par[6]  / sig_TiKa1e * TMath::Exp(-0.5 * arg[21] * arg[21]) // escape peaks ti
            + par[13] * par[6]  / sig_TiKa2e * TMath::Exp(-0.5 * arg[22] * arg[22])
            + par[13] * par[6]  / sig_TiKbe  * TMath::Exp(-0.5 * arg[23] * arg[23])

            + par[21] * par[10] * (1./2) * TMath::Erfc( (x[0]-par[9])/(SQRT2*sig_MnKa1) ) // mn ka1 shelf
            + par[21] * MnKa2_RI/ MnKa1_RI * par[10] * (1./2) * TMath::Erfc( (x[0]-mnka2m)/(SQRT2*sig_MnKa2) ) // ka2
            + par[21] * par[12] * par[10] * (1./2) * TMath::Erfc( (x[0]-par[11])/(SQRT2*sig_MnKb) ) // kb
 
            + par[21] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[5])/(SQRT2*sig_TiKa1) ) // ti ka1 shelf
            + par[21] * TiKa2_RI/ TiKa1_RI * par[6] * (1./2) * TMath::Erfc( (x[0]-tika2m)/(SQRT2*sig_TiKa2) )
            + par[21] * par[8] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[7])/(SQRT2*sig_TiKb) );

    }
    return timn_function;
}


Double_t TiCuFullFitFunc(Double_t *x, Double_t *par)
{  
    // With Fe peak included. The mean value of Fe Peak decided by Ti Cu, not free.  
    Double_t ticu_function = 0.;
    Double_t arg[20] = {0};    
    Double_t slope         = (par[9] - par[5])/(CuKa1 - TiKa1);    
    Double_t cuka2m        = par[9] - (CuKa1 - CuKa2) * slope;
    Double_t pcuka1m       = par[9] + par[15] * slope;
    Double_t pcuka2m       = cuka2m + par[15] * slope;
    Double_t pcukbm        = par[11]+ par[15] * slope;
    Double_t tcuka1m       = par[9] - par[17] * slope;
    Double_t tcuka2m       = cuka2m - par[17] * slope;
    Double_t tcukbm        = par[11]- par[17] * slope;
    Double_t sig_CuKa1     = sqrt( slope * par[9] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKa2     = sqrt( slope * cuka2m * SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKb      = sqrt( slope * par[11]* SiW * par[3] + par[4] * par[4] );
    Double_t sig_pCuKa1    = par[21] * sqrt( slope * pcuka1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_pCuKa2    = par[21] * sqrt( slope * pcuka2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_pCuKb     = par[21] * sqrt( slope * pcukbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tCuKa1    = sqrt( slope * tcuka1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tCuKa2    = sqrt( slope * tcuka2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tCuKb     = sqrt( slope * tcukbm * SiW * par[3] + par[4] * par[4] );

    Double_t tika2m        = par[5] - (TiKa1 - TiKa2) * slope;           // TiKa2 mean channel    
    Double_t ptika1m       = par[5] + par[15] * slope * 0.8;             // added 0.8 factor 23/02/2010
    Double_t ptika2m       = tika2m + par[15] * slope * 0.8;
    Double_t ptikbm        = par[7] + par[15] * slope * 0.8;
    Double_t ttika1m       = par[5] - par[17] * slope;
    Double_t ttika2m       = tika2m - par[17] * slope;
    Double_t ttikbm        = par[7] - par[17] * slope;
    Double_t sig_TiKa1     = sqrt( slope * par[5] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa2     = sqrt( slope * tika2m * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKb      = sqrt( slope * par[7] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_pTiKa1    = par[21] * sqrt( slope * ptika1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_pTiKa2    = par[21] * sqrt( slope * ptika2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_pTiKb     = par[21] * sqrt( slope * ptikbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa1    = sqrt( slope * ttika1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa2    = sqrt( slope * ttika2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKb     = sqrt( slope * ttikbm * SiW * par[3] + par[4] * par[4] );

    Double_t mnm           = (par[9]+par[5])/2 - ( (CuKa1+TiKa1)/2 - MnKa1 ) * slope; 
    Double_t sig_Mn        = sqrt( slope * mnm * SiW * par[3] + par[4] * par[4] );
    Double_t fem           = (par[9]+par[5])/2 - ( (CuKa1+TiKa1)/2 - FeKa1 ) * slope; 
    Double_t sig_Fe        = sqrt( slope * fem * SiW * par[3] + par[4] * par[4] );

    if( sig_TiKa1 * sig_TiKa2 * sig_TiKb * sig_CuKa1 * sig_CuKa2 * sig_CuKb  
            * sig_pTiKa1 * sig_pTiKa2 * sig_pTiKb * sig_pCuKa1 * sig_pCuKa2 * sig_pCuKb  
            * sig_tTiKa1 * sig_tTiKa2 * sig_tTiKb * sig_tCuKa1 * sig_tCuKa2 * sig_tCuKb  
            * sig_Mn   * sig_Fe != 0){
        arg[0]  = (x[0] - par[5])/sig_TiKa1;
        arg[1]  = (x[0] - tika2m)/sig_TiKa2;
        arg[2]  = (x[0] - par[7])/sig_TiKb;
        arg[3]  = (x[0] - par[9])/sig_CuKa1;
        arg[4]  = (x[0] - cuka2m)/sig_CuKa2;
        arg[5]  = (x[0] - par[11])/sig_CuKb;
        arg[6]  = (x[0] - mnm)/sig_Mn;
        arg[7]  = (x[0] - fem)/sig_Fe;

        arg[8]  = (x[0] - ptika1m)/sig_pTiKa1;
        arg[9]  = (x[0] - ptika2m)/sig_pTiKa2;
        arg[10] = (x[0] - ptikbm) /sig_pTiKb;
        arg[11] = (x[0] - ttika1m)/sig_tTiKa1/SQRT2;
        Double_t argtika1B1 = (x[0] - ttika1m)/sig_tTiKa1/par[20];
        Double_t argtika1B2 = 1./par[20]/SQRT2;
        Double_t argtika1B3 = 1./2./par[20]/par[20];
        Double_t normtika1  = 1./2./sig_tTiKa1/par[20]*TMath::Exp(argtika1B3);
        arg[12] = (x[0] - ttika2m)/sig_tTiKa2/SQRT2;
        Double_t argtika2B1 = (x[0] - ttika2m)/sig_tTiKa2/par[20];
        Double_t argtika2B2 = 1./par[20]/SQRT2;
        Double_t argtika2B3 = 1./2./par[20]/par[20];
        Double_t normtika2  = 1./2./sig_tTiKa2/par[20]*TMath::Exp(argtika2B3);
        arg[13] = (x[0] - ttikbm) /sig_tTiKb/SQRT2;
        Double_t argtikbB1 = (x[0] - ttikbm)/sig_tTiKb/par[20];
        Double_t argtikbB2 = 1./par[20]/SQRT2;
        Double_t argtikbB3 = 1./2./par[20]/par[20];
        Double_t normtikb  = 1./2./sig_tTiKb/par[20]*TMath::Exp(argtikbB3);

        arg[14] = (x[0] - pcuka1m)/sig_pCuKa1;
        arg[15] = (x[0] - pcuka2m)/sig_pCuKa2;
        arg[16] = (x[0] - pcukbm) /sig_pCuKb;
        arg[17] = (x[0] - tcuka1m)/sig_tCuKa1/SQRT2;
        Double_t argcuka1B1 = (x[0] - tcuka1m)/sig_tCuKa1/par[20];
        Double_t argcuka1B2 = 1./par[20]/SQRT2;
        Double_t argcuka1B3 = 1./2./par[20]/par[20];
        Double_t normcuka1  = 1./2./sig_tCuKa1/par[20]*TMath::Exp(argcuka1B3);
        arg[18] = (x[0] - tcuka2m)/sig_tCuKa2/SQRT2;
        Double_t argcuka2B1 = (x[0] - tcuka2m)/sig_tCuKa2/par[20];
        Double_t argcuka2B2 = 1./par[20]/SQRT2;
        Double_t argcuka2B3 = 1./2./par[20]/par[20];
        Double_t normcuka2  = 1./2./sig_tCuKa2/par[20]*TMath::Exp(argcuka2B3);
        arg[19] = (x[0] - tcukbm) /sig_tCuKb/SQRT2;
        Double_t argcukbB1 = (x[0] - tcukbm)/sig_tCuKb/par[20];
        Double_t argcukbB2 = 1./par[20]/SQRT2;
        Double_t argcukbB3 = 1./2./par[20]/par[20];
        Double_t normcukb  = 1./2./sig_tCuKb/par[20]*TMath::Exp(argcukbB3);

        
        ticu_function = backFunc(x,par)                            // Function containing tails for each peak
            + par[6]  / sig_TiKa1* TMath::Exp(-0.5 * arg[0] * arg[0])
            + TiKa2_RI/ TiKa1_RI * par[6] / sig_TiKa2 * TMath::Exp(-0.5 * arg[1] * arg[1])
            + par[8]  * par[6]   / sig_TiKb * TMath::Exp(-0.5 * arg[2] * arg[2])
            + par[16] * par[6]   / sig_pTiKa1 * TMath::Exp(-0.5 * arg[8] * arg[8])
            + par[16] * par[6]   * TiKa2_RI / TiKa1_RI / sig_pTiKa2 * TMath::Exp(-0.5 * arg[9] * arg[9])
            + par[16] * par[6]   * par[8] / sig_pTiKb * TMath::Exp(-0.5 * arg[10] * arg[10])
            + par[6]  * par[18]  * normtika1 * TMath::Exp(argtika1B1) * TMath::Erfc(arg[11] + argtika1B2)
            + par[6]  * par[18]  * TiKa2_RI / TiKa1_RI * normtika2 * TMath::Exp(argtika2B1) * TMath::Erfc(arg[12] + argtika2B2)
            + par[6]  * par[18]  * par[8] * normtikb * TMath::Exp(argtikbB1) * TMath::Erfc(arg[13] + argtikbB2)

            + par[10] / sig_CuKa1* TMath::Exp(-0.5 * arg[3] * arg[3])
            + CuKa2_RI/ CuKa1_RI * par[10] / sig_CuKa2 * TMath::Exp(-0.5 * arg[4] * arg[4])
            + par[12] * par[10]  / sig_CuKb * TMath::Exp(-0.5 * arg[5] * arg[5])
            + par[16] * par[10]  / sig_pCuKa1 * TMath::Exp(-0.5 * arg[14] * arg[14])
            + par[16] * par[10]  * CuKa2_RI / CuKa1_RI / sig_pCuKa2 * TMath::Exp(-0.5 * arg[15] * arg[15])
            + par[16] * par[10]  * par[12] / sig_pCuKb * TMath::Exp(-0.5 * arg[16] * arg[16])
            + par[10] * par[19]  * normcuka1 * TMath::Exp(argcuka1B1) * TMath::Erfc(arg[17] + argcuka1B2)
            + par[10] * par[19]  * CuKa2_RI / CuKa1_RI * normcuka2 * TMath::Exp(argcuka2B1) * TMath::Erfc(arg[18] + argcuka2B2)
            + par[10] * par[19]  * par[12] * normcukb * TMath::Exp(argcukbB1) * TMath::Erfc(arg[19] + argcukbB2)

            + par[13] * par[10]  / sig_Mn * TMath::Exp(-0.5 * arg[6] * arg[6])
            + par[14] * par[10]  / sig_Fe * TMath::Exp(-0.5 * arg[7] * arg[7]);
    }
    return ticu_function;
}


Double_t TiMnCuFullFitFunc(Double_t *x, Double_t *par)
{   
    Double_t timncu_function = 0.;
    Double_t arg[33] = {0.};    


    Double_t slope         = (par[9] - par[5])/(MnKa1 - TiKa1); // par[9] = init_mean_mn; par[5]=init_mean_ti; slope in ch / eV  

    Double_t mnka2m        = par[9] - (MnKa1 - MnKa2) * slope; // mn ka2 start value channel
    Double_t mnka1e        = par[9] - SiKa * slope; // channel with a mnka1 photon absorbed by Si -> lost energy not recorded
    Double_t mnka2e        = mnka2m - SiKa * slope; // same for mn-ka2
    Double_t mnkbe         = par[11]- SiKa * slope; // par[11]=init_mean_mnkb, same for mn-kb
   // Double_t pmnka1m       = par[9] + par[15] * slope; // par[15]=pileup mean shift in eV (60 eV), peak shift due to pulse pile up - at our rate prob = 0 
   // Double_t pmnka2m       = mnka2m + par[15] * slope; // same for mn-ka2
   // Double_t pmnkbm        = par[11]+ par[15] * slope; // same for mn-kb
    Double_t tmnka1m       = par[9] - par[14] * slope; // par[14]=Tail mean shift [eV]
    Double_t tmnka2m       = mnka2m - par[14] * slope; // same for mn-ka2
    Double_t tmnkbm        = par[11]- par[14] * slope; // same for mn-kb

    Double_t sig_MnKa1     = sqrt( slope * par[9] * SiW * par[3] + par[4] * par[4] ); // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_MnKa2     = sqrt( slope * mnka2m * SiW * par[3] + par[4] * par[4] ); 
    Double_t sig_MnKb      = sqrt( slope * par[11]* SiW * par[3] + par[4] * par[4] );
    Double_t sig_MnKa1e    = sqrt( slope * mnka1e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_MnKa2e    = sqrt( slope * mnka2e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_MnKbe     = sqrt( slope * mnkbe  * SiW * par[3] + par[4] * par[4] );
  //  Double_t sig_pMnKa1    = par[18] * sqrt( slope * pmnka1m* SiW * par[3] + par[4] * par[4] ); // par[18] is sigma broadening of pileup
  //  Double_t sig_pMnKa2    = par[18] * sqrt( slope * pmnka2m* SiW * par[3] + par[4] * par[4] );
  //  Double_t sig_pMnKb     = par[18] * sqrt( slope * pmnkbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKa1    = sqrt( slope * tmnka1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKa2    = sqrt( slope * tmnka2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKb     = sqrt( slope * tmnkbm * SiW * par[3] + par[4] * par[4] );

    Double_t tika2m        = par[5] - (TiKa1 - TiKa2) * slope;           // TiKa2 mean channel    
    Double_t tika1e        = par[5] - SiKa * slope;
    Double_t tika2e        = tika2m - SiKa * slope;
    Double_t tikbe         = par[7] - SiKa * slope; // par[7]=Ti-Kb channel
  //  Double_t ptika1m       = par[5] + par[15] * slope * 0.8;             // added 0.8 factor 23/02/2010
  //  Double_t ptika2m       = tika2m + par[15] * slope * 0.8;
  //  Double_t ptikbm        = par[7] + par[15] * slope * 0.8;
    Double_t ttika1m       = par[5] - par[14] * slope;
    Double_t ttika2m       = tika2m - par[14] * slope;
    Double_t ttikbm        = par[7] - par[14] * slope;

    Double_t sig_TiKa1     = sqrt( slope * par[5] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa2     = sqrt( slope * tika2m * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKb      = sqrt( slope * par[7] * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa1e    = sqrt( slope * tika1e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKa2e    = sqrt( slope * tika2e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_TiKbe     = sqrt( slope * tikbe  * SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKa1    = par[18] * sqrt( slope * ptika1m* SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKa2    = par[18] * sqrt( slope * ptika2m* SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKb     = par[18] * sqrt( slope * ptikbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa1    = sqrt( slope * ttika1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa2    = sqrt( slope * ttika2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKb     = sqrt( slope * ttikbm * SiW * par[3] + par[4] * par[4] );

    Double_t cuka2m        = par[22] - (CuKa1 - CuKa2) * slope; // cu ka2 start value channel
    Double_t cuka1e        = par[22] - SiKa * slope; // channel with a cuka1 photon absorbed by Si -> lost energy not recorded
    Double_t cuka2e        = cuka2m - SiKa * slope; // same for cu-ka2
    Double_t cukbe         = par[24]- SiKa * slope; // par[24]=init_mean_cukb, same for cu-kb
    Double_t tcuka1m       = par[22] - par[14] * slope; // par[14]=Tail mean shift [eV]
    Double_t tcuka2m       = mnka2m - par[14] * slope; // same for cu-ka2
    Double_t tcukbm        = par[24]- par[14] * slope; // same for cu-kb

    Double_t sig_CuKa1     = sqrt( slope * par[22] * SiW * par[3] + par[4] * par[4] ); // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_CuKa2     = sqrt( slope * cuka2m * SiW * par[3] + par[4] * par[4] ); 
    Double_t sig_CuKb      = sqrt( slope * par[24]* SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKa1e    = sqrt( slope * cuka1e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKa2e    = sqrt( slope * cuka2e * SiW * par[3] + par[4] * par[4] );
    Double_t sig_CuKbe     = sqrt( slope * cukbe  * SiW * par[3] + par[4] * par[4] );

    Double_t sig_tCuKa1    = sqrt( slope * tcuka1m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tCuKa2    = sqrt( slope * tcuka2m* SiW * par[3] + par[4] * par[4] );
    Double_t sig_tCuKb     = sqrt( slope * tcukbm * SiW * par[3] + par[4] * par[4] );

    if( sig_TiKa1 * sig_TiKa2 * sig_TiKb * sig_MnKa1 * sig_MnKa2 * sig_MnKb 
            * sig_CuKa1 * sig_CuKa2 * sig_CuKb 
            * sig_tTiKa1 * sig_tTiKa2 * sig_tTiKb * sig_tMnKa1 * sig_tMnKa2 * sig_tMnKb  
            * sig_MnKa1e * sig_MnKa2e * sig_MnKbe * sig_TiKa1e * sig_TiKa2e * sig_TiKbe 
            != 0){
       
        arg[0]  = (x[0] - par[5])/sig_TiKa1;
        arg[1]  = (x[0] - tika2m)/sig_TiKa2;
        arg[2]  = (x[0] - par[7])/sig_TiKb;
        arg[21] = (x[0] - tika1e)/sig_TiKa1e;   // Escape peak
        arg[22] = (x[0] - tika2e)/sig_TiKa2e;
        arg[23] = (x[0] - tikbe)/sig_TiKbe;
        
        arg[3]  = (x[0] - par[9])/sig_MnKa1;
        arg[4]  = (x[0] - mnka2m)/sig_MnKa2;
        arg[5]  = (x[0] - par[11])/sig_MnKb;
        arg[6]  = (x[0] - mnka1e)/sig_MnKa1e;  // Escape peak
        arg[7]  = (x[0] - mnka2e)/sig_MnKa2e;
        arg[20] = (x[0] - mnkbe) /sig_MnKbe;

        arg[24]  = (x[0] - par[22])/sig_CuKa1; // cu ka1
        arg[25]  = (x[0] - cuka2m)/sig_CuKa2;
        arg[26]  = (x[0] - par[24])/sig_CuKb;
        arg[27]  = (x[0] - cuka1e)/sig_CuKa1e;  // Escape peak
        arg[28]  = (x[0] - cuka2e)/sig_CuKa2e;
        arg[29] = (x[0] - cukbe) /sig_CuKbe;

        arg[11] = (x[0] - ttika1m)/sig_tTiKa1/SQRT2;
        Double_t argtika1B1 = (x[0] - ttika1m)/sig_tTiKa1/par[19]; //par[19]=tail beta slope
        Double_t argtika1B2 = 1./par[19]/SQRT2;
        Double_t argtika1B3 = 1./2./par[19]/par[19];
        Double_t normtika1  = 1./2./sig_tTiKa1/par[19]*TMath::Exp(argtika1B3);
        arg[12] = (x[0] - ttika2m)/sig_tTiKa2/SQRT2;
        Double_t argtika2B1 = (x[0] - ttika2m)/sig_tTiKa2/par[19];
        Double_t argtika2B2 = 1./par[19]/SQRT2;
        Double_t argtika2B3 = 1./2./par[19]/par[19];
        Double_t normtika2  = 1./2./sig_tTiKa2/par[19]*TMath::Exp(argtika2B3);
        arg[13] = (x[0] - ttikbm) /sig_tTiKb/SQRT2;
        Double_t argtikbB1 = (x[0] - ttikbm)/sig_tTiKb/par[19];
        Double_t argtikbB2 = 1./par[19]/SQRT2;
        Double_t argtikbB3 = 1./2./par[19]/par[19];
        Double_t normtikb  = 1./2./sig_tTiKb/par[19]*TMath::Exp(argtikbB3);
/*
        arg[14] = (x[0] - pmnka1m)/sig_pMnKa1;
        arg[15] = (x[0] - pmnka2m)/sig_pMnKa2;
        arg[16] = (x[0] - pmnkbm) /sig_pMnKb;
*/
        arg[17] = (x[0] - tmnka1m)/sig_tMnKa1/SQRT2;
        Double_t argmnka1B1 = (x[0] - tmnka1m)/sig_tMnKa1/par[20];
        Double_t argmnka1B2 = 1./par[20]/SQRT2;
        Double_t argmnka1B3 = 1./2./par[20]/par[20];
        Double_t normmnka1  = 1./2./sig_tMnKa1/par[20]*TMath::Exp(argmnka1B3);
        arg[18] = (x[0] - tmnka2m)/sig_tMnKa2/SQRT2;
        Double_t argmnka2B1 = (x[0] - tmnka2m)/sig_tMnKa2/par[20];
        Double_t argmnka2B2 = 1./par[20]/SQRT2;
        Double_t argmnka2B3 = 1./2./par[20]/par[20];
        Double_t normmnka2  = 1./2./sig_tMnKa2/par[20]*TMath::Exp(argmnka2B3);
        arg[19] = (x[0] - tmnkbm) /sig_tMnKb/SQRT2;
        Double_t argmnkbB1 = (x[0] - tmnkbm)/sig_tMnKb/par[20];
        Double_t argmnkbB2 = 1./par[20]/SQRT2;
        Double_t argmnkbB3 = 1./2./par[20]/par[20];
        Double_t normmnkb  = 1./2./sig_tMnKb/par[20]*TMath::Exp(argmnkbB3);


        arg[30] = (x[0] - tcuka1m)/sig_tCuKa1/SQRT2;
        Double_t argcuka1B1 = (x[0] - tcuka1m)/sig_tCuKa1/par[20];
        Double_t argcuka1B2 = 1./par[20]/SQRT2;
        Double_t argcuka1B3 = 1./2./par[20]/par[20];
        Double_t normcuka1  = 1./2./sig_tCuKa1/par[20]*TMath::Exp(argcuka1B3);
        arg[31] = (x[0] - tcuka2m)/sig_tCuKa2/SQRT2;
        Double_t argcuka2B1 = (x[0] - tcuka2m)/sig_tCuKa2/par[20];
        Double_t argcuka2B2 = 1./par[20]/SQRT2;
        Double_t argcuka2B3 = 1./2./par[20]/par[20];
        Double_t normcuka2  = 1./2./sig_tCuKa2/par[20]*TMath::Exp(argcuka2B3);
        arg[32] = (x[0] - tcukbm) /sig_tCuKb/SQRT2;
        Double_t argcukbB1 = (x[0] - tcukbm)/sig_tCuKb/par[20];
        Double_t argcukbB2 = 1./par[20]/SQRT2;
        Double_t argcukbB3 = 1./2./par[20]/par[20];
        Double_t normcukb  = 1./2./sig_tCuKb/par[20]*TMath::Exp(argcukbB3);


        
        timncu_function = backFunc(x,par)
            + par[6]  / sig_TiKa1* TMath::Exp(-0.5 * arg[0] * arg[0]) //par[6] = ti ka1 gain   Ti ka1 gauss
            + TiKa2_RI/ TiKa1_RI * par[6] / sig_TiKa2 * TMath::Exp(-0.5 * arg[1] * arg[1]) // 100/50 for ti
            + par[8]  * par[6]   / sig_TiKb * TMath::Exp(-0.5 * arg[2] * arg[2]) // par[8]=ti kb to ka gain ratio
            + par[6]  * par[15]  * normtika1 * TMath::Exp(argtika1B1) * TMath::Erfc(arg[11] + argtika1B2)// par[15]=ti tail gain ratio ka1
            + par[6]  * par[15]  * TiKa2_RI / TiKa1_RI * normtika2 * TMath::Exp(argtika2B1) * TMath::Erfc(arg[12] + argtika2B2) // ka2
            + par[6]  * par[16]  * par[8] * normtikb * TMath::Exp(argtikbB1) * TMath::Erfc(arg[13] + argtikbB2) // kb

            + par[10] / sig_MnKa1* TMath::Exp(-0.5 * arg[3] * arg[3]) // par[10] = mn ka1 gain
            + MnKa2_RI/ MnKa1_RI * par[10] / sig_MnKa2 * TMath::Exp(-0.5 * arg[4] * arg[4])
            + par[12] * par[10]  / sig_MnKb * TMath::Exp(-0.5 * arg[5] * arg[5]) // par[12]=mn kb to ka gain ratio (0.7 in one fit for lngs)
            + par[10] * par[17]  * normmnka1 * TMath::Exp(argmnka1B1) * TMath::Erfc(arg[17] + argmnka1B2)
            + par[10] * par[17]  * MnKa2_RI / MnKa1_RI * normmnka2 * TMath::Exp(argmnka2B1) * TMath::Erfc(arg[18] + argmnka2B2)
            + par[10] * par[18]  * par[12] * normmnkb * TMath::Exp(argmnkbB1) * TMath::Erfc(arg[19] + argmnkbB2) 

            + par[23] / sig_CuKa1* TMath::Exp(-0.5 * arg[24] * arg[24]) // cu ka1 peak
            + CuKa2_RI/ CuKa1_RI * par[23] / sig_CuKa2 * TMath::Exp(-0.5 * arg[25] * arg[25]) // cu ka2
            + par[25] * par[23]  / sig_CuKb * TMath::Exp(-0.5 * arg[26] * arg[26]) // par[12]=mn kb to ka gain ratio (0.7 in one fit for lngs) - cu kb
            + par[23] * par[26]  * normcuka1 * TMath::Exp(argcuka1B1) * TMath::Erfc(arg[30] + argcuka1B2) // cu ka1 tail
            + par[23] * par[26]  * CuKa2_RI / CuKa1_RI * normcuka2 * TMath::Exp(argcuka2B1) * TMath::Erfc(arg[31] + argcuka2B2) // cu ka2 tail
            + par[23] * par[25]  * par[26] * normcukb * TMath::Exp(argcukbB1) * TMath::Erfc(arg[32] + argcukbB2) // cu kb tail

            + par[13] * par[10]  / sig_MnKa1e * TMath::Exp(-0.5 * arg[6] * arg[6]) // Escape peak ... mn ka1
            + par[13] * par[10]  / sig_MnKa2e * TMath::Exp(-0.5 * arg[7] * arg[7]) // mn ka 2
            + par[13] * par[10]  / sig_MnKbe  * TMath::Exp(-0.5 * arg[20]* arg[20]) // mn kb


            + par[13] * par[6]  / sig_TiKa1e * TMath::Exp(-0.5 * arg[21] * arg[21]) // escape peak ti ka1
            + par[13] * par[6]  / sig_TiKa2e * TMath::Exp(-0.5 * arg[22] * arg[22]) // ka2
            + par[13] * par[6]  / sig_TiKbe  * TMath::Exp(-0.5 * arg[23] * arg[23]) // kb

            + par[21] * par[10] * (1./2) * TMath::Erfc( (x[0]-par[9])/(SQRT2*sig_MnKa1) ) // mn ka1 shelf
            + par[21] * MnKa2_RI/ MnKa1_RI * par[10] * (1./2) * TMath::Erfc( (x[0]-mnka2m)/(SQRT2*sig_MnKa2) ) // ka2
            + par[21] * par[12] * par[10] * (1./2) * TMath::Erfc( (x[0]-par[11])/(SQRT2*sig_MnKb) ) // kb
 
            + par[21] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[5])/(SQRT2*sig_TiKa1) ) // ti ka1 shelf
            + par[21] * TiKa2_RI/ TiKa1_RI * par[6] * (1./2) * TMath::Erfc( (x[0]-tika2m)/(SQRT2*sig_TiKa2) )
            + par[21] * par[8] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[7])/(SQRT2*sig_TiKb) );

    }
    return timncu_function;
}

Double_t TiCuZrFullFitFunc(Double_t *x, Double_t *par)
{   
    Double_t ticuzr_function = 0.;
    Double_t arg[33] = {0.};    


    Double_t slope         = (par[10] - par[2])/(ZrKa1 - TiKa1); // par[10] = init_mean_zr; par[2]=init_mean_ti; slope in ch / eV 
// slope calculated from zirconium and titanium atm 


    Double_t zrka2m        = par[10] - (ZrKa1 - ZrKa2) * slope; // zr ka2 start value channel
    Double_t zrkb2m        = par[12] - (ZrKb1 - ZrKb2) * slope; // zr kb2 start value channel
    Double_t zrkb3m        = par[12] - (ZrKb1 - ZrKb3) * slope; // zr kb3 start value channel

  /*   Double_t zrka1e        = par[11] - SiKa * slope; // channel with a zrka1 photon absorbed by Si -> lost energy not recorded
    Double_t zrka2e        = zrka2m - SiKa * slope; // same for mn-ka2
    Double_t zrkbe         = par[13]- SiKa * slope; // par[13]=init_mean_zrkb, same for zr-kb*/

    Double_t tzrka1m        = par[10] - par[14] * slope; // par[14]=Tail mean shift [eV] ... tail mean of zr ka1
    Double_t tzrka2m        = zrka2m - par[14] * slope; // same for zr-ka2
    Double_t tzrkb1m        = par[12]- par[14] * slope; // same for zr-kb1
    Double_t tzrkb2m        = zrkb2m- par[14] * slope; // same for zr-kb2
    Double_t tzrkb3m        = zrkb3m- par[14] * slope; // same for zr-kb3

 //   Double_t sig_brems     = sqrt( slope * par[18] * SiW * par[0] + par[1] * par[1] );  // sigma for bremsstrahlung - peak ... sigma in channel

    Double_t sig_ZrKa1     = sqrt( slope * par[10] * SiW * par[0] + par[1] * par[1] ); // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_ZrKa2     = sqrt( slope * zrka2m * SiW * par[0] + par[1] * par[1] ); 
    Double_t sig_ZrKb1     = sqrt( slope * par[12]* SiW * par[0] + par[1] * par[1] );
    Double_t sig_ZrKb2     = sqrt( slope * zrkb2m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_ZrKb3     = sqrt( slope * zrkb3m* SiW * par[0] + par[1] * par[1] );

 /*   Double_t sig_ZrKa1e    = sqrt( slope * zrka1e * SiW * par[0] + par[1] * par[1] );
    Double_t sig_ZrKa2e    = sqrt( slope * zrka2e * SiW * par[0] + par[1] * par[1] );
    Double_t sig_ZrKbe     = sqrt( slope * zrkbe  * SiW * par[0] + par[1] * par[1] );*/
    Double_t sig_tZrKa1    = sqrt( slope * tzrka1m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_tZrKa2    = sqrt( slope * tzrka2m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_tZrKb1     = sqrt( slope * tzrkb1m * SiW * par[0] + par[1] * par[1] );
    Double_t sig_tZrKb2     = sqrt( slope * tzrkb2m * SiW * par[0] + par[1] * par[1] );
    Double_t sig_tZrKb3     = sqrt( slope * tzrkb3m * SiW * par[0] + par[1] * par[1] );

    Double_t tika2m        = par[2] - (TiKa1 - TiKa2) * slope;           // TiKa2 mean channel  par[2] = Tika1 mean ch  

    Double_t ttika1m       = par[2] - par[22] * slope;// ti ka1 tail mean
    Double_t ttika2m       = tika2m - par[22] * slope;
    Double_t ttikbm        = par[4] - par[22] * slope;

    Double_t sig_TiKa1     = sqrt( slope * par[2] * SiW * par[0] + par[1] * par[1] );
    Double_t sig_TiKa2     = sqrt( slope * tika2m * SiW * par[0] + par[1] * par[1] );
    Double_t sig_TiKb      = sqrt( slope * par[4] * SiW * par[0] + par[1] * par[1] );
/*    Double_t sig_TiKa1e    = sqrt( slope * tika1e * SiW * par[0] + par[1] * par[1] );
    Double_t sig_TiKa2e    = sqrt( slope * tika2e * SiW * par[0] + par[1] * par[1] );
    Double_t sig_TiKbe     = sqrt( slope * tikbe  * SiW * par[0] + par[1] * par[1] );*/
    Double_t sig_tTiKa1    = sqrt( slope * ttika1m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_tTiKa2    = sqrt( slope * ttika2m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_tTiKb     = sqrt( slope * ttikbm * SiW * par[0] + par[1] * par[1] );

    Double_t cuka2m        = par[6] - (CuKa1 - CuKa2) * slope; // cu ka2 start value channel par[7] = cu ka1 
 /*   Double_t cuka1e        = par[7] - SiKa * slope; // channel with a cuka1 photon absorbed by Si -> lost energy not recorded
    Double_t cuka2e        = cuka2m - SiKa * slope; // same for cu-ka2
    Double_t cukbe         = par[29]- SiKa * slope; // par[24]=init_mean_cukb, same for cu-kb*/
    Double_t tcuka1m       = par[6] - par[23] * slope; // par[14]=Tail mean shift [eV] ... mean for tail of cu ka1 peak
    Double_t tcuka2m       = cuka2m - par[23] * slope; // same for cu-ka2
    Double_t tcukbm        = par[8] - par[23] * slope; // same for cu-kb

    Double_t sig_CuKa1     = sqrt( slope * par[6] * SiW * par[0] + par[1] * par[1] ); // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_CuKa2     = sqrt( slope * cuka2m * SiW * par[0] + par[1] * par[1] ); 
    Double_t sig_CuKb      = sqrt( slope * par[8]* SiW * par[0] + par[1] * par[1] );
 /*   Double_t sig_CuKa1e    = sqrt( slope * cuka1e * SiW * par[0] + par[1] * par[1] );
    Double_t sig_CuKa2e    = sqrt( slope * cuka2e * SiW * par[0] + par[1] * par[1] );
    Double_t sig_CuKbe     = sqrt( slope * cukbe  * SiW * par[0] + par[1] * par[1] );*/
    Double_t sig_tCuKa1    = sqrt( slope * tcuka1m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_tCuKa2    = sqrt( slope * tcuka2m* SiW * par[0] + par[1] * par[1] );
    Double_t sig_tCuKb     = sqrt( slope * tcukbm * SiW * par[0] + par[1] * par[1] );

    if( sig_TiKa1 * sig_TiKa2 * sig_TiKb * sig_ZrKa1 * sig_ZrKa2 * sig_ZrKb1 *sig_ZrKb2 * sig_ZrKb3 
            * sig_CuKa1 * sig_CuKa2 * sig_CuKb 
            * sig_tTiKa1 * sig_tTiKa2 * sig_tTiKb
 	    //* sig_tZrKa1 * sig_tZrKa2 * sig_tZrKb 
            //* sig_ZrKa1e * sig_ZrKa2e * sig_ZrKbe
	    //* sig_TiKa1e * sig_TiKa2e * sig_TiKbe 
            != 0){

 // ----------- Titanium ------------------
      
        arg[0]  = (x[0] - par[2])/sig_TiKa1; // ti main peak gauss arguments
        arg[1]  = (x[0] - tika2m)/sig_TiKa2;
        arg[2]  = (x[0] - par[4])/sig_TiKb;

        arg[14] = (x[0] - ttika1m)/sig_tTiKa1/SQRT2; // par[20] = tail beta slope for titanium ... for ti ka1 tail
        Double_t argtika1B1 = (x[0] - ttika1m)/sig_tTiKa1/par[20];
        Double_t argtika1B2 = 1./par[20]/SQRT2;
        Double_t argtika1B3 = 1./2./par[20]/par[20];
        Double_t normtika1  = 1./2./sig_tTiKa1/par[20]*TMath::Exp(argtika1B3);

        arg[15] = (x[0] - ttika2m)/sig_tTiKa2/SQRT2; // for ti ka2 tail
        Double_t argtika2B1 = (x[0] - ttika2m)/sig_tTiKa2/par[20]; 
        Double_t argtika2B2 = 1./par[20]/SQRT2;
        Double_t argtika2B3 = 1./2./par[20]/par[20];
        Double_t normtika2  = 1./2./sig_tTiKa2/par[20]*TMath::Exp(argtika2B3);

        arg[16] = (x[0] - ttikbm)/sig_tTiKb/SQRT2; // for ti kb1 tail
        Double_t argtikbB1 = (x[0] - ttikbm)/sig_tTiKb/par[20]; 
        Double_t argtikbB2 = 1./par[20]/SQRT2;
        Double_t argtikbB3 = 1./2./par[20]/par[20];
        Double_t normtikb  = 1./2./sig_tTiKb/par[20]*TMath::Exp(argtikbB3);
   
// --------- Cupper ---------------------

        arg[24]  = (x[0] - par[6])/sig_CuKa1; // cu ka1
        arg[25]  = (x[0] - cuka2m)/sig_CuKa2;
        arg[26]  = (x[0] - par[8])/sig_CuKb;

        arg[17] = (x[0] - tcuka1m)/sig_tCuKa1/SQRT2; // par[21] = tail beta slope for copper ... for cu ka1 tail
        Double_t argcuka1B1 = (x[0] - tcuka1m)/sig_tCuKa1/par[21];
        Double_t argcuka1B2 = 1./par[21]/SQRT2;
        Double_t argcuka1B3 = 1./2./par[21]/par[21];
        Double_t normcuka1  = 1./2./sig_tCuKa1/par[21]*TMath::Exp(argcuka1B3);

        arg[18] = (x[0] - tcuka2m)/sig_tCuKa2/SQRT2; // for cu ka2 tail
        Double_t argcuka2B1 = (x[0] - tcuka2m)/sig_tCuKa2/par[21]; 
        Double_t argcuka2B2 = 1./par[21]/SQRT2;
        Double_t argcuka2B3 = 1./2./par[21]/par[21];
        Double_t normcuka2  = 1./2./sig_tCuKa2/par[21]*TMath::Exp(argcuka2B3);

        arg[19] = (x[0] - tcukbm)/sig_tCuKb/SQRT2; // for cu kb1 tail
        Double_t argcukbB1 = (x[0] - tcukbm)/sig_tCuKb/par[21]; 
        Double_t argcukbB2 = 1./par[21]/SQRT2;
        Double_t argcukbB3 = 1./2./par[21]/par[21];
        Double_t normcukb  = 1./2./sig_tCuKb/par[21]*TMath::Exp(argcukbB3);


// -------------- Zirconium --------------------

        arg[3]  = (x[0] - par[10])/sig_ZrKa1; // Zr main peak gauss arguments 
        arg[4]  = (x[0] - zrka2m)/sig_ZrKa2;
        arg[5]  = (x[0] - par[12])/sig_ZrKb1;
        arg[6]  = (x[0] - zrkb2m)/sig_ZrKb2;
        arg[7]  = (x[0] - zrkb3m)/sig_ZrKb3;

        arg[8] = (x[0] - tzrka1m)/sig_tZrKa1/SQRT2; // par[16] = tail beta slope for zirconium ... for zr ka1 tail
        Double_t argzrka1B1 = (x[0] - tzrka1m)/sig_tZrKa1/par[16];
        Double_t argzrka1B2 = 1./par[16]/SQRT2;
        Double_t argzrka1B3 = 1./2./par[16]/par[16];
        Double_t normzrka1  = 1./2./sig_tZrKa1/par[16]*TMath::Exp(argzrka1B3);

        arg[9] = (x[0] - tzrka2m)/sig_tZrKa2/SQRT2; // for zr ka2 tail
        Double_t argzrka2B1 = (x[0] - tzrka2m)/sig_tZrKa2/par[16]; 
        Double_t argzrka2B2 = 1./par[16]/SQRT2;
        Double_t argzrka2B3 = 1./2./par[16]/par[16];
        Double_t normzrka2  = 1./2./sig_tZrKa2/par[16]*TMath::Exp(argzrka2B3);

        arg[10] = (x[0] - tzrkb1m)/sig_tZrKb1/SQRT2; // for zr kb1 tail
        Double_t argzrkb1B1 = (x[0] - tzrkb1m)/sig_tZrKb1/par[16]; 
        Double_t argzrkb1B2 = 1./par[16]/SQRT2;
        Double_t argzrkb1B3 = 1./2./par[16]/par[16];
        Double_t normzrkb1  = 1./2./sig_tZrKb1/par[16]*TMath::Exp(argzrkb1B3);

        arg[11] = (x[0] - tzrkb2m)/sig_tZrKb2/SQRT2; // for zr kb2 tail
        Double_t argzrkb2B1 = (x[0] - tzrkb2m)/sig_tZrKb2/par[16]; 
        Double_t argzrkb2B2 = 1./par[16]/SQRT2;
        Double_t argzrkb2B3 = 1./2./par[16]/par[16];
        Double_t normzrkb2  = 1./2./sig_tZrKb2/par[16]*TMath::Exp(argzrkb2B3);

        arg[12] = (x[0] - tzrkb3m)/sig_tZrKb3/SQRT2; // for zr kb3 tail
        Double_t argzrkb3B1 = (x[0] - tzrkb3m)/sig_tZrKb3/par[16]; 
        Double_t argzrkb3B2 = 1./par[16]/SQRT2;
        Double_t argzrkb3B3 = 1./2./par[16]/par[16];
        Double_t normzrkb3  = 1./2./sig_tZrKb3/par[16]*TMath::Exp(argzrkb3B3);

        arg[13] = (x[0] - par[25])/par[28]/SQRT2; // for bremsstrahlung - peak
        Double_t argbremsB1 = (x[0] - par[25])/par[28]/par[27]; 
        Double_t argbremsB2 = 1./par[27]/SQRT2;
        Double_t argbremsB3 = 1./2./par[27]/par[27];
        Double_t normbrems  = 1./2./par[28]/par[27]*TMath::Exp(argbremsB3);



/*
        arg[30] = (x[0] - tcuka1m)/sig_tCuKa1/SQRT2;
        Double_t argcuka1B1 = (x[0] - tcuka1m)/sig_tCuKa1/par[20];
        Double_t argcuka1B2 = 1./par[20]/SQRT2;
        Double_t argcuka1B3 = 1./2./par[20]/par[20];
        Double_t normcuka1  = 1./2./sig_tCuKa1/par[20]*TMath::Exp(argcuka1B3);
        arg[31] = (x[0] - tcuka2m)/sig_tCuKa2/SQRT2;
        Double_t argcuka2B1 = (x[0] - tcuka2m)/sig_tCuKa2/par[20];
        Double_t argcuka2B2 = 1./par[20]/SQRT2;
        Double_t argcuka2B3 = 1./2./par[20]/par[20];
        Double_t normcuka2  = 1./2./sig_tCuKa2/par[20]*TMath::Exp(argcuka2B3);
        arg[32] = (x[0] - tcukbm) /sig_tCuKb/SQRT2;
        Double_t argcukbB1 = (x[0] - tcukbm)/sig_tCuKb/par[20];
        Double_t argcukbB2 = 1./par[20]/SQRT2;
        Double_t argcukbB3 = 1./2./par[20]/par[20];
        Double_t normcukb  = 1./2./sig_tCuKb/par[20]*TMath::Exp(argcukbB3);

*/
        
        ticuzr_function = 

 // ----------- Titanium -----------
              par[3]  / sig_TiKa1* TMath::Exp(-0.5 * arg[0] * arg[0]) //par[4] = ti ka1 gain   Ti ka1 gauss
            + TiKa2_RI/ TiKa1_RI * par[3] / sig_TiKa2 * TMath::Exp(-0.5 * arg[1] * arg[1]) // 100/50 for ti ... ti ka2 gauss
            + par[5]  * par[3]   / sig_TiKb * TMath::Exp(-0.5 * arg[2] * arg[2]) // par[5]=ti kb to ka gain ratio ... ti kb gauss

            + par[3]  * par[18]  * normtika1 * TMath::Exp(argtika1B1) * TMath::Erfc(arg[14] + argtika1B2)// par[15]=ti tail gain ratio ka1 ... ti ka1 tail
            + par[3]  * par[18]  * TiKa2_RI / TiKa1_RI * normtika2 * TMath::Exp(argtika2B1) * TMath::Erfc(arg[15] + argtika2B2) // ti ka2 tail
            + par[3]  * par[18]  * par[5] * normtikb * TMath::Exp(argtikbB1) * TMath::Erfc(arg[16] + argtikbB2) // kb tail

// -----------  Zirconium ------------

            + par[11] / sig_ZrKa1* TMath::Exp(-0.5 * arg[3] * arg[3]) // par[11] = zr ka1 gain   Zr ka 1 gauss
            + ZrKa2_RI/ ZrKa1_RI * par[11] / sig_ZrKa2 * TMath::Exp(-0.5 * arg[4] * arg[4])
            + par[13] * par[11]  / sig_ZrKb1 * TMath::Exp(-0.5 * arg[5] * arg[5]) // par[13]=zr kb1 to ka1 gain ratio (0.7 in one fit for lngs)
	    + ZrKb2_RI / ZrKb1_RI * par[13] * par[11] / sig_ZrKb2 * TMath::Exp(-0.5 * arg[6] * arg[6]) // zr kb2 gauss
	    + ZrKb3_RI / ZrKb1_RI * par[13] * par[11] / sig_ZrKb3 * TMath::Exp(-0.5 * arg[7] * arg[7]) // zr kb3 gauss

            + par[11] * par[15]  * normzrka1 * TMath::Exp(argzrka1B1) * TMath::Erfc(arg[8] + argzrka1B2) // par[15] = tail gain zr ka ... tail for zr ka1
            + ZrKa2_RI/ ZrKa1_RI * par[11] * par[15]  * normzrka2 * TMath::Exp(argzrka2B1) * TMath::Erfc(arg[9] + argzrka2B2) // tail for zr ka2
            + par[13] * par[11] * par[17]  * normzrkb1 * TMath::Exp(argzrkb1B1) * TMath::Erfc(arg[10] + argzrkb1B2) // tail for zr kb1
            + ZrKb2_RI / ZrKb1_RI * par[13] * par[11] * par[17]  * normzrkb2 * TMath::Exp(argzrkb2B1) * TMath::Erfc(arg[11] + argzrkb2B2) //  tail for zr kb2
            + ZrKb3_RI / ZrKb1_RI * par[13] * par[11] * par[17]  * normzrkb3 * TMath::Exp(argzrkb3B1) * TMath::Erfc(arg[12] + argzrkb3B2) // tail for zr kb3


	    + par[11] * par[24] * (1./2) * TMath::Erfc( (x[0]-par[10])/(SQRT2*sig_ZrKa1) ) // zr ka1 shelf 

// --------- Copper ----------------

            + par[7] / sig_CuKa1* TMath::Exp(-0.5 * arg[24] * arg[24]) // par[8] = cu ka1 gain.. cu ka1 gauss
            + CuKa2_RI/ CuKa1_RI * par[7] / sig_CuKa2 * TMath::Exp(-0.5 * arg[25] * arg[25]) // cu ka2
            + par[7] * par[9]  / sig_CuKb * TMath::Exp(-0.5 * arg[26] * arg[26]) // par[12]=mn kb to ka gain ratio (0.7 in one fit for lngs) - cu kb

 	    + par[7] * par[19]  * normcuka1 * TMath::Exp(argcuka1B1) * TMath::Erfc(arg[17] + argcuka1B2) // cu ka1 tail
            + par[7] * par[19]  * CuKa2_RI / CuKa1_RI * normcuka2 * TMath::Exp(argcuka2B1) * TMath::Erfc(arg[18] + argcuka2B2) // cu ka2 tail
            + par[7] * par[19]  * par[9] * normcukb * TMath::Exp(argcukbB1) * TMath::Erfc(arg[19] + argcukbB2) // cu kb tail

// ---------  Bremsstrahlung --------

	    + par[26]  * normbrems * TMath::Exp(argbremsB1) * TMath::Erfc(arg[13] + argbremsB2); // bremsstrahlungs peak
   
    }
    return ticuzr_function;
}

/* Two point energy channel calibration line */
Double_t peakLine(Double_t *x, Double_t *par){
    Double_t calib_line = par[0] + par[1] * x[0];
    return calib_line;
}


Double_t TiMnCuFullFitFunc_Energy(Double_t *x, Double_t *par)
{   
    Double_t timncu_function = 0.;
    Double_t arg[33] = {0.};    
    
//    par0 = tika1sig
//    par1 = tikbsig
//    par2 = mnka1sig
//    par3 = mnkbsig
//    par4 = background constant
//    
//    par5 = tika1
//    par6 = tika1 gain
//    par7 = tikb
//    par8 = tika to ti kb ratio
//    par[9] = mnka1
//    par10 = mnka1 gain
//    par11 = mnkb
//    par12 = mn ka to kb ratio
//    par13 = escape gain
//    par14 = tail mean shift mn and ti
//    par15 = ti ka tail gain ratio
//    par16 = ti kb tail gain ratio
//    par17 = mn ka tail gain ratio
//    par18 = mn kb tail gain ratio
//    par[19]=tail beta slope ti     
//    par[20]=tail beta slope mn and cu
//    par21 = shelf gain
//    par 22 = cuka1
//    par23 = cuka1 gain       
//    par24 = cukb
//    par25 = ratio cu ka to kb
//    par26 = tail gain cuka
//    par27 = sigcuka
//    par[28] = sigcukb
            

    //Double_t slope         = (par[9] - par[5])/(MnKa1 - TiKa1); // par[9] = init_mean_mn; par[5]=init_mean_ti; slope in ch / eV  

    Double_t mnka2m        = par[9] - (MnKa1 - MnKa2); // mn ka2 start value channel
    Double_t mnka1e        = par[9] - SiKa; // channel with a mnka1 photon absorbed by Si -> lost energy not recorded
    Double_t mnka2e        = mnka2m - SiKa; // same for mn-ka2
    Double_t mnkbe         = par[11]- SiKa; // par[11]=init_mean_mnkb, same for mn-kb
   // Double_t pmnka1m       = par[9] + par[15] * slope; // par[15]=pileup mean shift in eV (60 eV), peak shift due to pulse pile up - at our rate prob = 0 
   // Double_t pmnka2m       = mnka2m + par[15] * slope; // same for mn-ka2
   // Double_t pmnkbm        = par[11]+ par[15] * slope; // same for mn-kb
    Double_t tmnka1m       = par[9] - par[14]; // par[14]=Tail mean shift [eV]
    Double_t tmnka2m       = mnka2m - par[14]; // same for mn-ka2
    Double_t tmnkbm        = par[11]- par[14]; // same for mn-kb

    Double_t sig_MnKa1     = par[2]; // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_MnKa2     = par[2]; 
    Double_t sig_MnKb      = par[3];
    Double_t sig_MnKa1e    = par[2] - 20 ;
    Double_t sig_MnKa2e    = par[2] - 20 ;
    Double_t sig_MnKbe     = par[3] - 20;
  //  Double_t sig_pMnKa1    = par[18] * sqrt( slope * pmnka1m* SiW * par[3] + par[4] * par[4] ); // par[18] is sigma broadening of pileup
  //  Double_t sig_pMnKa2    = par[18] * sqrt( slope * pmnka2m* SiW * par[3] + par[4] * par[4] );
  //  Double_t sig_pMnKb     = par[18] * sqrt( slope * pmnkbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tMnKa1    = par[2];
    Double_t sig_tMnKa2    = par[2];
    Double_t sig_tMnKb     = par[3];

    Double_t tika2m        = par[5] - (TiKa1 - TiKa2);           // TiKa2 mean channel    
    Double_t tika1e        = par[5] - SiKa;
    Double_t tika2e        = tika2m - SiKa;
    Double_t tikbe         = par[7] - SiKa; // par[7]=Ti-Kb channel
  //  Double_t ptika1m       = par[5] + par[15] * slope * 0.8;             // added 0.8 factor 23/02/2010
  //  Double_t ptika2m       = tika2m + par[15] * slope * 0.8;
  //  Double_t ptikbm        = par[7] + par[15] * slope * 0.8;
    Double_t ttika1m       = par[5] - par[14];
    Double_t ttika2m       = tika2m - par[14];
    Double_t ttikbm        = par[7] - par[14];

    Double_t sig_TiKa1     = par[0];
    Double_t sig_TiKa2     = par[0];
    Double_t sig_TiKb      = par[1];
    Double_t sig_TiKa1e    = par[0] - 20;
    Double_t sig_TiKa2e    = par[0] - 20;
    Double_t sig_TiKbe     = par[1] - 20;
 //   Double_t sig_pTiKa1    = par[18] * sqrt( slope * ptika1m* SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKa2    = par[18] * sqrt( slope * ptika2m* SiW * par[3] + par[4] * par[4] );
 //   Double_t sig_pTiKb     = par[18] * sqrt( slope * ptikbm * SiW * par[3] + par[4] * par[4] );
    Double_t sig_tTiKa1    = par[0];
    Double_t sig_tTiKa2    = par[0];
    Double_t sig_tTiKb     = par[1];

    Double_t cuka2m        = par[22] - (CuKa1 - CuKa2); // cu ka2 start value channel
    Double_t cuka1e        = par[22] - SiKa; // channel with a cuka1 photon absorbed by Si -> lost energy not recorded
    Double_t cuka2e        = cuka2m - SiKa; // same for cu-ka2
    Double_t cukbe         = par[24]- SiKa; // par[24]=init_mean_cukb, same for cu-kb
    Double_t tcuka1m       = par[22] - par[14]; // par[14]=Tail mean shift [eV]
    Double_t tcuka2m       = mnka2m - par[14]; // same for cu-ka2
    Double_t tcukbm        = par[24]- par[14]; // same for cu-kb

    Double_t sig_CuKa1     = par[27]; // SiW = ionization energy of Si @ 77 Kelvin = 3.81 eV; par[3]=fano factor; par[4]= constant noise; 8.15 eV should be ionization energy of Si; the band gap is 1.16 eV?? -> formula is true bc: var(Ne)=F*Ne; sigma(Ch)=slope*Eion*sigma(Ne); var=sigma^2; 
    Double_t sig_CuKa2     = par[27]; 
    Double_t sig_CuKb      = par[28];
    Double_t sig_CuKa1e    = par[27] - 20;
    Double_t sig_CuKa2e    = par[27] - 20;
    Double_t sig_CuKbe     = par[28] - 20;

    Double_t sig_tCuKa1    = par[27];
    Double_t sig_tCuKa2    = par[27];
    Double_t sig_tCuKb     = par[28];

    if( sig_TiKa1 * sig_TiKa2 * sig_TiKb * sig_MnKa1 * sig_MnKa2 * sig_MnKb 
            * sig_CuKa1 * sig_CuKa2 * sig_CuKb 
            * sig_tTiKa1 * sig_tTiKa2 * sig_tTiKb * sig_tMnKa1 * sig_tMnKa2 * sig_tMnKb  
            * sig_MnKa1e * sig_MnKa2e * sig_MnKbe * sig_TiKa1e * sig_TiKa2e * sig_TiKbe 
            != 0){
       
        arg[0]  = (x[0] - par[5])/sig_TiKa1;
        arg[1]  = (x[0] - tika2m)/sig_TiKa2;
        arg[2]  = (x[0] - par[7])/sig_TiKb;
        arg[21] = (x[0] - tika1e)/sig_TiKa1e;   // Escape peak
        arg[22] = (x[0] - tika2e)/sig_TiKa2e;
        arg[23] = (x[0] - tikbe)/sig_TiKbe;
        
        arg[3]  = (x[0] - par[9])/sig_MnKa1;
        arg[4]  = (x[0] - mnka2m)/sig_MnKa2;
        arg[5]  = (x[0] - par[11])/sig_MnKb;
        arg[6]  = (x[0] - mnka1e)/sig_MnKa1e;  // Escape peak
        arg[7]  = (x[0] - mnka2e)/sig_MnKa2e;
        arg[20] = (x[0] - mnkbe) /sig_MnKbe;

        arg[24]  = (x[0] - par[22])/sig_CuKa1; // cu ka1
        arg[25]  = (x[0] - cuka2m)/sig_CuKa2;
        arg[26]  = (x[0] - par[24])/sig_CuKb;
        arg[27]  = (x[0] - cuka1e)/sig_CuKa1e;  // Escape peak
        arg[28]  = (x[0] - cuka2e)/sig_CuKa2e;
        arg[29] = (x[0] - cukbe) /sig_CuKbe;

        arg[11] = (x[0] - ttika1m)/sig_tTiKa1/SQRT2;
        Double_t argtika1B1 = (x[0] - ttika1m)/sig_tTiKa1/par[19]; //par[19]=tail beta slope
        Double_t argtika1B2 = 1./par[19]/SQRT2;
        Double_t argtika1B3 = 1./2./par[19]/par[19];
        Double_t normtika1  = 1./2./sig_tTiKa1/par[19]*TMath::Exp(argtika1B3);
        arg[12] = (x[0] - ttika2m)/sig_tTiKa2/SQRT2;
        Double_t argtika2B1 = (x[0] - ttika2m)/sig_tTiKa2/par[19];
        Double_t argtika2B2 = 1./par[19]/SQRT2;
        Double_t argtika2B3 = 1./2./par[19]/par[19];
        Double_t normtika2  = 1./2./sig_tTiKa2/par[19]*TMath::Exp(argtika2B3);
        arg[13] = (x[0] - ttikbm) /sig_tTiKb/SQRT2;
        Double_t argtikbB1 = (x[0] - ttikbm)/sig_tTiKb/par[19];
        Double_t argtikbB2 = 1./par[19]/SQRT2;
        Double_t argtikbB3 = 1./2./par[19]/par[19];
        Double_t normtikb  = 1./2./sig_tTiKb/par[19]*TMath::Exp(argtikbB3);
/*
        arg[14] = (x[0] - pmnka1m)/sig_pMnKa1;
        arg[15] = (x[0] - pmnka2m)/sig_pMnKa2;
        arg[16] = (x[0] - pmnkbm) /sig_pMnKb;
*/
        arg[17] = (x[0] - tmnka1m)/sig_tMnKa1/SQRT2;
        Double_t argmnka1B1 = (x[0] - tmnka1m)/sig_tMnKa1/par[20];
        Double_t argmnka1B2 = 1./par[20]/SQRT2;
        Double_t argmnka1B3 = 1./2./par[20]/par[20];
        Double_t normmnka1  = 1./2./sig_tMnKa1/par[20]*TMath::Exp(argmnka1B3);
        arg[18] = (x[0] - tmnka2m)/sig_tMnKa2/SQRT2;
        Double_t argmnka2B1 = (x[0] - tmnka2m)/sig_tMnKa2/par[20];
        Double_t argmnka2B2 = 1./par[20]/SQRT2;
        Double_t argmnka2B3 = 1./2./par[20]/par[20];
        Double_t normmnka2  = 1./2./sig_tMnKa2/par[20]*TMath::Exp(argmnka2B3);
        arg[19] = (x[0] - tmnkbm) /sig_tMnKb/SQRT2;
        Double_t argmnkbB1 = (x[0] - tmnkbm)/sig_tMnKb/par[20];
        Double_t argmnkbB2 = 1./par[20]/SQRT2;
        Double_t argmnkbB3 = 1./2./par[20]/par[20];
        Double_t normmnkb  = 1./2./sig_tMnKb/par[20]*TMath::Exp(argmnkbB3);


        arg[30] = (x[0] - tcuka1m)/sig_tCuKa1/SQRT2;
        Double_t argcuka1B1 = (x[0] - tcuka1m)/sig_tCuKa1/par[20];
        Double_t argcuka1B2 = 1./par[20]/SQRT2;
        Double_t argcuka1B3 = 1./2./par[20]/par[20];
        Double_t normcuka1  = 1./2./sig_tCuKa1/par[20]*TMath::Exp(argcuka1B3);
        arg[31] = (x[0] - tcuka2m)/sig_tCuKa2/SQRT2;
        Double_t argcuka2B1 = (x[0] - tcuka2m)/sig_tCuKa2/par[20];
        Double_t argcuka2B2 = 1./par[20]/SQRT2;
        Double_t argcuka2B3 = 1./2./par[20]/par[20];
        Double_t normcuka2  = 1./2./sig_tCuKa2/par[20]*TMath::Exp(argcuka2B3);
        arg[32] = (x[0] - tcukbm) /sig_tCuKb/SQRT2;
        Double_t argcukbB1 = (x[0] - tcukbm)/sig_tCuKb/par[20];
        Double_t argcukbB2 = 1./par[20]/SQRT2;
        Double_t argcukbB3 = 1./2./par[20]/par[20];
        Double_t normcukb  = 1./2./sig_tCuKb/par[20]*TMath::Exp(argcukbB3);


        
        timncu_function = //backFunc(x,par)
            par[4]
            + par[6]  / sig_TiKa1* TMath::Exp(-0.5 * arg[0] * arg[0]) //par[6] = ti ka1 gain   Ti ka1 gauss
            + TiKa2_RI/ TiKa1_RI * par[6] / sig_TiKa2 * TMath::Exp(-0.5 * arg[1] * arg[1]) // 100/50 for ti
            + par[8]  * par[6]   / sig_TiKb * TMath::Exp(-0.5 * arg[2] * arg[2]) // par[8]=ti kb to ka gain ratio
            + par[6]  * par[15]  * normtika1 * TMath::Exp(argtika1B1) * TMath::Erfc(arg[11] + argtika1B2)// par[15]=ti tail gain ratio ka1
            + par[6]  * par[15]  * TiKa2_RI / TiKa1_RI * normtika2 * TMath::Exp(argtika2B1) * TMath::Erfc(arg[12] + argtika2B2) // ka2
            + par[6]  * par[16]  * par[8] * normtikb * TMath::Exp(argtikbB1) * TMath::Erfc(arg[13] + argtikbB2) // kb

            + par[10] / sig_MnKa1* TMath::Exp(-0.5 * arg[3] * arg[3]) // par[10] = mn ka1 gain
            + MnKa2_RI/ MnKa1_RI * par[10] / sig_MnKa2 * TMath::Exp(-0.5 * arg[4] * arg[4])
            + par[12] * par[10]  / sig_MnKb * TMath::Exp(-0.5 * arg[5] * arg[5]) // par[12]=mn kb to ka gain ratio (0.7 in one fit for lngs)
            + par[10] * par[17]  * normmnka1 * TMath::Exp(argmnka1B1) * TMath::Erfc(arg[17] + argmnka1B2)
            + par[10] * par[17]  * MnKa2_RI / MnKa1_RI * normmnka2 * TMath::Exp(argmnka2B1) * TMath::Erfc(arg[18] + argmnka2B2)
            + par[10] * par[18]  * par[12] * normmnkb * TMath::Exp(argmnkbB1) * TMath::Erfc(arg[19] + argmnkbB2) 

            + par[23] / sig_CuKa1* TMath::Exp(-0.5 * arg[24] * arg[24]) // cu ka1 peak
            + CuKa2_RI/ CuKa1_RI * par[23] / sig_CuKa2 * TMath::Exp(-0.5 * arg[25] * arg[25]) // cu ka2
            + par[25] * par[23]  / sig_CuKb * TMath::Exp(-0.5 * arg[26] * arg[26]) // par[12]=mn kb to ka gain ratio (0.7 in one fit for lngs) - cu kb
            + par[23] * par[26]  * normcuka1 * TMath::Exp(argcuka1B1) * TMath::Erfc(arg[30] + argcuka1B2) // cu ka1 tail
            + par[23] * par[26]  * CuKa2_RI / CuKa1_RI * normcuka2 * TMath::Exp(argcuka2B1) * TMath::Erfc(arg[31] + argcuka2B2) // cu ka2 tail
            + par[23] * par[25]  * par[26] * normcukb * TMath::Exp(argcukbB1) * TMath::Erfc(arg[32] + argcukbB2) // cu kb tail

            + par[13] * par[10]  / sig_MnKa1e * TMath::Exp(-0.5 * arg[6] * arg[6]) // Escape peak ... mn ka1
            + par[13] * par[10]  / sig_MnKa2e * TMath::Exp(-0.5 * arg[7] * arg[7]) // mn ka 2
            + par[13] * par[10]  / sig_MnKbe  * TMath::Exp(-0.5 * arg[20]* arg[20]) // mn kb


            + par[13] * par[6]  / sig_TiKa1e * TMath::Exp(-0.5 * arg[21] * arg[21]) // escape peak ti ka1
            + par[13] * par[6]  / sig_TiKa2e * TMath::Exp(-0.5 * arg[22] * arg[22]) // ka2
            + par[13] * par[6]  / sig_TiKbe  * TMath::Exp(-0.5 * arg[23] * arg[23]) // kb

            + par[21] * par[10] * (1./2) * TMath::Erfc( (x[0]-par[9])/(SQRT2*sig_MnKa1) ) // mn ka1 shelf
            + par[21] * MnKa2_RI/ MnKa1_RI * par[10] * (1./2) * TMath::Erfc( (x[0]-mnka2m)/(SQRT2*sig_MnKa2) ) // ka2
            + par[21] * par[12] * par[10] * (1./2) * TMath::Erfc( (x[0]-par[11])/(SQRT2*sig_MnKb) ) // kb
 
            + par[21] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[5])/(SQRT2*sig_TiKa1) ) // ti ka1 shelf
            + par[21] * TiKa2_RI/ TiKa1_RI * par[6] * (1./2) * TMath::Erfc( (x[0]-tika2m)/(SQRT2*sig_TiKa2) )
            + par[21] * par[8] * par[6] * (1./2) * TMath::Erfc( (x[0]-par[7])/(SQRT2*sig_TiKb) );

    }
    return timncu_function;
}