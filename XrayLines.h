/************************************
           XrayLines.h
           25  Mar  2011
           04  May  2012 
               Added kaonic Kapton 
               EM lines. 
               (W/O VP & recoil corrections)
           21  Sep  2012 
               Change the Kapton 
               lines into VP corrected
               values. 

 ************************************/
// "_RI" refers to Relative Intensity. 
#ifndef XRAYLINES_H
#define XRAYLINES_H


/******** KHe4 EM Lines ********/
// Koike calculation
static const Double_t KHe4La = 6463.5;
static const Double_t KHe4Lb = 8721.7;
static const Double_t KHe4Lr = 9766.8;

/******** KHe3 EM Lines ********/
// K.G. with Vacuum Polarization
static const Double_t KHe3La = 6224.6;
// K.G. only:
//static const Double_t KHe3La = 6211.02;
static const Double_t KHe3Lb = 8384.82;
static const Double_t KHe3Lr = 9390.96;

/******** Kp EM Lines *******/
// T.Ito et al PRC 58(1997)
// Calculated by Koike.
/*
static const Double_t KpKa = 6479.9;
static const Double_t KpKb = 7676.7;
static const Double_t KpKr = 8095.5;
static const Double_t KpKd = 8289.4;
static const Double_t KpKe = 8394.7;
static const Double_t KpKz = 8458.1;;
static const Double_t KpKy = 8499.3;
static const Double_t KpKinf = 8633.9;
*/
static const Double_t KpKa = 6481.2;
static const Double_t KpKb = 7678.0;
static const Double_t KpKr = 8096.8;
static const Double_t KpKd = 8290.6;
static const Double_t KpKe = 8395.9;
static const Double_t KpKz = 8459.4;;
static const Double_t KpKy = 8500.6;
static const Double_t KpKinf = 8633.9;

// Kp shift and width (Exp. results)
// KpX KEK PS E228
static const Double_t KpXEps    = -323.0;
static const Double_t KpXEpsSta = 63;
static const Double_t KpXEpsSys = 11;
static const Double_t KpXGm     = 407;
static const Double_t KpXGmSta  = 208;
static const Double_t KpXGmSys  = 100;
// DEAR 
static const Double_t DearEps    = -193;
static const Double_t DearEpsSta = 37;
static const Double_t DearEpsSys = 6;
static const Double_t DearGm     = 249;
static const Double_t DearGmSta  = 111;
static const Double_t DearGmSys  = 30;
// SIDDHARTA 
static const Double_t SidtEps    = -283;
static const Double_t SidtEpsSta = 36;
static const Double_t SidtEpsSys = 6;
static const Double_t SidtGm     = 541;
static const Double_t SidtGmSta  = 79;
static const Double_t SidtGmSys  = 22;

/******** Kd  EM Lines *******/
// K. G. with Vacuum Polarization
// using Z. Phys. D 15, 321-325 (1990), Eq.2
static const Double_t KdKa = 7834.0;
static const Double_t KdKb = 9280.2;
static const Double_t KdKr = 9786.2;

//****  K-Kapton lines *******/  04 May 2012 Added Vacuum Polarization and recoil.
// KC lines:
//static const Double_t KC54 = 10197.41;
//static const Double_t KC65 = 5539.15;
//static const Double_t KC76 = 3339.88;
//static const Double_t KC75 = 8879.03;
static const Double_t KC54 = 10216.5;
static const Double_t KC65 =  5544.9;
static const Double_t KC76 =  3341.7;
static const Double_t KC75 =  8885.8;

// KN lines:
//static const Double_t KN76 = 4573.68;
//static const Double_t KN65 = 7585.44;
//static const Double_t KN75 = 12159.13;
static const Double_t KN76 = 4577.68;
static const Double_t KN65 =  7595.7;
static const Double_t KN75 = 12171.1;

// KO lines: 
//static const Double_t KO76 = 6001.13;
//static const Double_t KO65 = 9952.93;
static const Double_t KO76 = 6006.8;
static const Double_t KO65 = 9968.7;
//*****************************/

// KAl lines:
//static const Double_t KAl76 = 16058.03;
//static const Double_t KAl87 = 10421.82;
//static const Double_t KAl98 = 7144.97;
static const Double_t KAl76 = 16088.3;
static const Double_t KAl87 = 10435.1;
static const Double_t KAl98 =  7150.4;

///////////////////////////////////



/******** Silicon *******/
// K-alpha 
static const Double_t SiKa = 1739. ;

/******** Iron **********/
// K-alpha1 
static const Double_t FeKa1 = 6403.84;
static const Double_t FeKa1_RI = 100 ;
// K-alpha2 
static const Double_t FeKa2 = 6390.84;
static const Double_t FeKa2_RI = 50  ;
// K-beta1,3
static const Double_t FeKb1_3 = 7058.0;
static const Double_t FeKb1_3_RI = 17.;

/******** Calcium **********/
// K-alpha1 
static const Double_t CaKa1 = 3691.7;
static const Double_t CaKa1_RI = 100 ;
// K-alpha2 
static const Double_t CaKa2 = 3688.1;
static const Double_t CaKa2_RI = 50  ;
// K-beta1,3
static const Double_t CaKb1_3 = 4012.7;
static const Double_t CaKb1_3_RI = 13.;

/******** Copper ********/
/* Pauli obeying transitions */
// K-alpha1
static const Double_t CuKa1 = 8047.78;
static const Double_t CuKa1_W  = 2.11;
static const Double_t CuKa1_RI = 100 ;
// K-alpha2
static const Double_t CuKa2 = 8027.83;
static const Double_t CuKa2_W  = 2.17;
static const Double_t CuKa2_RI = 51  ;
// K-beta1
static const Double_t CuKb1 = 8905.29;
static const Double_t CuKb1_RI = 17  ;

/* Pauli violating transitions */
// according to INFN-13-21/LNF, 21 November 2013
// K-alpha1
static const Double_t CuKa1_V = 7746.73;
// K-alpha2
static const Double_t CuKa2_V = 7728.83;
// K-beta1
static const Double_t CuKb1_V = 8531.69;


/******* Zinc ********/
// K-alpha1 
static const Double_t ZnKa1 = 8638.86;
static const Double_t ZnKa1_RI = 100 ;
// K-alpha2
static const Double_t ZnKa2 = 8615.78;
static const Double_t ZnKa2_RI = 51  ;
// K-beta1
static const Double_t ZnKb1 = 9572.0;
static const Double_t ZnKb1_RI = 17 ;
// L-alpha1
static const Double_t ZnLa1 = 1011.7;
// L-alpha2
static const Double_t ZnLa2 = 1011.7;

/******** Titanium ********/
// K-alpha1
static const Double_t TiKa1 = 4510.84;
static const Double_t TiKa1_W  = 1.16;
static const Double_t TiKa1_RI = 100;
// K-alpha2
static const Double_t TiKa2 = 4504.86;
static const Double_t TiKa2_W  = 1.18;
static const Double_t TiKa2_RI = 50 ;
// K-beta1
static const Double_t TiKb1 = 4931.81;
static const Double_t TiKb1_RI = 15 ;

/******* Manganese *******/
// K-alpha1 
static const Double_t MnKa1 = 5898.75;
static const Double_t MnKa1_RI = 100 ;
// K-alpha2 
static const Double_t MnKa2 = 5887.65;
static const Double_t MnKa2_RI = 50  ;
// K-beta1 
static const Double_t MnKb1 = 6490.45;
static const Double_t MnKb1_RI = 17  ;

/******* Nickel *********/ 
// K-alpha1 
static const Double_t NiKa1 = 7478.15;
static const Double_t NiKa1_RI = 100;

static const Double_t NiKa2 = 7460.89;
static const Double_t NiKa2_RI = 51;

static const Double_t NiKb1_3 = 8264.7; 
static const Double_t NiKb1_RI = 17;

/******* Aunium  **********/
// L-alpha1
static const Double_t AuLa1 = 9713.3;
static const Double_t AuLa1_RI = 100.;
// L-alpha2
static const Double_t AuLa2 = 9628.0;
static const Double_t AuLa2_RI = 11.;
// L-beta 1
static const Double_t AuLb1 = 11442.3;
static const Double_t AuLb1_RI = 67  ;
// L-beta 2
static const Double_t AuLb2 = 11584.7;
static const Double_t AuLb2_RI = 23  ;

/******* Zirconium ********/
// K-alpha1
static const Double_t ZrKa1 = 15775.1;
static const Double_t ZrKa1_RI = 100 ;
// K-alpha2 
static const Double_t ZrKa2 = 15690.9;
static const Double_t ZrKa2_RI = 52  ;
// K-beta1
static const Double_t ZrKb1 = 17667.8;
static const Double_t ZrKb1_RI = 15  ;
// K-beta2
static const Double_t ZrKb2 = 17970  ;
static const Double_t ZrKb2_RI = 3   ;
// K-beta3
static const Double_t ZrKb3 = 17654  ;
static const Double_t ZrKb3_RI = 8   ;

/****** Lead *******/
// L-alpha1
static const Double_t PbLa1 = 10551.5;
static const Double_t PbLa1_RI = 100.;
// L-alpha2
static const Double_t PbLa2 = 10449.5;
static const Double_t PbLa2_RI = 11.;

/****** Bromine *******/
// K-alpha1
static const Double_t BrKa1 = 11924.2;
static const Double_t BrKa1_RI  = 100.;
// K-alpha2
static const Double_t BrKa2 = 11877.6;
static const Double_t BrKa2_RI  = 52.;
// K-beta1
static const Double_t BrKb1 = 13291.4;
static const Double_t BrKb1_RI = 14.;
// K-beta2
static const Double_t BrKb3 = 13284.5;
static const Double_t BrKb3_RI = 7.;




/****** Silver *******/ 
// K-alpha1
static const Double_t AgKa1 = 22162.9; 
static const Double_t AgKa1_RI = 100.;
// K-alpha2
static const Double_t AgKa2 = 21990.3; 
static const Double_t AgKa2_RI = 53.;
// K-beta1
static const Double_t AgKb1 = 24942.4;
static const Double_t AgKb1_RI = 16.;
// K-beta3
static const Double_t AgKb3 = 24911.5;
static const Double_t AgKb3_RI = 8.;

/****** Antinomy *******/ // In alloy with Lead, cable sheath etc. 
// K-alpha1
static const Double_t SbKa1 = 26359.1;
static const Double_t SbKa1_RI = 100.;
// K-alpha2 
static const Double_t SbKa2 = 26110.8;
static const Double_t SbKa2_RI = 53.;

/****** Tin ******/
// K-alpha1 
static const Double_t SnKa1 = 25271.3;
static const Double_t SnKa1_RI = 100.;
// K-alpha2
static const Double_t SnKa2 = 25044.0;
static const Double_t SnKa2_RI = 53.;

/****** Pd ******/ // In the electronics and soldering material
// K-alpha1
static const Double_t PdKa1 = 21177.1;
static const Double_t PdKa1_RI = 100.;
// K-alpha2
static const Double_t PdKa2 = 21020.1;
static const Double_t PdKa2_RI = 53.;
// K-beta1 
static const Double_t PdKb1 = 23818.7;
static const Double_t PdKb1_RI = 16;
// K-beta2 
static const Double_t PdKb2 = 24299.1;
static const Double_t PdKb2_RI = 4.;
// K-beta3 
static const Double_t PdKb3 = 23791.1;
static const Double_t PdKb3_RI = 8.;

#endif
