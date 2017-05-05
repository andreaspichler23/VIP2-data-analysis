/*************************************
            CalibFunction.h
              May 16 2014
              H. Shi
              Modified for VIP2 use
            
**************************************/
#ifndef CalibFunction_H
#define CalibFunction_H

extern Double_t backFunc( Double_t *x, Double_t *par);
extern Double_t shelfFunc(Double_t *x, Double_t *par);

extern Double_t caka1Func(Double_t *x, Double_t *par);
extern Double_t tika1Func(Double_t *x, Double_t *par);
extern Double_t mnka1Func(Double_t *x, Double_t *par);
extern Double_t cuka1Func(Double_t *x, Double_t *par);
extern Double_t zrka1Func(Double_t *x, Double_t *par);

//extern Double_t peakFunc(Double_t *x,  Double_t *par);      
//extern Double_t kaFunc(Double_t *x,  Double_t *par);
extern Double_t peakLine(Double_t *x, Double_t *par);

extern Double_t catailFunc(Double_t *x, Double_t *par);
extern Double_t titailFunc(Double_t *x, Double_t *par);
extern Double_t mntailFunc(Double_t *x, Double_t *par);
extern Double_t cutailFunc(Double_t *x, Double_t *par);
extern Double_t zrtailFunc(Double_t *x, Double_t *par);

extern Double_t capileFunc(Double_t *x, Double_t *par);
extern Double_t tipileFunc(Double_t *x, Double_t *par);
extern Double_t mnpileFunc(Double_t *x, Double_t *par);
extern Double_t cupileFunc(Double_t *x, Double_t *par);
extern Double_t zrpileFunc(Double_t *x, Double_t *par);

extern Double_t tiescFunc(Double_t  *x, Double_t *par);
extern Double_t mnescFunc(Double_t  *x, Double_t *par);
extern Double_t cuescFunc(Double_t  *x, Double_t *par);


extern Double_t TiCuSourceFitFunc(Double_t *x, Double_t *par);
extern Double_t TiCuFullFitFunc(Double_t *x, Double_t *par);
extern Double_t TiMnFullFitFunc(Double_t *x, Double_t *par);
extern Double_t TiMnCuFullFitFunc(Double_t *x, Double_t *par);
extern Double_t TiCuZrFullFitFunc(Double_t *x, Double_t *par);

#endif
