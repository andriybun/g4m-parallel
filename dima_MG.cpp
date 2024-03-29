//******************************************************************************
// DIMA model
//******************************************************************************
#ifndef dima_h_
#define dima_h_
//#include <cmath>

double forVal(double *, double &, double &, double &);
double agrVal(double *);

// Value of Forestry during multiple rotation
//double forVal(double * compr, double &rotInt, double &priceTimber, double &forValNC)
double forVal(double * compr, double &rotInt, double &priceTimber, double &forValNC)
 {
// decompressing data  
// from separate files
  double priceIndex    = compr[4];
  double r             = compr[6];
// from the combined file
  double NPP           = compr[8];
  double forest        = compr[11];                                             // has to be calculated and cumulated every year
// coefficients
  double FCUptake      = compr[16];
  double ftimber       = compr[17];
  double harvloos      = compr[18];
  double priceC        = compr[19];
  double fraclongprod  = compr[20];
  double decRateL      = compr[21];
  double decRateS      = compr[22];
  double maxRotInter   = compr[26];
  double minRotInter   = compr[27];
  double baseline      = compr[28];
  double priceTimberMax0=compr[29];
  double priceTimberMin0=compr[30];
  double plantingCosts0= compr[31];
  double sPopDens      = compr[32];
  double priceIndex0   = compr[33];
// calculations
  if(forest < 0.) {forest = 0.;}
  if(forest > 1.) {forest = 1.;}
// Mean anual carbon uptake (t-C/ha/year)
  double cUptake = 10. * NPP * FCUptake; //kg/m2 -> t/ha
// Harvestable wood-volume increment (m3/ha/year)
  double vIncr = cUptake * ftimber;
// Rotation interval of a Forest in Years
  double rotInter = 100;                                                        // vykorystovujetsya dali
  { double harvestVolume = 600. - abs(vIncr -6.) * 50.;
  if(cUptake > 0.) {rotInter = harvestVolume/vIncr;}
  if(rotInter < minRotInter) {rotInter = minRotInter;}
  if(rotInter > maxRotInter) {rotInter = maxRotInter;} }
//Costs to plant 1 ha of forest
  double plantingCosts;
  { double plantrate = (vIncr-3.)/6.;
    if(plantrate > 1.) {plantrate = 1.;}
    if(plantrate < 0.) {plantrate = 0.;}
    plantingCosts = plantrate * plantingCosts0 * priceIndex / priceIndex0; }
//Timber price
  { double sfor = (1. - forest) * 9. + 1.;
    double c4 = (priceTimberMax0 - priceTimberMin0)/99.;
    double c3 = priceTimberMin0 - c4;
//************************************************
    priceTimber = (c3 + c4 * sPopDens * sfor) * priceIndex/priceIndex0; }
//Price to harvest the timber
  double priceHarvest = priceTimber * .0;
//Harvest volume of the timber
  double woodHarvestVol = vIncr * rotInter * (1. - harvloos);
//Fraction of carbon costs during harvest
  double beta = 1. - decRateL/(decRateL+r)*fraclongprod
	     - decRateS/(decRateS+r) * (1.-fraclongprod);
//Carbon benefit
  double cBenefit = priceC * cUptake * (1. - baseline) * (((1. - pow(1.+r,-rotInter) ) /r) -
	     rotInter * (1.-beta) * pow(1.+r, -rotInter));
// Value of Forestry during one rotation
  double forestValueOne = (-plantingCosts + (priceTimber - priceHarvest)
         *woodHarvestVol) + cBenefit;
// Return values
  rotInt = rotInter;
// Value of Forestry multiple rotation No Carbon Price
  double forestValueOneNC = -plantingCosts + (priceTimber
                          - priceHarvest)*woodHarvestVol;//*pow(1+r, -rotInter);
  forValNC = forestValueOneNC / (1. - pow(1. + r, -rotInter));
  
//cout<<"forVal= "<<  forestValueOne/(1.-pow(1.+r, -rotInter))<<endl;
//cout<<"forValNC= "<<  forValueNC<<endl;
  return(forestValueOne/(1.-pow(1.+r, -rotInter)));
 }

//Net present Value of Agriculture
double agrVal(double * compr)
{
// decompressing
// from separate files
  double priceIndex    = compr[4];
// from the combined file
  double sagrsuit       = compr[46];
// coefficients
  double priceLandMinR = compr[14];
  double priceLandMaxR = compr[15];
  double sPopDens      = compr[32];
  double priceIndex0   = compr[33];
// calculations
  double priceLevel = priceLandMinR * priceIndex / priceIndex0;
//Importance of Population density
  double popImp = (log(priceLandMaxR) - log(priceLandMinR))/(2. * log(10.));  
//Importance of the Suitable for Agriculture
  double agrImp = popImp;
//cout<<"sPopDens= "<<sPopDens<<endl;
//cout<<"sagrsuit= "<<sagrsuit<<endl;
//cout<<"priceLandMinR= "<<priceLandMinR<<endl;
//cout<<"priceLandMaxR= "<<priceLandMaxR<<endl;
//cout<<"popImp= "<<popImp<<endl;
//cout<<"agrVal= "<<priceLevel * pow(sPopDens,popImp) * pow(sagrsuit,agrImp)<<endl;

  return(priceLevel * pow(sPopDens,popImp) * pow(sagrsuit,agrImp));
}

#endif
