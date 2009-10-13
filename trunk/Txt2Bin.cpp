#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <ctime>

#include "misc.h"                       // ipol class for interpolation
#include "dataStruct_MG.h"                 // data structures
#include "griddata3_MG.h"                  // GridData class                       (v1.1, v1.2 or v1.3)
#include "forest_MG.h"                     // definitions
#include "readInput_MG.cpp"                // code for reading from files
#include "forestCalculations_MG.cpp"       // code for calculations
#include "countryData.h"                // data by countries

#define index(x,y) (x*25 + y)

int main ()
 {
cout << "Starting..." << endl;

  time_t start=time(NULL);
// resolution of model  
  ResLatitude = int(floor(180/GridStepLat));
  ResLongitude = int(floor(360/GridStepLon));
  
  //*******************************
//** Reading coefficients file **
//*******************************
cout << "Calling readCoeff(coeff)" << endl;
  readCoeff(coeff);
  g4m::coeffStruct coeff_bin;
//*******************************

//*******************************************
//** Reading detailed input data from file **
//*******************************************
  dataDetStruct plots ;                         // structure with data
  readInputDet(plots);                         // plots[<elNum>].<variable>[year]
  cout << "> Time is " << difftime(time(NULL),start) << " sec." << endl;
//******************************************************************************

  dataDetStruct::iterator it = plots.begin();

  int i=0;

  int xy[41969][2];
  double *alldata = new double[41969*25];
  if (alldata == 0) {
cout << "Cannot allocate memory for alldata" << endl;
  // error assigning memory for alldata. Take measures.
  };

  while (it != plots.end())
   {

   if ((it->FOREST[2000])+(it->CROP[2000])+(it->BUILTUP[2000])>1)
   {it->FOREST.insert(2000,(1-((it->CROP[2000])+(it->BUILTUP[2000]))));   }

   xy[i][0]=(it->x);
   xy[i][1]=(it->y);
//cout<<"x= "<<(it->x)<<"    y= "<<(it->y)<<endl;

   alldata[index(i,0)]=(it->COUNTRY[2000]);
   alldata[index(i,1)]=(it->POTVEG[2000]);
   alldata[index(i,2)]=(it->USED[2000]);
   alldata[index(i,3)]=(it->PROTECT[2000]);   

   alldata[index(i,4)]=(it->LANDAREA[2000])*100; // land area in ha
//cout<<"Landarea[2005] = "<<(it->LANDAREA[2005])<<endl;
   alldata[index(i,5)]=(it->NPP[2000]);
//cout<<"NPP[2005] = "<<(it->NPP[2005])<<endl;

   alldata[index(i,6)]=(it->POPDENS[2000]);
   alldata[index(i,7)]=(it->POPDENS[2005]);
//if (alldata[index(i,7)]!=0) {cout<<"POPDENS[2005]="<< alldata[index(i,7)]<<endl;}

   alldata[index(i,8)]=(it->SAGRSUIT[2000]);
   alldata[index(i,9)]=(it->AGRSUIT[2000]);

//cout<<"AGRSUIT[2000] = "<<(it->AGRSUIT[2000])<<endl;

   alldata[index(i,10)]=(it->PRICEINDEX[2000]);

   alldata[index(i,11)]=(it->BIOMASS[2000]);

   alldata[index(i,12)]=(it->FOREST[2000]);

   alldata[index(i,13)]=(it->R[2000]);

   alldata[index(i,14)]=(it->GDP[2000]);
   alldata[index(i,15)]=(it->GDP[2005]);

   alldata[index(i,16)]=(it->BUILTUP[2000]);
   alldata[index(i,17)]=(it->BUILTUP[2005]); 
//cout<<"BUILTUP[2000] = "<<(it->BUILTUP[2000])<<endl;

   alldata[index(i,18)]=(it->CROP[2000]);
   alldata[index(i,19)]=(it->CROP[2005]);

   alldata[index(i,20)]=(it->FRACLONGPROD[2000]);

   alldata[index(i,21)]=(it->CORRUPTION[2000]);

   alldata[index(i,22)]=(it->SLASHBURN[2000]);

   alldata[index(i,23)]=(it->SPOPDENS[2000]);//MG: added to use in forest_calculations (dima)
   alldata[index(i,24)]=(it->SPOPDENS[2005]);   

//cout<<"SPOPDENS[2000] = "<<(it->SPOPDENS[2000])<<endl;

i++;
it++;
   }  
// Reading coefficients into an array
double coefficients[28];  
  coefficients[0] = coeff.bYear;
  coefficients[1] = coeff.eYear;
  coefficients[2] = coeff.cellsInteract;
  coefficients[3] = coeff.inclAffor;
  coefficients[4] = coeff.noPay;
  coefficients[5] = coeff.uBiomass;
  coefficients[6] = coeff.litter;
  coefficients[7] = coeff.SOC;
  coefficients[8] = coeff.PriceLandMinR[0];
  coefficients[9] = coeff.PriceLandMaxR[0];
  coefficients[10] = coeff.FCuptake[0];
  coefficients[11] = coeff.FTimber[0];
  coefficients[12] = coeff.HarvLoos[0];
  coefficients[13] = coeff.PriceC[0];
  coefficients[14] = coeff.FracLongProd[0];
  coefficients[15] = coeff.decRateL[0];
  coefficients[16] = coeff.decRateS[0];
  coefficients[17] = coeff.SlashBurn[0];
  coefficients[18] = coeff.FreqAid[0];
  coefficients[19] = coeff.PriceCAid[0];
  coefficients[20] = coeff.MaxRotInter[0];
  coefficients[21] = coeff.MinRotInter[0];
  coefficients[22] = coeff.baseline[0];
  coefficients[23] = coeff.PriceTimberMaxR[0];
  coefficients[24] = coeff.PriceTimberMinR[0];
  coefficients[25] = coeff.PriceIndexE[0];
  coefficients[26] = coeff.PlantingCostsR[0];
  coefficients[27] = coeff.sPopDens[0];



//for (int i=0; i<28;i++){   cout << "i= "<<i<<"    Coefficient "<< coefficients[i]<<endl;}

int size_coeff = sizeof(coefficients);
int size_plots = sizeof(plots);
int size_xy = sizeof(xy);
int size_alldata = sizeof(*alldata)*41969*25;

cout << "size_coeff " << size_coeff << endl;
cout << "size_plots " << size_plots << endl;
cout << "xy " << size_xy << endl;
cout << "size_alldata " << size_alldata << endl;

//
//  double totFor1=0;
//  while (it != plots.end())
//   {
//   if ((it->FOREST[2000])+(it->CROP[2000])+(it->BUILTUP[2000])>1)
//   {it->FOREST.insert(2000,(1-((it->CROP[2000])+(it->BUILTUP[2000])))); 
//   }
//    forestCover2000.set(it->x,it->y,it->FOREST[2000]); 
//it++;
//   }  


 
 
   fstream xy_out;
  xy_out.open("C:/MGusti/CurrentWork/GFM/georgPrgs/dima/DeforAforCCurves_growth/Optimise/xy.bin", ios::out|ios::binary);
  if(xy_out.is_open()) {
    xy_out.write ((char*)&xy, size_xy);
    xy_out.close();
  } 
  else {cout << "Cannot create xy.bin!" << endl;}

  fstream coeffout;
  coeffout.open("C:/MGusti/CurrentWork/GFM/georgPrgs/dima/DeforAforCCurves_growth/Optimise/coef.bin", ios::out|ios::binary);
  if(coeffout.is_open()) {
    coeffout.write ((char*)&coefficients, size_coeff);
    coeffout.close();
  } 
  else {cout << "Cannot create coef.bin!" << endl;}


//---------------
// Writing alldata to the binary file data.bin    
  fstream dataout;
//  dataout.open("C:/MGusti/CurrentWork/GFM/georgPrgs/dima/DeforAforCCurves_growth/Optimise/data.bin", ios::out|ios::binary);
  dataout.open("C:/MGusti/CurrentWork/GFM/georgPrgs/dima/DeforAforCCurves_growth/Optimise/data_cramerNPP.bin", ios::out|ios::binary);
  if(dataout.is_open()) {
    dataout.write ((char*)alldata, size_alldata);  
    dataout.close();
  } 
  else {cout << "Cannot create data.bin!" << endl;    }
//----------------------------------------------

//  fstream coeffin;
//  coeffin.open("C:/MGusti/CurrentWork/GFM/georgPrgs/dima/DeforAforCCurves_growth/Optimise/coef.bin", ios::in|ios::binary);
//  if(coeffin.is_open()) {
//    coeffin.read ((char*)&coeff_bin, sizeof(size_coeff));
//    coeffin.close();
//
//cout << "eYear=" << coeff_bin.cellsInteract << endl;
//
//} else
//{cout << "Cannot open coef.bin!" << endl;}

delete []alldata;

  cout << "> Time is " << difftime(time(NULL),start) << " sec." << endl;
}      
           
