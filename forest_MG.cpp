//   Name:      g4m
//   Author:    Andriy Bun, based on works of ...
//   Date:      9.September.2009

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

#include "misc.h"                          // ipol class for interpolation
#include "increment.cpp"
#include "ageStruct.cpp"
#include "interpol.h"
#include "dima.h"
#include "dataStruct_MG.h"                 // data structures
#include "griddata3_MG.h"                  // GridData class                       (v1.1, v1.2 or v1.3)
#include "griddata2.h"                     // data of any type on x-y grid by Andriy Bun
#include "fillContainer.cpp"

//#include "countryData.h"
//#include "countryData_MG1.h"                   // data by countries
#include "countryData.h"                   // data by countries with Poles regions and printToFile

#include "forest_MG.h"                     // definitions
#include "cPrices_NatDefra.h"
#include "MAI_country.cpp"
#include "hurdle_and_deforaffor_annex1_FM.h"
#include "readInput_MG.cpp"                // code for reading from files
//#include "forestCalculations_MG.cpp"       // code for calculations
#include "initManagedForest4.cpp"
#include "calc.cpp"
#include "woodHarvestStatCountry.cpp" // FAO stat on wood harvest by countries
#include "adjustManagedForest.cpp"  // Adjusting FM every year to match wood production (external)

//******************************************************************************
//***********************************  MAIN  ***********************************
//******************************************************************************
int main ()
 {
  time_t start=time(NULL);
// resolution of model  
  ResLatitude = int(floor(180/GridStepLat));
  ResLongitude = int(floor(360/GridStepLon));
// Setting years for output
  set<int> years;
  years.insert(2000);
  years.insert(2005);
  years.insert(2010);
  years.insert(2020);
  years.insert(2030);
  years.insert(2040);
  years.insert(2050);

// Setting regions for calculations:
  regions.insert(1);
  regions.insert(3);  
  regions.insert(5);  
  regions.insert(6);
  regions.insert(7);  
  regions.insert(8);  
  regions.insert(9);
  regions.insert(10);  
  regions.insert(12);  
  regions.insert(18);
  regions.insert(27);  

//  regions.insert(2);
//  regions.insert(4);  
//  regions.insert(11);  
//  regions.insert(13);  
//  regions.insert(14);  
//  regions.insert(15);  
//  regions.insert(16);  
//  regions.insert(17);  
//  regions.insert(19);  
//  regions.insert(20);  
//  regions.insert(21);  
//  regions.insert(22);  
//  regions.insert(23);  
//  regions.insert(24);  
//  regions.insert(25);  
//  regions.insert(26);  
   
//*******************************
//** Reading coefficients file **
//*******************************
  readCoeff(coeff);
// starting and ending years of calculations
  byear = coeff.bYear;
  eyear = coeff.eYear;
//*****************************************************
//** Reading input files (with resolution 0.5 x 0.5) **
//*****************************************************
//*******************************************
//** Reading detailed input data from file **
//*******************************************
  dataDetStruct plots;                         // structure with data
  int numRecords = readInputDet(plots);        // plots[<elNum>].<variable>[year]
  cout << "> Time is " << difftime(time(NULL),start) << " sec." << endl;
//******************************************************************************
//***************************** start calculations *****************************
//******************************************************************************
//** Initializing forest cover array by gridcells **
//**************************************************
  griddata harvestGrid = griddata(ResLongitude,ResLatitude,0);
  griddata maiForest = griddata(ResLongitude,ResLatitude,0);
  griddata rotationForest = griddata(ResLongitude,ResLatitude,0);
  griddata rotationForestNew = griddata(ResLongitude,ResLatitude,0);

  griddata2<char> decisionGrid = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> managedForest = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> thinningForest = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> rotationType = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> thinningForestNew = griddata2<char>(ResLongitude,ResLatitude,0);
  
  // Setup forest increment table
  g4m::incrementTab fi(-4.25, -1.39, 0.329, 169., 81., 4.48, -1.39, 0.780,
                       -0.176, 1./4.5, 0.5, -1.08, -7.7956, -0.2362, -0.6316,
                        24.62079, 0.46573, -0.02426, 2.58884, -0.30364, -0.5782,
		                0.8513, 0.0289, 1.4720, 0.5, 2., 1.6, 10., 0.25, 600., 1.);

  MAI_country();
  cPrices();
  hurdle_aff_deff();
  woodHarvestStatCountry();
  for (int i=0; i < NumberOfCountries; i++) {  
    EmissionsCurCountry[i]=0.; 
    EmissionsCurAfforCountry[i]=0.;
  }
  ageStructVector cohort_all;
  ageStructVector newCohort_all;
  cohort_all.reserve(numRecords);
  newCohort_all.reserve(numRecords);
  datGlobal dat_all;
  initManagedForest(plots, fi, dat_all, maiForest, thinningForest, rotationType, 
                    managedForest, rotationForest, harvestGrid);
//*************************
//**** loop by prices *****
//*************************
  for (int i=0;i<=0;i++) {         // MG: PriceC loop
cout << "Price C value " << i << endl;    
    initLoop(i, plots, fi, cohort_all, newCohort_all, dat_all, maiForest, thinningForest, rotationForest);
//************************
//**** loop by years *****
//************************
    int year = byear;
    do {
      cout << "Processing year " << year << endl;
//      decision.setYear(year);
      int Age = year-byear;
      if (year > byear)
       {
cout << "Adjusting FM .."<< endl;
       adjustManagedForest(plots, fi, cohort_all, 
              newCohort_all, dat_all, maiForest, 
              thinningForest, rotationForest, managedForest,
              rotationForestNew, thinningForestNew,
              rotationType, harvestGrid, year);
cout << "finished Adjusting FM .."<< endl;
       }


//************************************
//** processing data from all plots **
//************************************
      dataDetStruct::iterator it = plots.begin();
      while (it != plots.end()) {
        if (it->PROTECT[2000]==0) {
      if (regions.find(it->POLESREG[2000]) != regions.end()) { // Test only some regions
            int asID = it->asID;
//if (year > byear) cout << "start Calc"<< endl;
//if ((cohort_all[asID].getArea(0) == 0.)||(cohort_all[asID].getArea(0) >= 1./500.))
{
            calc(*it, fi, *cohort_all[asID], *newCohort_all[asID], dat_all[asID], managedForest,
                 maiForest, rotationForest, rotationForestNew, thinningForest, thinningForestNew, 
                 harvestGrid,year, i, asID);
}                 
//if (year > byear) cout<< "finished Calc"<<endl;                 
          }                        // End if Country ...
        }                          // End if not protected
        it++;
      }                            // End loop by plots
    year++;
    } while (year <= eyear);       // End loop by years
     CountriesNforCover.PrintToFile("newForestHa",coeff.bYear,coeff.eYear,1);
     CountriesAfforHaYear.PrintToFile("afforHaYear",coeff.bYear,coeff.eYear,1);          
     CountriesNforTotC.PrintToFile("afforSinkC",coeff.bYear,coeff.eYear,1); 
     CountriesAfforCYear.PrintToFile("afforSinkCYear",coeff.bYear,coeff.eYear,1);  
  
     CountriesOforCover.PrintToFile("oldForestHa",coeff.bYear,coeff.eYear,1); 
     CountriesOforBiomassC.PrintToFile("oldForestBiomassC",coeff.bYear,coeff.eYear,1); 
     CountriesOforTotC.PrintToFile("oldForestEmissC",coeff.bYear,coeff.eYear,1); 
     CountriesDeforHaYear.PrintToFile("deforHaYear",coeff.bYear,coeff.eYear,1); 
     CountriesDeforCYear.PrintToFile("deforCYear",coeff.bYear,coeff.eYear,1);    

     CountriesWoodHarvestM3Year.PrintToFile("harvest_m3_year_Iloop",coeff.bYear,coeff.eYear,1); 
     CountriesWoodLoosCYear.PrintToFile("harvest_lostCYear_Iloop",coeff.bYear,coeff.eYear,1);
     
     CountriesManagedForHa.PrintToFile("managedForHa",coeff.bYear,coeff.eYear,1);
     CountriesManagedCount.PrintToFile("managedCount",coeff.bYear,coeff.eYear,1);
  }                                // End loop by prices
  cout << "> Working time is " << difftime(time(NULL),start) << " sec." << endl;
  system("pause");
 }
//******************************************************************************
// end main
//******************************************************************************
