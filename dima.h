#ifndef DIMA_H_
#define DIMA_H_

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;
//////////
// DIMA //
//////////
class dima {
public:
  double forVal();      //Value of Forestry during multiple rotation
  double forValNC();      //Value of Forestry during multiple rotation
  double agrVal();      //Net present Value of Agriculture
  double amenVal();      //Net present Value of Agriculture
  double minPriceC();      //Get the minimum Carbonprice where Forest=Argic
  int setYear(int);
  double setForest(double);
  double forValExt();   //MG: Value of Forestry using External Timber price
  double forValNCExt();      //MG: Value of Forestry during multiple rotation with 0 C price and using External wood price
  double agrVal2000();   //MG: Net present Value of Agriculture in the year 2000
  dima(int ayear
       , interpol anpp          //npp [kg-C/m2/year]
       , interpol asPopDens   //Standardised (1-10) population density
       , interpol asAgrSuit   //Standardised (1-10) agricultural suitability
       , interpol apriceIndex  //Priceindex
       , g4m::ipol<double,double> apriceIndex0 //Priceindex ref Country
       , interpol ar          //discount rate
       , g4m::ipol<double,double> apriceC      //carbon price [$/tC]
       , g4m::ipol<double,double> aplantingCosts0 //Planting costs ref Country [$/ha]
       , g4m::ipol<double,double> apriceLandMin0   //Minimum Landprice in ref country [$/ha]
       , g4m::ipol<double,double> apriceLandMax0  //Maximum Landprice in ref country [$/ha]
       , g4m::ipol<double,double> amaxRotInter   //Maximal rotation intervall [years]
       , g4m::ipol<double,double> aminRotInter       //Minimal rotation intervall [years]
       , g4m::ipol<double,double> adecLongProd  //Decay rate for long lived prod
       , g4m::ipol<double,double> adecShortProd  //Decay rate for short lived pro
       //Fraction of carbon stored in longterm products
       , interpol afracLongProd
       //Fraction of carbon substracted due to baseline considerations
       , g4m::ipol<double,double> abaseline
       //Comercial timbervolume per ton of carbon [m3/tC]
//       , g4m::ipol<double,double> aftimber
       , interpol aftimber
//         , int aftimber
       //Maximal timber price in reference country [$/m3]
       , g4m::ipol<double,double> apriceTimberMax0
       //Minimal timber price in reference country {$/m3]
       , g4m::ipol<double,double> apriceTimberMin0
       , g4m::ipol<double,double> afcuptake  //Factor of carbon uptake from npp
       , interpol agdp       //Gross domestic production
       , g4m::ipol<double,double> aharvLoos  //Harvesting losses
       , double aforest      //Share of landarea used by forest [0-1]
       , interpol awpricecorr //MG: Added for wood price correction
       , double awpricecorr0 //MG: Added for wood price correction in year 2000 when price of carbon=0
       , double arotInterM  //MG: Added rotation interval estimated from Georg's Forest Management Tool
       , double aharvWood   //MG: Added harvestable wood from Georg's Forest Management Tool
       ) : year(ayear)
	   , npp(anpp)
	   , sPopDens(asPopDens)
	   , sAgrSuit(asAgrSuit)
	   , priceIndex(apriceIndex)
	   , priceIndex0(apriceIndex0)
	   , r(ar)
	   , priceC(apriceC)
	   , plantingCosts0(aplantingCosts0)
	   , priceLandMin0(apriceLandMin0)
	   , priceLandMax0(apriceLandMax0)
	   , maxRotInter(amaxRotInter)
	   , minRotInter(aminRotInter)
	   , decLongProd(adecLongProd)
	   , decShortProd(adecShortProd)
	   , fracLongProd(afracLongProd)
	   , baseline(abaseline)
	   , ftimber(aftimber)
	   , priceTimberMax0(apriceTimberMax0)
	   , priceTimberMin0(apriceTimberMin0)
	   , fcuptake(afcuptake)
	   , gdp(agdp)
	   , harvLoos(aharvLoos)
	   , forest(aforest)
	   , wpricecorr(awpricecorr)
	   , wpricecorr0(awpricecorr0)
	   , rotInterM(arotInterM)
	   , harvWood(aharvWood)
  {
    ;
  }
  double priceTimber();      //Timber price
  double priceTimberExt();  //MG: Timber price External
  //private:
  int year;
  double rotInter();         //Rotation interval of a Forest in Years
  double woodHarvestVol();   //Harvest volume of the timber
  double cUptake();          //mean anual carbon uptake
  double vIncr();            //Harvestable wood-volum increment (m3/ha/year)
private:
  double beta();             //Fraction of carbon costs during harvest
  double cBenefit();         //Carbon benefit
  double priceHarvest();     //Price to harvest the timber
  double forestValueOne();   //Value of Forestry during one rotation
  double forestValueOneExt();//MG: Value of Forestry during one rotation  with External Timber price
  double forestValueOneNC();   //Value of Forestry during one rotation
  double forestValueOneNCExt();   //MG: Value of Forestry during one rotation with External Timber price
  double plantingCosts();    //Costs to plant 1 ha of forest

  g4m::ipol<double,double> plantingCosts0;     //Costs to Plant 1ha of forests in ref country
  interpol r;                  //Discount rate
  interpol sAgrSuit;           //Standardised Agricultural Suitability
  interpol sPopDens;           //Standardised Population Density (1-10)
  interpol npp;                //neto primary production
  interpol priceIndex;         //Price Index
  g4m::ipol<double,double> priceIndex0;        //Price Index of reference country
  g4m::ipol<double,double> priceC;             //Carbon Price
  g4m::ipol<double,double> priceLandMin0;      //Minimum Landprice in ref country
  g4m::ipol<double,double> priceLandMax0;      //Maximum Landprice in ref country
  g4m::ipol<double,double> maxRotInter;        //Maximal rotation intervall in years
  g4m::ipol<double,double> minRotInter;        //Minimal rotation intervall
  g4m::ipol<double,double> decLongProd;        //Decay rate for long lived products
  g4m::ipol<double,double> decShortProd;       //Decay rate for short lived products
  interpol fracLongProd;       //Fraction of carbon stored in longterm products
  //Fraction of carbon substracted due to baseline considerations
  g4m::ipol<double,double> baseline;

//  g4m::ipol<double,double> ftimber;          //Comercial timbervolume per ton of carbon (m3/tC)
    interpol ftimber;
//      int ftimber;

  g4m::ipol<double,double> priceTimberMax0;    //Maximal timber price in reference country
  g4m::ipol<double,double> priceTimberMin0;    //Minimal timber price in reference country
  g4m::ipol<double,double> fcuptake;            //Factor of carbon uptake from npp
  interpol gdp;                 //Gross domestic production
  g4m::ipol<double,double> harvLoos;
  double forest;                //Share of landarea used by forest [0-1]
  interpol wpricecorr;          //MG: Added for wood price correction
  double wpricecorr0;          //MG: Added for wood price correction in year 2000
  double rotInterM;            //MG: Added rotation interval estimated from Georg's Forest Management Tool
  double harvWood;             //MG: Added harvestable wood estimated from Georg's Forest Management Tool
};

double dima::setForest(double x) {
  forest = x;
  if(forest < 0.) {forest = 0.;}
  if(forest > 1.) {forest = 1.;}
  return(forest);
}

int dima::setYear(int i) {
  year = i;
  return(year);
}

double dima::cUptake() {     //mean anual carbon uptake (t-C/ha/year)
  //This value should depend on the NPP which is reduced by a factor
  //This factor should not be the same all over the world
  //Tropical takes not so much NPP because of:
  //   *)high rotation of leaves and litle brances
  //   *)Insects eat the leaves

  //  return(1.02*.13*(npp+.1)/10.);
//cout << "npp.v(year)= "<<npp.v(year)<<endl;
  return(npp.v(year) * 10. * fcuptake.v(year)); //kg/m2 -> t/ha

}

double dima::vIncr() {  //Harvestable wood-volume increment (m3/ha/year)
//cout << "ftimber.v(year)= "<<ftimber.v(year)<<endl;
  return(cUptake() * ftimber.v(year));
//  return(cUptake() * ftimber);

}

// MG: Use Georg's Optimal Rotation Time
double dima::rotInter() {          //Rotation interval of a Forest in Years

//  double t = 100;
//  double harvestVolume = 600. - abs(vIncr() -6.) * 50.;
//  if(cUptake() > 0.) {t = harvestVolume/vIncr();}
//  if(t < minRotInter.v(year)) {t = minRotInter.v(year);}
//  if(t > maxRotInter.v(year)) {t = maxRotInter.v(year);}

double t=rotInterM;  // MG: Use Georg's Optimal Rotation Time
//cout<<"t= "<<t<<endl;
return(t);

}

double dima::beta() {              //Fraction of carbon costs during harvest
  //Depends on fraction of short and long term products
  return(1. - decLongProd.v(year)/(decLongProd.v(year)+r.v(year))*fracLongProd.v(year)
	 - decShortProd.v(year)/(decShortProd.v(year)+r.v(year))
	 *(1.-fracLongProd.v(year)));
}

double dima::cBenefit() {       //Carbon benefit (Eq. 3)
//	cout << "priceC.v(year)"<<"\t"<<priceC.v(year)<<endl;
  return(priceC.v(year) * cUptake() * (1. - baseline.v(year)) *
	 ( ((1. - pow(1.+r.v(year),-rotInter()) ) /r.v(year)) -
	   rotInter() * (1.-beta()) * pow(1.+r.v(year), -rotInter())));
}

double dima::woodHarvestVol() {  //Harvest volume of the timber during 1 rotation period
//  return(vIncr() * rotInter() * (1. - harvLoos.v(year)));
//cout << harvWood * rotInter() <<"   "<<vIncr() * rotInter() * (1. - harvLoos.v(year))<<endl;
return(harvWood * rotInter());  // MG: Georg's harvestable wood

  //Disturbances can also be mentioned
}

double dima::priceTimber() {     //Timber price internal
  //double c4 = (priceTimberMax0.v(year) - priceTimberMin0.v(year))/9.;
  //double c3 = priceTimberMin0.v(year) - c4;
  //return((c3 + c4 * sPopDens.v(year))
  //	 * priceIndex.v(year)/priceIndex0.v(year));

  double sfor = (1. - forest) * 9. + 1.;
  double c4 = (priceTimberMax0.v(year) - priceTimberMin0.v(year))/99.; 
  double c3 = priceTimberMin0.v(year) - c4;
  return((c3 + c4 * sPopDens.v(year) * sfor)
	 * priceIndex.v(year)/priceIndex0.v(year));

}

double dima::priceTimberExt() {     //MG: Timber price external
  //double c4 = (priceTimberMax0.v(year) - priceTimberMin0.v(year))/9.;
  //double c3 = priceTimberMin0.v(year) - c4;
  //return((c3 + c4 * sPopDens.v(year))
  //	 * priceIndex.v(year)/priceIndex0.v(year));

  double sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  double c4 = (priceTimberMax0.v(2000) - priceTimberMin0.v(2000))/99.;
  double c3 = priceTimberMin0.v(2000) - c4;
  return((c3 + c4 * sPopDens.v(2000) * sfor)
	 * priceIndex.v(2000)/priceIndex0.v(2000)
	 * wpricecorr.v(year)/wpricecorr0);
}

  //double shareFuelwood = 1./(1 + exp(-1.09525 + 0.30313*gdp.v(year)
  //			     + 1.50520*forest));
  //double harvestAmountYear = woodHarvestVol()/rotInter();
  //double harvFuel = harvestAmountYear * shareFuelwood;
  //if(harvFuel > 10.) {harvFuel = 10.;}
  //double harvInd = harvestAmountYear - harvFuel;
  //if(harvInd > 4.) {harvInd = 4.;}
  //double indWoodPrice = priceTimberMin0.v(year) + gdp.v(year)
  //  * exp(3.7619 - 1.098*forest - 0.5928 * harvInd);
  //double fuelWoodPrice = priceTimberMin0.v(year) + gdp.v(year)
  //  * exp(1.38493 - 0.97638*forest + 0.05744*harvFuel);
  //double price = fuelWoodPrice * shareFuelwood
  //  + indWoodPrice * (1. - shareFuelwood);
  //return(price);


double dima::priceHarvest() {     //Price to harvest the timber
  return(priceTimber() * .0);
//Beside harvesting costs also thinning costs, branc-removal,... can be
//considered
}

double dima::plantingCosts() {    //Costs to plant 1 ha of forest
  //return(plantingCosts0.v(year)*priceIndex.v(year)/priceIndex0.v(year));
  //Maybe these costs do not ocure on the second rotation intervall
  //because of *)natural regeneration *)coppice forests
  double plantrate = (vIncr()-3.)/6.;
  if(plantrate > 1.) {plantrate = 1.;}
  if(plantrate < 0.) {plantrate = 0.;}
  return(plantrate*plantingCosts0.v(year)*priceIndex.v(year)/
	 priceIndex0.v(year));
}


double dima::forestValueOneNC() {//Value of Forestry one rotation NoCarbonPrice
  //return(-plantingCosts() +
  // (priceTimber()
  // -priceHarvest())*woodHarvestVol()*pow(1+r.v(year), -rotInter()));
  return(-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimber()
    -priceHarvest())*woodHarvestVol());  // MG: I deleted *pow(1+r.v(year), -rotInter()) to make it similar to forestValueOne
}

double dima::forestValueOneNCExt() {//MG: Value of Forestry one rotation NoCarbonPrice using External wood price
  //return(-plantingCosts() +
  // (priceTimber()
  // -priceHarvest())*woodHarvestVol()*pow(1+r.v(year), -rotInter()));
  return(-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimberExt()
    -priceHarvest())*woodHarvestVol());  // MG: I deleted *pow(1+r.v(year), -rotInter()) to make it similar to forestValueOne
}

double dima::forValNC() { //Value of Forestry multiple rotation No Carbon Price
  return(forestValueOneNC()/(1.-pow(1.+r.v(year), -rotInter())));
//cout<<"rotInter()= "<<rotInter()<<endl;
}

double dima::forValNCExt() { //MG: Value of Forestry multiple rotation No Carbon Price using External wood price
  return(forestValueOneNCExt()/(1.-pow(1.+r.v(year), -rotInter())));
}

double dima::forestValueOne() { //Value of Forestry during one rotation (Eq.1)
  //return(-plantingCosts() +
  // (priceTimber()
  //  -priceHarvest())*woodHarvestVol()*pow(1+r.v(year), -rotInter())
  // +cBenefit());
  return(-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimber()
    	  -priceHarvest())*woodHarvestVol()
  	 +cBenefit()
);
}

double dima::forestValueOneExt() { //MG: Value of Forestry during one rotation External
  //return(-plantingCosts() +
  // (priceTimber()
  //  -priceHarvest())*woodHarvestVol()*pow(1+r.v(year), -rotInter())
  // +cBenefit());
  return(-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimberExt()
    	  -priceHarvest())*woodHarvestVol()
  	 +cBenefit()
);


}

double dima::forVal() { //Value of Forestry during multiple rotation (Eq.4)
  return(forestValueOne()/(1.-pow(1.+r.v(year), -rotInter())));
//  return(forestValueOne()/(1.-pow(1.+r.v(year), -rotInter()))+cBenefit());  
}

double dima::forValExt() { //MG: Value of Forestry during multiple rotation with External Timber price

//cout<<"forValExt()= "<<forestValueOneExt()/(1.-pow(1.+r.v(year), -rotInter()))<<endl;
  return(forestValueOneExt()/(1.-pow(1.+r.v(year), -rotInter())));
//  return(forestValueOne()/(1.-pow(1.+r.v(year), -rotInter()))+cBenefit());  
}

double dima::agrVal() {   //Net present Value of Agriculture (Eq.5)
  double priceLevel = priceLandMin0.v(year)
    * priceIndex.v(year)/priceIndex0.v(year);
  //Importance of Population density
  double popImp = (log(priceLandMax0.v(year))
		   - log(priceLandMin0.v(year)))/(2. * log(10.));
  //Importance of the Suitable for Agriculture
  double agrImp = popImp;
//cout << "priceLandMin0.v(year)= "<<priceLandMin0.v(year)<<"\t priceIndex.v(year)= "<<priceIndex.v(year)<<"\t priceIndex0.v(year)= "<<priceIndex0.v(year);
//cout << "\t priceLevel= "<<priceLevel<<"\t sPopDens.v(year)= "<<sPopDens.v(year)<<"\t popImp= "<<popImp<<"\t sAgrSuit.v(year)= "<<sAgrSuit.v(year); 
//cout << "\t agrVal()= "<<priceLevel * pow(sPopDens.v(year),popImp) * pow(sAgrSuit.v(year),agrImp)<<endl;
  return(priceLevel * pow(sPopDens.v(year),popImp)
	 * pow(sAgrSuit.v(year),agrImp));

}

////MG: Attantion: agrVal changed! Only 2000 value is estimated here
double dima::agrVal2000() {   //MG: Net present Value of Agriculture in the year 2000
  double priceLevel = priceLandMin0.v(2000)
    * priceIndex.v(2000)/priceIndex0.v(2000);
  //Importance of Population density
  double popImp = (log(priceLandMax0.v(2000))
		   - log(priceLandMin0.v(2000)))/(2. * log(10.));
  //Importance of the Suitable for Agriculture
  double agrImp = popImp;
  return(priceLevel * pow(sPopDens.v(2000),popImp)
	 * pow(sAgrSuit.v(year),agrImp));

}

double dima::amenVal() {   //Value of amenity
  double priceLevel = priceLandMin0.v(year)
    * priceIndex.v(year)/priceIndex0.v(year);
  //Importance of Population density
  double popImp = (log(priceLandMax0.v(year))
		   - log(priceLandMin0.v(year)))/(2. * log(10.));
  //Importance of GDP
  double gdpImp = popImp;
  return(priceLevel * pow(sPopDens.v(year),popImp)
	 * pow(gdp.v(year),gdpImp));
}

double dima::minPriceC() {
  //Get the minimum Carbonprice where Forest=Argic (Eq.6)
  return((agrVal() * (1-pow(1+r.v(year),-rotInter())) +
	  plantingCosts() -
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimber()	  
	  
	   -priceHarvest())*woodHarvestVol()*pow(1+r.v(year), -rotInter())
	  )
	 /
	 (cUptake() * (1. - baseline.v(year)) *
	  ( ((1. - pow(1.+r.v(year),-rotInter()) ) /r.v(year)) -
	    rotInter() * (1.-beta()) * pow(1.+r.v(year), -rotInter()))));
}
////////////// END DIMA //////////////////////////////////////

#endif
