// remove int asID

void calc(g4m::dataStruct &it, g4m::incrementTab &fi, g4m::ageStruct &cohort, g4m::ageStruct &newCohort,
          dat &singleCell, griddata2<char> &managedForest, griddata &maiForest, griddata &rotationForest,
          griddata &rotationForestNew, griddata2<char> &thinningForest, griddata2<char> &thinningForestNew, 
          griddata &harvestGrid, int year, int i, int asID)
// griddata &harvestGrid
 {
//  decision.setYear(year);
//MG: setup forest Age
//******************************************************************************
if (year >= 1991) {
  cout << asID << endl;
  if (asID >= 63) {
    cout << "1" << endl;
    system("pause");
  }
}
//******************************************************************************
  int Age = year-byear;
  int xi = (it.x);
  int yi = (it.y);
  int Country = (int)it.COUNTRY[2000];
// --------------Applying FM policy -----------------------------------
//******************************************************************************
if (year >= 1991) {
  if (asID >= 63) {
    cout << "2 before get: " << xi << "\t" << yi << endl;
    system("pause");
  }
}
//******************************************************************************
  int rotationTimeCurr = rotationForest.get(xi,yi);
//******************************************************************************
if (year >= 1991) {
  if (asID >= 63) {
    cout << "3 after get" << endl;
    system("pause");
  }
}
//******************************************************************************
//  cohort.setRotPeriod(rotationTimeCurr);       
//  cohort.setStockingdegree(thinningForest.get(xi,yi));
//  newCohort.setRotPeriod(rotationForestNew.get(xi,yi));       
//  newCohort.setStockingdegree(thinningForestNew.get(xi,yi));
  

  if (managedForest.get(xi,yi) > 0) {
// Changing rotation time to the optimal one
//    if ((year == 2010) && (rotationTimeCurr < rotMaxBmTh)) {
//      rotationTimeCurr = rotMaxBmTh; //MG: change rotation time to the scenario specific
//      cohort.setRotPeriod(rotationTimeCurr);
//    }
  } 
//  else 
//  {
//    if (year == 1990) rotationTimeNew = optRotUnmanaged+1;
//    if (year == 1990) {
//      rotationTimeNew = optRot+1;
//      cohort.setRotPeriod(rotationTimeNew);
//    }
//  }
int fl = 0;
//if (asID > 14) fl = 1;
//  g4m::ageStruct::v res = cohort.aging1(fl); //MG
//if (asID > 13) {cohort.prn(); system("pause");}
//******************************************************************************
if (year >= 1991) {
  if (asID >= 63) {
    cout << "4 B4 aging" << endl;
    system("pause");
  }
}
//******************************************************************************
  g4m::ageStruct::v res = cohort.aging(); //MG
//******************************************************************************
if (year >= 1991) {
  if (asID >= 63) {
    cout << "5 after aging" << endl;
    system("pause");
  }
}
//******************************************************************************
  g4m::ageStruct::v newRes = newCohort.aging();
//******************************************************************************
//******************************************************************************
if (year >= 1991) {
  if (asID >= 63) {
    cout << "6 after aging" << endl;
    system("pause");
  }
}
//******************************************************************************
if (year >= 1991) {
  if (asID >= 63) {
    cout << "2" << endl;
    system("pause");
  }
}
//******************************************************************************
  map<string, interpol> data;
  fillContainer(it,data);     
  double sawnW = 0.;
  double restW = 0.;
  double sawnThW = 0.;
  double restThW = 0.;
  double sawnWlost = 0.;
  double restWlost = 0.;
  double sawnThWlost = 0.;
  double restThWlost = 0.;
  double sawnWnew = 0.;
  double restWnew = 0.;
  double sawnThWnew = 0.;
  double restThWnew = 0.;
  double sawnWlostNew = 0.;
  double restWlostNew = 0.;
  double sawnThWlostNew = 0.;
  double restThWlostNew = 0.;
  double forNPV_CH = 0.;
  double abBiomassO = 0.;
  if (managedForest.get(xi,yi)>0) {  
    sawnW = res.enSw ;           // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
    restW = res.enRw ;           // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
    sawnThW = res.vnSw ;         // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
    restThW = res.vnRw ;         // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
    sawnWnew = newRes.enSw;      // MG: get harvestable sawnwood for the set (new) forest tC/ha for final cut. New (planted) forest
    restWnew = newRes.enRw;      // MG: get harvestable restwood for the set (new) forest tC/ha for final cut. 
    sawnThWnew = newRes.vnSw;    // MG: get harvestable sawnwood for the set (new) forest tC/ha for thinning. 
    restThWnew = newRes.vnRw;    // MG: get harvestable restwood for the set (new) forest tC/ha for thinning.
  } else {
    sawnWlost = res.enSw ;           // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
    restWlost = res.enRw ;           // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
    sawnThWlost = res.vnSw ;         // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
    restThWlost = res.vnRw ;         // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
    sawnWlostNew = newRes.enSw;      // MG: get harvestable sawnwood for the set (new) forest tC/ha for final cut. 
    restWlostNew = newRes.enRw;      // MG: get harvestable restwood for the set (new) forest tC/ha for final cut. 
    sawnThWlostNew = newRes.vnSw;    // MG: get harvestable sawnwood for the set (new) forest tC/ha for thinning. 
    restThWlostNew = newRes.vnRw;    // MG: get harvestable restwood for the set (new) forest tC/ha for thinning.
  }
  abBiomassO = cohort.getBm();
  double harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
  double harvWoodLost = (sawnWlost + restWlost + sawnThWlost + restThWlost); // Total current "lost" wood in the cell, tC (in remote forests)
  double harvWoodNew = (sawnWnew + restWnew + sawnThWnew + restThWnew) * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
  double harvWoodLostNew = (sawnWlostNew + restWlostNew + sawnThWlostNew + restThWlostNew); // Total current "lost" wood in the cell, tC (in remote forests)
// Rotation time fitted to get certain biomass under certain MAI (w/o thinning)
  int biomasRot = 0;
  int rotMAI = 0;
  int PriceCi = PriceCiS[i];
  if (singleCell.forestShare > 0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi)> 0) {
//    biomasRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);
    rotMAI = fi.gTopt(maiForest.get(xi,yi), 1);
  }
//  if (biomasRot < 0) biomasRot = 0;
//  if ((biomasRot == 0)||(harvWood < 0)) harvWood = 0.;
  if ((singleCell.Rotation <= 0)||(harvWood < 0)) harvWood = 0.;
// MG: setup price of carbon
  coeff.PriceC.clear();
  coeff.PriceC.insert(0, PriceCi * LinPrice2020[Age] * data["CORRUPTION"].v());
  string regprice, regprice0;
  char prstr[5];
  char regstr[3];
  int2str((int)data["POLESREG"].v(),regstr);
  if (PriceCi!=0) {
    int2str(PriceCi,prstr);
    regprice = "re"+string(regstr)+"price"+ string(prstr);
  } else {
    regprice = "re"+string(regstr)+"price0";
  }
  regprice0 = "re"+string(regstr)+"price0";
  {
// ---- Checking NPV of forestry for current harvest------------       
    dima decision(year
         , data["NPP"]
         , data["SPOPDENS"]
         , data["SAGRSUIT"]
         , data["PRICEINDEX"]
		 , coeff.PriceIndexE
		 , data["R"]
		 , coeff.PriceC
		 , coeff.PlantingCostsR
		 , coeff.PriceLandMinR
		 , coeff.PriceLandMaxR
		 , coeff.MaxRotInter
		 , coeff.MinRotInter
		 , coeff.decRateL
		 , coeff.decRateS
		 , data["FRACLONGPROD"]
		 , coeff.baseline
		 , data["FTIMBER"]
		 , coeff.PriceTimberMaxR
		 , coeff.PriceTimberMinR
		 , coeff.FCuptake
         , data["GDP"]
         , coeff.HarvLoos
         , data["FOREST"].v(year)
         , wprice[regprice]
         , wprice[regprice0].v(2000)
         , rotationTimeCurr
         , harvWood + harvWoodLost * data["FTIMBER"].v() );
    if (year<=2005) {
      forNPV_CH = decision.forVal(); // MG: use internal NPV
    } else {
      forNPV_CH = decision.forValExt();
    } // MG: use External NPV
//-------------------------------------------------------------------------------
  }
//  dima decision(year
//		       , data["NPP"]
//		       , data["SPOPDENS"]
//		       , data["SAGRSUIT"]
//		       , data["PRICEINDEX"]
//		       , coeff.PriceIndexE
//		       , data["R"]
//		       , coeff.PriceC
//		       , coeff.PlantingCostsR
//		       , coeff.PriceLandMinR
//		       , coeff.PriceLandMaxR
//		       , coeff.MaxRotInter
//		       , coeff.MinRotInter
//		       , coeff.decRateL
//		       , coeff.decRateS
//		       , data["FRACLONGPROD"]
//		       , coeff.baseline
//		       , data["FTIMBER"].v()
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v(year)
//		       , wprice[regprice]
//		       , wprice[regprice0].v(2000)
//               , singleCell.Rotation
//               , singleCell.potHarvest);
//             

double harvMAI = maiForest.get(xi,yi)*data["FTIMBER"].v()*(1-coeff.HarvLoos.v());

  dima decision(year
		       , data["NPP"]
		       , data["SPOPDENS"]
		       , data["SAGRSUIT"]
		       , data["PRICEINDEX"]
		       , coeff.PriceIndexE
		       , data["R"]
		       , coeff.PriceC
		       , coeff.PlantingCostsR
		       , coeff.PriceLandMinR
		       , coeff.PriceLandMaxR
		       , coeff.MaxRotInter
		       , coeff.MinRotInter
		       , coeff.decRateL
		       , coeff.decRateS
		       , data["FRACLONGPROD"]
		       , coeff.baseline
		       , data["FTIMBER"]
		       , coeff.PriceTimberMaxR
		       , coeff.PriceTimberMinR
		       , coeff.FCuptake
		       , data["GDP"]
		       , coeff.HarvLoos
		       , data["FOREST"].v(year)
		       , wprice[regprice]
		       , wprice[regprice0].v(2000)
               , rotMAI
               , harvMAI);


  double EmissionsCur=0.;
  double EmissionsProductCur= 0.;  
  double EmissionsLitterCur = 0.;  
  double EmissionsSOCCur = 0.;      
  double EmissionsSlashBurnCur = 0.;
  double EmissionsDeadBurnCur = 0.;
  double EmissionsCRootBurnCur = 0.;
  double EmissionsSOCAfforCur = 0.; 
  double EmissionsLitterAfforCur = 0.;
  double EmissionsAfforCur = 0.;
  double EmissionsCurNP=0.;
  double EmissionsProductCurNP= 0.;  
  double EmissionsLitterCurNP = 0.;  
  double EmissionsSOCCurNP = 0.;      
  double EmissionsSlashBurnCurNP = 0.;
  double EmissionsDeadBurnCurNP = 0.;
  double EmissionsCRootBurnCurNP = 0.;
  double EmissionsSOCAfforCurNP = 0.; 
  double EmissionsLitterAfforCurNP = 0.;
  double EmissionsAfforCurNP = 0.; 
  double PlantPhytHaBmGr = 0.;
  double PlantPhytHaBmGrNP =0.;
  double PlantPhytHaBmGrG =0.;
  double PlantPhytHaBlGr = 0.;
  double PlantPhytHaBlGrNP = 0.;
//Timber price
  double TimberPrice = 0.;  
  if (year<=2005) {
    TimberPrice = decision.priceTimber();
  } else { // MG: use internal TimberPrice
    TimberPrice = decision.priceTimberExt(); // MG: use External TimberPrice
  }
//Forestry Value
  double fval = 0.;  
  if (year<=2005) {
    fval = decision.forVal(); // MG: use internal NPV
  } else {
    fval = decision.forValExt();
  } // MG: use External NPV
//Forestry Value No carbon price
  double fvalNC = 0.;  
  if (year<=2005) {
    fvalNC = decision.forValNC(); // MG: use internal NPV
  } else {
    fvalNC = decision.forValNCExt(); // MG: use External NPV
  }
//Agricultural Value
  double aval = 0.;
  if (year<=2005) {
    aval = decision.agrVal(); // MG: use internal G4M land price
  } else {
    aval = decision.agrVal2000() * lprice[regprice].v(year)/lprice[regprice0].v(2000); // MG: use GLOBIOM price
  }
  double crop = data["CROP"].v(year);
  double builtup = data["BUILTUP"].v(year);
  double freeLand = 1 - (singleCell.forestShare + crop + builtup);
  double spopdens = data["SPOPDENS"].v(year);
//defrorestation speed
  double defShare = 0.;
//MG: Afforestation speed
  double affShare = 0.;
  double OforestShare = singleCell.OforestShare;
  double AforestShare = singleCell.AforestShare;
  double maxfor = 1. - (data["BUILTUP"].v(year)
	    		  + data["CROP"].v(year));
  if (maxfor < 0.) maxfor = 0;
    double gdp = 1644.; // MG: gdp definition
    if (data["POPDENS"].v(year) > 0.) {
      gdp = 200000. * data["GDP"].v(year) / data["POPDENS"].v(year);
      gdp = gdp * deflator; // Global GDP deflator GDP(1995)/GDP(2000)=8.807/10 (World Bank)
    }
//MG: Deforestation
    if(data["AGRSUIT"].v(year) > 0 && OforestShare > 0) {
      defShare = 0.05/(1. + exp(1.799e+00
				+ 2.200e-01/OforestShare
				+ 1.663e-01/data["AGRSUIT"].v(year)
				- 4.029e-02 * data["POPDENS"].v(year)
				+ 5.305e-04 * data["POPDENS"].v(year) * data["POPDENS"].v(year)
				+ 1.282e-04 * gdp ));
//cout<<"defshare= "<<defShare<<"\t deforRate_opt= "<<deforRate_opt[Country-1]<<"\t Country= "<<Country<<endl;
      defShare *= deforRate_opt[Country-1] ;
    }
//MG: Afforestation
    
    if ((OforestShare + AforestShare < maxfor) && (data["POTVEG"].v()<9) ) { // MG: We afforest only places, where potential vegetation is forest
      affShare = 0.01/(1+exp(1.+0.1/data["AGRSUIT"].v(year)+1/(0.001*gdp)));
      if (affShare < 1./(data["LANDAREA"].v()*1000000.)) affShare = 0.;   // minimun one tree (1m^2)
//cout<<"affshare= "<<affShare<<"\t afforRate_opt= "<<afforRate_opt[Country-1]<<endl;
      affShare *= afforRate_opt[Country-1]; 
    }
//cout<<"defshare= "<<defShare<<"\t affShare= "<<affShare<<endl;

// MG: End of gdp definition
  double defIncome = 0.;
//Pay if Carbon get's to air (Harvest goes to products)
  double pDefIncome = abBiomassO *
	 				(TimberPrice * data["FTIMBER"].v()
					* (1. -coeff.HarvLoos.v())
					- coeff.PriceC.v(year) * (1 + data["R"].v(year))
					* (data["FRACLONGPROD"].v(year) * coeff.decRateL.v(year)
					/ (coeff.decRateL.v(year) + data["R"].v(year))
					+ data["FRACLONGPROD"].v(year) * coeff.decRateS.v(year) //MG: mistake: must be DECRATES but not DECRATEL
					/ (coeff.decRateS.v(year) + data["R"].v(year))));
//*******************************************************************************************************
/*if (asID >= 14) {
cout << data["R"].v(year) << endl;
cout << data["CLITTERHA"].v() << endl;
cout << data["DECWOOD"].v(year) << endl;
cout << data["DECHERB"].v(year) << endl;
cout << data["CBELOWHA"].v() << endl;
cout << data["SOCHA"].v() << endl;
cout << data["DECSOC"].v(year) << endl;
system("pause");
}*/
//*******************************************************************************************************
  double pDefLoss = - coeff.PriceC.v(year) * (1 + data["R"].v(year))
                  * (data["CLITTERHA"].v()*(0.3*data["DECWOOD"].v(year)/(data["DECWOOD"].v(year)+data["R"].v(year))
                  + 0.7*data["DECHERB"].v(year)/(data["DECHERB"].v(year)+data["R"].v(year))) 
                  + data["CBELOWHA"].v()*0.3*data["DECHERB"].v(year)/(data["DECHERB"].v(year)+data["R"].v(year))
                  + data["SOCHA"].v()*data["DECSOC"].v(year)/(data["DECSOC"].v(year)+data["R"].v(year)));
//Immediate Pay if deforested (Slash and Burn)
  double sDefIncome = abBiomassO *
					(TimberPrice * data["FTIMBER"].v()
					* (1. -coeff.HarvLoos.v())
					- coeff.PriceC.v(year));
  double sDefLoss = - coeff.PriceC.v(year) * (data["CDEADHA"].v()+ data["CBELOWHA"].v(year)*0.7);
  defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
			+ sDefIncome * data["SLASHBURN"].v()
			+ pDefLoss + sDefLoss;
  if ((aval + defIncome) > ((fval ) * Hurdle_opt[Country-1] ) && (OforestShare > 0)) {  // MG: adjust the multiplier to account for a forest saving policy in some countries
    OforestShare -= defShare; // Decrease Forest share
    if (OforestShare > maxfor) OforestShare = maxfor;
    if (OforestShare < 0.) OforestShare = 0.;
  }
//cout<<"NPP= "<< data["NPP"].v()<<"\t aval= "<<aval<<"\t fval= "<<fval<<"\t maxfor= "<<maxfor<<"\t potveg= "<<data["POTVEG"].v()<<"\t OforestShare= "<<OforestShare <<endl;
  if (data["NPP"].v()>0) {
    if ((aval < (fval + coeff.PriceC.v()* (((data["SOCHA"].v()*0.4/(8*decision.rotInter()))
       + 5/(1.053*decision.rotInter()))/data["R"].v())) * Hurdle_opt[Country-1])
       && (OforestShare + AforestShare < maxfor) && data["POTVEG"].v()<9) { // MG: We afforest only places, where potential vegetation is forest
      AforestShare += affShare;
//cout<< "AforestShare0= "<<AforestShare<<endl;
      if (AforestShare > maxfor) AforestShare = maxfor;
      if (AforestShare < 0.) {
        AforestShare = 0.;
        affShare = 0.;
      }	    
      if (OforestShare + AforestShare > maxfor) AforestShare = maxfor - OforestShare;
      if (AforestShare < 0.) AforestShare = 0.;
//cout<< "AforestShare2= "<<AforestShare<<endl;

    }
  }
  if (OforestShare > maxfor) {OforestShare = maxfor;}
  if (OforestShare < 0.) {OforestShare = 0.;}	       
  if (AforestShare < 0.) {AforestShare = 0.;}	       
  if (AforestShare > maxfor) {AforestShare = maxfor;}	       
  if (OforestShare + AforestShare > maxfor) {AforestShare = maxfor - OforestShare;}
  if (AforestShare < 0.) {AforestShare = 0.;}	
//cout<< "AforestShare3= "<<AforestShare<<endl;  	    
  singleCell.forestAgeShare[Age] = AforestShare - singleCell.AforestSharePrev;
  if (singleCell.forestAgeShare[Age] < 0.) {singleCell.forestAgeShare[Age]=0.;}
//cout<<"singleCell.AforestSharePrev= "<<singleCell.AforestSharePrev<<"\t singleCell.forestAgeShare[Age]= "<<singleCell.forestAgeShare[Age]<<endl;
  newCohort.afforest(singleCell.forestAgeShare[Age]); // MG: Afforest
  singleCell.deforestA[Age] = singleCell.prevOForShare - OforestShare;
  if (singleCell.deforestA[Age] < 0.) {singleCell.deforestA[Age]=0.;}
 
    if (year < 2000) { // MG: We start deforestation in 2000, because we have initial data for 2000 (discussion with Hannes Bottcher 10.09.09)
    OforestShare = singleCell.OforestShare;
    AforestShare = 0.;  
    }
  res = cohort.deforest(singleCell.deforestA[Age],0); //MG: deforest set forest and get deforested biomass. Devide by refForShare to get per ha value
  double defBiomass = (res.enSw + res.enRw + res.vnSw + res.vnRw);///deforestA[Age]; // deforested biomass per ha
  singleCell.OforestShare = OforestShare;
  singleCell.AforestShare = AforestShare;  
  singleCell.forestShare = OforestShare + AforestShare;
  double deforestHa = singleCell.deforestA[Age]*singleCell.LandAreaHa;
  singleCell.deforestHaTot += deforestHa;
  double afforestHa = singleCell.forestAgeShare[Age]*singleCell.LandAreaHa;
  singleCell.afforestHaTot += afforestHa;
  singleCell.ProdLongA[Age] = defBiomass*data["FRACLONGPROD"].v()*(1-coeff.HarvLoos.v())*deforestHa;
  singleCell.ProdShortA[Age] = singleCell.ProdLongA[Age];
  singleCell.LitterA[Age] = data["CLITTERHA"].v()*deforestHa;
  singleCell.FineRootA[Age] = data["CBELOWHA"].v()*0.3*deforestHa;
  singleCell.SOCA[Age] = data["SOCHA"].v()*deforestHa;
  singleCell.SOCaffor[Age] = data["SOCHA"].v() * afforestHa;
  singleCell.LitterAffor[Age] = data["CLITTERHA"].v() * afforestHa;
//********* Emissions from deforestation *************
  for (int i=0;i<=Age;i++) {
// Calculate Emissions from deforestation in current cell for current year caused by decomposition
    EmissionsProductCur = EmissionsProductCur + singleCell.ProdLongA[i]*coeff.decRateL.v()
                        + singleCell.ProdShortA[i]*coeff.decRateS.v();
    EmissionsLitterCur = EmissionsLitterCur + singleCell.LitterA[i]*0.3*data["DECWOOD"].v(year)
                       + singleCell.LitterA[i]*0.7*data["DECHERB"].v() + singleCell.FineRootA[i]
                       * data["DECHERB"].v(year);
    if (singleCell.SOCA[i] >= data["SOCHA"].v()*0.6*singleCell.deforestA[i]*singleCell.LandAreaHa) {                                                 
      EmissionsSOCCur = EmissionsSOCCur + singleCell.SOCA[i] * data["DECSOC"].v(year);
    }  
// MG: Recalculate carbon pools
    singleCell.ProdLongA[i] = singleCell.ProdLongA[i]*(1-coeff.decRateL.v());
    singleCell.ProdShortA[i] = singleCell.ProdShortA[i]*(1-coeff.decRateS.v());
    singleCell.LitterA[i] = singleCell.LitterA[i]*0.3*(1-data["DECWOOD"].v())
                          + singleCell.LitterA[i]*0.7*(1-data["DECHERB"].v());
    singleCell.FineRootA[i] = singleCell.FineRootA[i]*(1-data["DECHERB"].v());
    singleCell.SOCA[i] = singleCell.SOCA[i]*(1-data["DECSOC"].v(year));
    if (singleCell.SOCA[i]<data["SOCHA"].v()*0.6*singleCell.deforestA[i]) {
      singleCell.SOCA[i] = data["SOCHA"].v()*0.6*singleCell.deforestA[i]*singleCell.LandAreaHa; // SOCagr;
    }
  }
//Emissions from deforestation in current cell for current year caused by burning  
  EmissionsSlashBurnCur = EmissionsSlashBurnCur + abBiomassO*data["SLASHBURN"].v(year)* deforestHa;
  EmissionsDeadBurnCur = EmissionsDeadBurnCur + data["CDEADHA"].v()* deforestHa;
  EmissionsCRootBurnCur = EmissionsCRootBurnCur + data["CBELOWHA"].v(year)*0.7 * deforestHa;
//Emissions for current cell summed over years
  singleCell.EmissionsProduct += EmissionsProductCur;
  singleCell.EmissionsLitter += EmissionsLitterCur;
  singleCell.EmissionsSOC += EmissionsSOCCur;
  singleCell.EmissionsSlashBurn += EmissionsSlashBurnCur;
  singleCell.EmissionsDeadBurn += EmissionsDeadBurnCur;
  singleCell.EmissionsCRootBurn += EmissionsCRootBurnCur;
//Total emissions in current cell for current year
  EmissionsCur = EmissionsProductCur+EmissionsLitterCur+EmissionsSOCCur+EmissionsSlashBurnCur
               + EmissionsDeadBurnCur+EmissionsCRootBurnCur;
//Total emissions in current cell summed over years       
  singleCell.EmissionsTot += EmissionsCur;    
  EmissionsCurCountry[Country] += EmissionsCur;
//*************** END Emissions from deforestation ****************
//*************** Afforestation "negative" emissions block ********
  for (int ia=0; ia<Age; ia++) {
    if (singleCell.forestAgeShare[ia]>0) {
      double abovePhCur = newCohort.getBm(Age-ia) * singleCell.forestAgeShare[ia]; // afforested biomass per ha * afforShare
      PlantPhytHaBmGr += abovePhCur;
      if (singleCell.LitterAffor[ia]<5*singleCell.forestAgeShare[ia]*singleCell.LandAreaHa) {
        double CurEmissionsLitterAfforCur = 0.95 * pow((1.-exp(-0.1*abovePhCur)),3)
                                * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa;
        EmissionsLitterAfforCur += CurEmissionsLitterAfforCur;
        singleCell.LitterAffor[ia] += CurEmissionsLitterAfforCur;
      }
      if (singleCell.SOCaffor[ia]<=data["SOCHA"].v()*1.4* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa) {
        if (data["POTVEG"].v()==4||data["POTVEG"].v()==6) {
          EmissionsSOCAfforCur += 0.04 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
                                / (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
                                * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa;
          singleCell.SOCaffor[ia] += 0.04 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
                                / (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3)
                                * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa;                   
        }
        if (data["POTVEG"].v()==8) {
          EmissionsSOCAfforCur += 0.2 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
                                / (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
                                * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa;
          singleCell.SOCaffor[ia] += 0.2 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
                                   / (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
                                   * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa; 
        } else {
          EmissionsSOCAfforCur += 0.35 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
                                / (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
                                * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa; 
          singleCell.SOCaffor[ia] += 0.35 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
                                   / (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
                                   * singleCell.forestAgeShare[ia]*singleCell.LandAreaHa;
        }  
      }           // End for if SOCaffor <= ...
    }             // End for if forestShare > 0
  }               // End for age loop
//     Belowground phytomass of planted forest = 18% for tropical, 22 for temperate and 25 for boreal of aboveground phytomass
//  PlantPhytHaBmGr = newCohort.getBm();
  if (data["POTVEG"].v()==1 || data["POTVEG"].v()==2){
    PlantPhytHaBlGr = PlantPhytHaBmGr*0.18;
  }
  if (data["POTVEG"].v()==6 || data["POTVEG"].v()==7) {PlantPhytHaBlGr = PlantPhytHaBmGr*0.25;}
  if (data["POTVEG"].v()==3 || data["POTVEG"].v()==4 || data["POTVEG"].v()==5||(data["POTVEG"].v()==8)) {
    PlantPhytHaBlGr = PlantPhytHaBmGr*0.22;
  }
  double CurPlantPhytHaBmGr=(PlantPhytHaBmGr-singleCell.prevPlantPhytHaBmGr) * singleCell.LandAreaHa;
  double CurPlantPhytHaBlGr=(PlantPhytHaBlGr-singleCell.prevPlantPhytHaBlGr) * singleCell.LandAreaHa;
// Emissions in current cell summed over years 
  singleCell.EmLitterAffor += EmissionsLitterAfforCur;
  singleCell.EmSOCAffor += EmissionsSOCAfforCur;
//     Total (negative) emissions from afforestation for current year and current cell
  EmissionsAfforCur = (CurPlantPhytHaBmGr+CurPlantPhytHaBlGr) + EmissionsLitterAfforCur+EmissionsSOCAfforCur;
  EmissionsCurAfforCountry[Country] += EmissionsAfforCur;
//Total emissions in current cell summed over years       
  singleCell.EmissionsAffor += EmissionsAfforCur;
//////////////////// END Afforestation "negative" emissions block                


// Output results for countries
    double harvestTotM3 = (harvWood*OforestShare+harvWoodNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa;

     CountriesNforCover.inc(Country,year,AforestShare * singleCell.LandAreaHa);
     CountriesAfforHaYear.inc(Country,year,afforestHa);          
     CountriesNforTotC.inc(Country,year,singleCell.EmissionsAffor);
     CountriesAfforCYear.inc(Country,year,EmissionsAfforCur);  
  
     CountriesOforCover.inc(Country,year,OforestShare * singleCell.LandAreaHa);
     CountriesOforBiomassC.inc(Country,year,abBiomassO*OforestShare*singleCell.LandAreaHa);
     CountriesOforTotC.inc(Country,year,singleCell.EmissionsTot);
     CountriesDeforHaYear.inc(Country,year,deforestHa);  
     CountriesDeforCYear.inc(Country,year,EmissionsCur);    

     CountriesWoodHarvestM3Year.inc(Country,year,harvestTotM3);
     CountriesWoodLoosCYear.inc(Country,year,(harvWoodLost*OforestShare+harvWoodLostNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa);

     harvestGrid.set(xi,yi,harvestTotM3);



/*
// MG: Writing to file at specified time points

//	    if ((year==2000) || (year==2001)|| (year==2002)|| (year==2003)|| (year==2004) || (year==2005) || (year==2006)
//		if ((year==2020) || (year==2030) || (year==2050)) {
//	    if ((year==2000) || (year==2005) 
//                         || (year==2010) || (year==2015) || (year==2020) || (year==2030) || (year==2050)
//                         || (year==2100)
//         )
    {
      ireportYear++;
      double deltaReportYear = (year-prevReportYear);
      if (deltaReportYear == 0.) {deltaReportYear = 1.;};
	  double DefCurHa = (prevOForShareRP - OforestShare) * LandAreaHa / deltaReportYear;
      if (DefCurHa < 0) {DefCurHa = 0;}
      double AffCurHa = (AforestShare - prevAForShareRP) * LandAreaHa / deltaReportYear;
      if (AffCurHa < 0) {AffCurHa = 0;}
      double EmissionsDeforCur = (EmissionsTot - EmissionsTotPrev)/deltaReportYear;
      double EmissionsAfforCur1 = (EmissionsAffor - EmissionsAfforPrev)/deltaReportYear;

// Country Arrays
//   CountriesNforCover 
//   CountriesNforBiomass 
//   CountriesAfforHaYear   
//   CountriesAfforCYear  
//   CountriesOforCover
//   CountriesOforBiomass
//   CountriesDeforHaYear 
//   CountriesDeforCYear  

      CountriesNforCover.inc(Country,year,AforestShare * LandAreaHa);
      CountriesAfforHaYear.inc(Country,year,AffCurHa);          
      prevReportYear = year;
      prevOForShareRP = OforestShare;
      prevAForShareRP = AforestShare;
      EmissionsTotPrev = EmissionsTot;
      EmissionsAfforPrev = EmissionsAffor;
		  
//		  prevOForShareRPNP = OfsNoPay;
//		  prevAForShareRPNP = AfsNoPay;		  
//		  prevPlantPhytHaBmGr = PlantPhytHaBmGr;

//cohortNew.aging(); // MG: we grow planted forest

    } // MG: End of output years loop

    AforestSharePrev = AforestShare;
    prevOForShare = OforestShare;
    
    prevPlantPhytHaBmGr = PlantPhytHaBmGr;
    prevPlantPhytHaBlGr = PlantPhytHaBlGr;
    
    cropland = crop;
    builtupland = builtup;


//        cout<<Country+1<<"\t"<<Age+1<<"\t"<<iprice<<"\t"<<EmissionsCurAfforCountry[Country]<<endl;
//      arrdef.set(Country,Age+1,iprice+1,EmissionsCurCountry[Country]);
//      arraff.set(Country,Age+1,iprice+1,EmissionsCurAfforCountry[Country]);
      
//      arrdefNP.set(Country,Age+1,iprice+1,EmissionsCurCountryNP[Country]);
//      arraffNP.set(Country,Age+1,iprice+1,EmissionsCurAfforCountryNP[Country]);
      


//cout<< "vIncr= " <<decision.vIncr()<< "    gTopt= "<< rotInterM<< "   rotInter()="<<decision.rotInter()<<endl;
//cout<<"harvWood= "<< harvWood<<"     HarvestVol()= "<<decision.woodHarvestVol()<<"     HarvWood_gHbm="<<fi.gHbm(rotInterM,decision.vIncr())<<endl;  	  
//cout<<"    NPV= "<<decision.forVal() <<endl;

//cout << "forFlag = " << forFlag << endl;


if (forFlag > 0)
{

cout <<xi<<"\t"<<yi;
cout <<"\t"<<Country;
cout <<"\t"<<year;
cout <<"\t" << rotationTimeCurr << "\t"<< Rotation <<"\t" << MAIRot  <<"\t";

//cout <<sawnThW*OforestShare*LandAreaHa<<"\t"<<restThW*OforestShare*LandAreaHa<<"\t";
//cout <<sawnW*OforestShare*LandAreaHa<<"\t"<<restW*OforestShare*LandAreaHa<<"\t"<<abBiomassO*OforestShare*LandAreaHa<<"\t";
//cout <<data["CABOVEHA"].v()*OforestShare*LandAreaHa<<"\t";

cout <<sawnThW*OforestShare*LandAreaHa*4;
cout <<"\t"<<restThW*OforestShare*LandAreaHa*4<<"\t";
cout <<sawnW*OforestShare*LandAreaHa*4;
cout <<"\t"<<restW*OforestShare*LandAreaHa*4<<"\t";

cout << sawnThWpot*OforestShare*LandAreaHa*4<<"\t" << restThWpot*OforestShare*LandAreaHa*4 << "\t";
cout << sawnWpot*OforestShare*LandAreaHa*4 <<"\t" << restWpot*OforestShare*LandAreaHa*4;


cout << "\t" << abBiomassO*OforestShare*LandAreaHa <<"\t";
cout << data["CABOVEHA"].v()*refForShare*LandAreaHa <<"\t";

cout << maiForest.get(xi,yi) << "\t" <<fval << "\t" << forNPV_CH;
cout << "\t" << TimberPrice << "\t";


cout <<sawnThW<<"\t"<<restThW<<"\t";
cout <<sawnW<<"\t"<<restW<<"\t"<<abBiomassO<<"\t";

cout <<data["CABOVEHA"].v()<<"\t";

cout << refForShare*LandAreaHa<<"\t";

cout<< int(managedForest.get(xi,yi))<<"\t";

cout<< int(thinningForest.get(xi,yi))<<"\t";

cout<< EmissionsTot<<"\t";
cout<<EmissionsAffor<<"\t";

cout<<harvWoodLost*OforestShare*LandAreaHa;
cout << "\t" << OforestShare * LandAreaHa;  // Remained old forest, ha
cout << "\t" << AforestShare * LandAreaHa; // Planted forest, ha (accumulated)

cout << "\t" << PlantPhytHaBmGr*LandAreaHa;
cout << "\t" << EmLitterAffor;
cout << "\t" << EmSOCAffor;
cout << "\t" << harvWoodNew;
cout << "\t" << harvWoodLostNew;
cout<< endl;
}
*/
 }



