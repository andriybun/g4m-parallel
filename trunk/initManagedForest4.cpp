void initManagedForest(dataDetStruct &data_all, g4m::incrementTab &fi, datGlobal &dat_all,
                       griddata &maiForest, griddata2<char> &thinningForest,
                       griddata2<char> &rotationType, griddata2<char> &managedForest,
                       griddata &rotationForest, griddata &harvestGrid) 
 {
  double woodHarvest[209];
  double woodLost[209];
//  double woodHarvestStat[209];
  int managedCount[209];

  for (int i=0; i<=208; i++){
    woodHarvest[i]=0.; 
    woodLost[i]=0.;
//    woodHarvestStat[i]=0.;
    managedCount[i]=0;
  }

//  woodHarvestStat[11-1]=17318000.;     // Austria
//  woodHarvestStat[25-1]=368706000.;    // Brasil
//  woodHarvestStat[33-1]=195869000.;    // Canada
//  woodHarvestStat[38-1]=159081000.;    // China
//  woodHarvestStat[61-1]=47203000.;     // Finland
//  woodHarvestStat[62-1]=55621000.;     // France
//  woodHarvestStat[69-1]=42177000.;     // Germany
//  woodHarvestStat[156-1]=336527000.;   // Russia
//  woodHarvestStat[165-1]=5545000.;     // Slovakia
//  woodHarvestStat[179-1]=58140000.;    // Sweden
//  woodHarvestStat[197-1]=596920000.;   // US
  
  double sawnW = 0.;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
  double restW = 0.;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
  double sawnThW = 0.;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
  double restThW = 0.;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
  double defIncome = 0.;

  dataDetStruct::iterator iter = data_all.begin();
  while (iter != data_all.end()) {
    if (iter->PROTECT[2000] == 0) {
      if (  
            (iter->COUNTRY[2000] == 10) 
          || (iter->COUNTRY[2000] == 11)
          || (iter->COUNTRY[2000] == 17)          
          || (iter->COUNTRY[2000] == 25)
          || (iter->COUNTRY[2000] == 27)
          || (iter->COUNTRY[2000] == 30)
          || (iter->COUNTRY[2000] == 33)        
          || (iter->COUNTRY[2000] == 43)  
          || (iter->COUNTRY[2000] == 46)    
          || (iter->COUNTRY[2000] == 47)  
          || (iter->COUNTRY[2000] == 56)                        
          || (iter->COUNTRY[2000] == 61)
          || (iter->COUNTRY[2000] == 62)
          || (iter->COUNTRY[2000] == 69) 
          || (iter->COUNTRY[2000] == 71) 
          || (iter->COUNTRY[2000] == 82)    
          || (iter->COUNTRY[2000] == 83)   
          || (iter->COUNTRY[2000] == 89)  
          || (iter->COUNTRY[2000] == 92)    
          || (iter->COUNTRY[2000] == 96)   
          || (iter->COUNTRY[2000] == 107) 
          || (iter->COUNTRY[2000] == 112)    
          || (iter->COUNTRY[2000] == 113)
          || (iter->COUNTRY[2000] == 114)          
          || (iter->COUNTRY[2000] == 128)          
          || (iter->COUNTRY[2000] == 135)           
          || (iter->COUNTRY[2000] == 137)            
          || (iter->COUNTRY[2000] == 142)    
          || (iter->COUNTRY[2000] == 150)   
          || (iter->COUNTRY[2000] == 151) 
          || (iter->COUNTRY[2000] == 155)                                                                                                             
          || (iter->COUNTRY[2000] == 156)             
          || (iter->COUNTRY[2000] == 165) 
          || (iter->COUNTRY[2000] == 166)       
          || (iter->COUNTRY[2000] == 170)                         
          || (iter->COUNTRY[2000] == 179) 
          || (iter->COUNTRY[2000] == 180) 
          || (iter->COUNTRY[2000] == 190)        
          || (iter->COUNTRY[2000] == 194)    
          || (iter->COUNTRY[2000] == 196)                                               
          || (iter->COUNTRY[2000] == 197) ) { // Test only some countries
        map<string, interpol> data;
        fillContainer(*iter,data);     
//        xi = int((data["X"].v() - 0.25 + 180)*2);
//        yi = int((data["Y"].v() - 0.25 + 90)*2);    
    	int xi = (iter->x);
    	int yi = (iter->y);
        double X = (iter->x)*0.5+0.25-180;
        double Y = (iter->y)*0.5+0.25-90;
        int Country0 = (int)data["COUNTRY"].v();
        double LandAreaHa = data["LANDAREA"].v()*100;
        double forestArea0 = LandAreaHa * data["FOREST"].v();
        int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
        int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)    
        double harvWood=0.; //MG: harvestable wood, m3
        double abBiomassO = 0.;
        double MAI = 0.;
//*******************************************************************************
//********************************************************************************
//**********************************************************************************
//Setting up forest with initial biomass in the grid, "potential" rotation time
// for getting this biomass
// using Georg's FM tool
//------------------
//MG: From Georg's "Optimal rotation time"
    //Get optimal rotation time
    //0 .. Rotation time Unmanaged forests       // Very long RotTime - untill the forest breaks
    //1 .. Highest average increment             // Short RotTime - approx = age of max MAI
    //2 .. Maximum avarage Biomass               // Very long RotTime - approx = age when forest accumulates max biomass
    //3 .. Maximum average Biomass with thinning // Very long (longer than 2) RotTime - approx = age when forest accumulates max biomass
    //4 .. Maximum harvest at final cut          // Long RotTime (approx 2 times longer than at 1 but almost 2 times shorter than at 2)
    //5 .. Maximum average harvest with final cut // Very short RotTime (shorter than at 1)

        int optimUnmanaged = 0;
        int optimMAI = 1;
        int optimMaxBm = 2;
        int optimMaxBmTh = 3;
        int optimHarvFin = 4;
        int optimHarvAve = 5;
        int rotUnmanaged = 0;
        int rotMAI = 0;
        int rotMaxBm = 0;
        int rotMaxBmTh = 0;
        int rotHarvFin = 0;
        int rotHarvAve = 0;

        g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
        g4m::ipol<double,double> hlv;  //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
        g4m::ipol<double,double> hle;  //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
        g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
        g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)

        sws.insert(10, .0);
        sws.insert(30, .6);
        hlv.insert(0, .0);
        hlv.insert(25, .7);
        hle.insert(0, .0);
        hle.insert(25, .7);
        dbv.insert(0, 2);
        dbe.insert(0, 3);
 
        g4m::ageStruct::v res; // MG: results vector for the set (old) forest    

// -- Initialise DIMA to get MAI
        data["PRICEC"].insert(0,0.);
        {
        dima decision((int)data["BYEAR"].v()
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
		       , data["FOREST"].v((int)data["BYEAR"].v())
		       , wprice["regprice"]
		       , wprice["regprice0"].v()
               , biomasRot+1
               , harvWood);
               MAI = decision.vIncr() / data["FTIMBER"].v(); //MG: mean annual increment in tC/ha/year
           }
        if (MAI < 0.){MAI = 0.;}     
        maiForest.set(xi,yi,MAI);
        double harvMAI = MAI*data["FTIMBER"].v()*(1-coeff.HarvLoos.v());


        double forFlag = 0.;          //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
        if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && MAI > 0) {
          biomasRot = fi.gU(data["CABOVEHA"].v(), MAI, 1);        // rotation time to get current biomass (without thinning)
          biomasRotTh = fi.gUt(data["CABOVEHA"].v(), MAI, 1);     // rotation time to get current biomass (with thinning)     
          rotMAI = fi.gTopt(MAI, optimMAI);
          rotMaxBm = fi.gTopt(MAI, optimMaxBm);                        
          rotMaxBmTh = fi.gTopt(MAI, optimMaxBmTh);
          forFlag = 1.0;
        }     
        double Thinning = -1.;
        double ThinningInit = -1.;
        int Rotation = 0;
        int RotationInit = 0;
        if ((data["POPDENS"].v(1990) >0) && (data["GDP"].v(1990) > 0)) {
          ThinningInit = 1.;
          thinningForest.set(xi,yi,1.);
          RotationInit = biomasRotTh+1;
          rotationType.set(xi,yi,11);
          g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, MAI, RotationInit, ThinningInit,forFlag, 0.75);
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
          dima decision(1990
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
		       , data["FOREST"].v(1990)
		       , wprice["regprice"]
		       , wprice["regprice0"].v(2000)
               , rotMAI
               , harvMAI); 
//               , Rotation
//               , harvWood );
          abBiomassO = cohort.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * data["FTIMBER"].v()
		         * (1. -coeff.HarvLoos.v()));
	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
		            + sDefIncome * data["SLASHBURN"].v();

          if (MAI > MAI_CountryUprotect[Country0-1]) {
            if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome)) {
              if ((woodHarvest[Country0-1] < woodHarvestStat[Country0-1])) {
                managedForest.set(xi,yi,3.);
                thinningForest.set(xi,yi,1.);
                Rotation = rotMAI;
                rotationType.set(xi,yi,1);
              } else {
                managedForest.set(xi,yi,3.);
                Rotation = rotMaxBmTh;
                rotationType.set(xi,yi,3);
              }
            } else {
              if ((woodHarvest[Country0-1] < woodHarvestStat[Country0-1])) {
                managedForest.set(xi,yi,2.);
                Rotation = biomasRotTh+1;
                rotationType.set(xi,yi,11);
              } else {
                managedForest.set(xi,yi,2.);
                Rotation = rotMaxBmTh;
                rotationType.set(xi,yi,3);
              }
            }
          } else {
            if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome)) {
              if ((woodHarvest[Country0-1] < woodHarvestStat[Country0-1])) {
                managedForest.set(xi,yi,2.);
                thinningForest.set(xi,yi,1.);
                Rotation = biomasRotTh+1;
                rotationType.set(xi,yi,11);
              } else {
                managedForest.set(xi,yi,2.);
                Rotation = rotMaxBmTh;
                rotationType.set(xi,yi,3);
              }
            } else {
              if ((woodHarvest[Country0-1] < woodHarvestStat[Country0-1])) {
                managedForest.set(xi,yi,1.);
                Rotation = biomasRotTh+1;
                rotationType.set(xi,yi,11);
              } else {
                managedForest.set(xi,yi,1.);
                Rotation = rotMaxBmTh;
                rotationType.set(xi,yi,3);
              }
            }
          }
          cohort.setRotPeriod(Rotation);
//          cohort.setStockingdegree(thinningForest.get(xi,yi));
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
/*          coeff.PriceC.insert(0,0.);
          dima decision(1990
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
		       , data["FTIMBER"].v() 
		       , coeff.PriceTimberMaxR
		       , coeff.PriceTimberMinR
		       , coeff.FCuptake
		       , data["GDP"]
		       , coeff.HarvLoos
		       , data["FOREST"].v(1990)
		       , wprice["regprice"]
		       , wprice["regprice0"].v(2000)
               , Rotation
               , harvWood );
          abBiomassO = cohort.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * coeff.FTimber.v()* (1. -coeff.HarvLoos.v()));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * coeff.FTimber.v()
		         * (1. -coeff.HarvLoos.v()));
	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
		            + sDefIncome * data["SLASHBURN"].v(); */
        } else {
          ThinningInit = -1.;
          thinningForest.set(xi,yi,-1.);
          RotationInit = biomasRot+1;
          rotationType.set(xi,yi,10);
          g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, MAI, RotationInit, ThinningInit,forFlag, 0.75);  
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
          dima decision(1990
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
		       , data["FOREST"].v(1990)
		       , wprice["regprice"]
		       , wprice["regprice0"].v(2000)
               , rotMAI
               , harvMAI); 
//               , Rotation
//               , harvWood );
          abBiomassO = cohort.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * data["FTIMBER"].v()
		         * (1. -coeff.HarvLoos.v()));
	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
		            + sDefIncome * data["SLASHBURN"].v();

          if (MAI > MAI_CountryUprotect[Country0-1]) {
            if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome)) {
              managedForest.set(xi,yi,0.);
              Rotation = biomasRot+1;
              rotationType.set(xi,yi,1);
            } else {
              managedForest.set(xi,yi,-1.);
              Rotation = biomasRot+1;
              rotationType.set(xi,yi,10);
            }
          } else {
            if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome)) {
              managedForest.set(xi,yi,-1.);
              Rotation = biomasRot+1;
              rotationType.set(xi,yi,10);
            } else {      
              managedForest.set(xi,yi,-2.);
              Rotation = biomasRot+1;
              rotationType.set(xi,yi,10);
            }
          }
          cohort.setRotPeriod(Rotation);
//          cohort.setStockingdegree(thinningForest.get(xi,yi));
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
/*          coeff.PriceC.insert(0,0.);
          dima decision(1990
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
		       , data["FTIMBER"].v() 
		       , coeff.PriceTimberMaxR
		       , coeff.PriceTimberMinR
		       , coeff.FCuptake
		       , data["GDP"]
		       , coeff.HarvLoos
		       , data["FOREST"].v(1990)
		       , wprice["regprice"]
		       , wprice["regprice0"].v(2000)
               , Rotation
               , harvWood );
          abBiomassO = cohort.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * coeff.FTimber.v()* (1. -coeff.HarvLoos.v()));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * coeff.FTimber.v()
		         * (1. -coeff.HarvLoos.v()));
	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
		            + sDefIncome * data["SLASHBURN"].v(); */
        }
//        if (managedForest.get(xi,yi) > 0) {
        if (thinningForest.get(xi,yi) > 0) {          
          woodHarvest[Country0-1] += harvWood * forestArea0;
          managedCount[Country0-1] +=1;
        } else {
          woodLost[Country0-1] += harvWood * forestArea0;
        }
        rotationForest.set(xi,yi,Rotation);	
      }        // End for COUNTRY test
    }          // End if (iter->PROTECT[2000] == 0)
    iter++;
  }            // End while(iter != data_all.end())
  cout << "end of first pass" << endl;
//******************************************************************************
//**************************Second Pass********************
//******************************************************************************

  iter = data_all.begin();
  while (iter != data_all.end()) {
    int Country0 = 1;
    double HarvestTmp = 0;
    double newHarvestTmp = 0;
    double forestArea0 = 0.;
    double abBiomassO = 0.;  	 
    if (iter->PROTECT[2000]==0) {   // We consider only unprotected land
      if (  
            (iter->COUNTRY[2000] == 10) 
          || (iter->COUNTRY[2000] == 11)
          || (iter->COUNTRY[2000] == 17)          
          || (iter->COUNTRY[2000] == 25)
          || (iter->COUNTRY[2000] == 27)
          || (iter->COUNTRY[2000] == 30)
          || (iter->COUNTRY[2000] == 33)        
          || (iter->COUNTRY[2000] == 43)  
          || (iter->COUNTRY[2000] == 46)    
          || (iter->COUNTRY[2000] == 47)  
          || (iter->COUNTRY[2000] == 56)                        
          || (iter->COUNTRY[2000] == 61)
          || (iter->COUNTRY[2000] == 62)
          || (iter->COUNTRY[2000] == 69) 
          || (iter->COUNTRY[2000] == 71) 
          || (iter->COUNTRY[2000] == 82)    
          || (iter->COUNTRY[2000] == 83)   
          || (iter->COUNTRY[2000] == 89)  
          || (iter->COUNTRY[2000] == 92)    
          || (iter->COUNTRY[2000] == 96)   
          || (iter->COUNTRY[2000] == 107) 
          || (iter->COUNTRY[2000] == 112)    
          || (iter->COUNTRY[2000] == 113)
          || (iter->COUNTRY[2000] == 114)          
          || (iter->COUNTRY[2000] == 128)          
          || (iter->COUNTRY[2000] == 135)           
          || (iter->COUNTRY[2000] == 137)            
          || (iter->COUNTRY[2000] == 142)    
          || (iter->COUNTRY[2000] == 150)   
          || (iter->COUNTRY[2000] == 151) 
          || (iter->COUNTRY[2000] == 155)                                                                                                             
          || (iter->COUNTRY[2000] == 156)             
          || (iter->COUNTRY[2000] == 165) 
          || (iter->COUNTRY[2000] == 166)       
          || (iter->COUNTRY[2000] == 170)                         
          || (iter->COUNTRY[2000] == 179) 
          || (iter->COUNTRY[2000] == 180) 
          || (iter->COUNTRY[2000] == 190)        
          || (iter->COUNTRY[2000] == 194)    
          || (iter->COUNTRY[2000] == 196)                                               
          || (iter->COUNTRY[2000] == 197) ) { // Test only some countries
        do {
          map<string, interpol> data;
          fillContainer(*iter, data);
          int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
          int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
          double harvWood=0; //MG: harvestable wood, m3
          int xi = (iter->x);
          int yi = (iter->y);
          double X = (iter->x)*0.5+0.25-180;
          double Y = (iter->y)*0.5+0.25-90;
          Country0 = (int)data["COUNTRY"].v();
          double LandAreaHa = data["LANDAREA"].v()*100;
          forestArea0 = LandAreaHa * data["FOREST"].v();
          double harvMAI = maiForest.get(xi,yi)*data["FTIMBER"].v()*(1-coeff.HarvLoos.v());
          int optimUnmanaged = 0;
          int optimMAI = 1;
          int optimMaxBm = 2;
          int optimMaxBmTh = 3;
          int optimHarvFin = 4;
          int optimHarvAve = 5;
          int rotUnmanaged = 0;
          int rotMAI = 0;
          int rotMaxBm = 0;
          int rotMaxBmTh = 0;
          int rotHarvFin = 0;
          int rotHarvAve = 0;

          g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
          g4m::ipol<double,double> hlv;  //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
          g4m::ipol<double,double> hle;  //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
          g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
          g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)

          sws.insert(10, .0);
          sws.insert(30, .6);
          hlv.insert(0, .0);
          hlv.insert(25, .7);
          hle.insert(0, .0);
          hle.insert(25, .7);
          dbv.insert(0, 2);
          dbe.insert(0, 3);
          g4m::ageStruct::v res;  // MG: results vector for the set (old) forest 
          int Rotation = 0;
          double forFlag = 0.;    //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
          if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0) {
            biomasRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);       // rotation time to get current biomass (without thinning)
            biomasRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
            if (thinningForest.get(xi,yi) == 1) {
              Rotation = biomasRotTh+1; 
            }  else {
              Rotation = biomasRot+1;
            }
            rotMAI = fi.gTopt(maiForest.get(xi,yi), optimMAI);
            rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), optimMaxBmTh);         
            forFlag = 1.0;
          }     
          g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75); 
          cohort.aging();
          cohort.setRotPeriod(rotationForest.get(xi,yi));
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v() ;
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
          coeff.PriceC.insert(0,0.);
          dima decision(1990
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
		       , data["FOREST"].v(1990)
		       , wprice["regprice"]
		       , wprice["regprice0"].v(2000)
               , rotMAI
               , harvMAI); 
//               , Rotation
//               , harvWood );
          abBiomassO = cohort.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * data["FTIMBER"].v() 
		         * (1. -coeff.HarvLoos.v()));
	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
		            + sDefIncome * data["SLASHBURN"].v();
		            
          harvestGrid.set(xi,yi,harvWood * forestArea0);   
          if (managedForest.get(xi,yi) > 0) {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}
            if (woodHarvest[Country0-1] <= (0.9 * woodHarvestStat[Country0-1])) {
              if ((Rotation > rotMAI) && (rotationForest.get(xi,yi) > Rotation)) { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,Rotation);		   
                cohort.setRotPeriod(Rotation);
                res = cohort.aging();
                sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v() ;
                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//                coeff.PriceC.insert(0,0.);
//                dima decision(1990
//	                 , data["NPP"]
//		             , data["SPOPDENS"]
//		             , data["SAGRSUIT"]
//		             , data["PRICEINDEX"]
//		             , coeff.PriceIndexE
//		             , data["R"]
//		             , coeff.PriceC
//		             , coeff.PlantingCostsR
//		             , coeff.PriceLandMinR
//		             , coeff.PriceLandMaxR
//		             , coeff.MaxRotInter
//		             , coeff.MinRotInter
//		             , coeff.decRateL
//		             , coeff.decRateS
//		             , data["FRACLONGPROD"]
//		             , coeff.baseline
//		             , data["FTIMBER"] 
//		             , coeff.PriceTimberMaxR
//		             , coeff.PriceTimberMinR
//		             , coeff.FCuptake
//		             , data["GDP"]
//		             , coeff.HarvLoos
//		             , data["FOREST"].v(1990)
//		             , wprice["regprice"]
//		             , wprice["regprice0"].v(2000)
//                     , Rotation
//                     , harvWood );
//                abBiomassO = cohort.getBm();
//	            double pDefIncome = abBiomassO * 
//                         (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	            //Immediate Pay if deforested (Slash and Burn)
//                double sDefIncome = abBiomassO *
//		                 (decision.priceTimber() * data["FTIMBER"].v() 
//                       * (1. -coeff.HarvLoos.v()));
//                defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		                  + sDefIncome * data["SLASHBURN"].v();
         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              } else if (rotationForest.get(xi,yi) > rotMAI) { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,rotMAI);		   
                cohort.setRotPeriod(rotMAI);
                res = cohort.aging();
                sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v() ;
                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//                coeff.PriceC.insert(0,0.);
//                dima decision(1990
//	                 , data["NPP"]
//		             , data["SPOPDENS"]
//		             , data["SAGRSUIT"]
//		             , data["PRICEINDEX"]
//		             , coeff.PriceIndexE
//		             , data["R"]
//		             , coeff.PriceC
//		             , coeff.PlantingCostsR
//		             , coeff.PriceLandMinR
//		             , coeff.PriceLandMaxR
//		             , coeff.MaxRotInter
//		             , coeff.MinRotInter
//		             , coeff.decRateL
//		             , coeff.decRateS
//		             , data["FRACLONGPROD"]
//		             , coeff.baseline
//		             , data["FTIMBER"] 
//		             , coeff.PriceTimberMaxR
//		             , coeff.PriceTimberMinR
//		             , coeff.FCuptake
//		             , data["GDP"]
//		             , coeff.HarvLoos
//		             , data["FOREST"].v(1990)
//		             , wprice["regprice"]
//		             , wprice["regprice0"].v(2000)
//                     , Rotation
//                     , harvWood );
//                abBiomassO = cohort.getBm();
//	            double pDefIncome = abBiomassO * 
//                         (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	            //Immediate Pay if deforested (Slash and Burn)
//                double sDefIncome = abBiomassO *
//		                 (decision.priceTimber() * data["FTIMBER"].v() 
//                       * (1. -coeff.HarvLoos.v()));
//                defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		                  + sDefIncome * data["SLASHBURN"].v();
         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              }            
            } else if (woodHarvest[Country0-1] >= 1.1 * woodHarvestStat[Country0-1]) {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}                   
              if (rotationForest.get(xi,yi) < Rotation) { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,Rotation);		   
                cohort.setRotPeriod(Rotation);
                res = cohort.aging();
                sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v() ;
                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//                coeff.PriceC.insert(0,0.);
//                dima decision(1990
//	                 , data["NPP"]
//		             , data["SPOPDENS"]
//		             , data["SAGRSUIT"]
//		             , data["PRICEINDEX"]
//		             , coeff.PriceIndexE
//		             , data["R"]
//		             , coeff.PriceC
//		             , coeff.PlantingCostsR
//		             , coeff.PriceLandMinR
//		             , coeff.PriceLandMaxR
//		             , coeff.MaxRotInter
//		             , coeff.MinRotInter
//		             , coeff.decRateL
//		             , coeff.decRateS
//		             , data["FRACLONGPROD"]
//		             , coeff.baseline
//		             , data["FTIMBER"] 
//		             , coeff.PriceTimberMaxR
//		             , coeff.PriceTimberMinR
//		             , coeff.FCuptake
//		             , data["GDP"]
//		             , coeff.HarvLoos
//		             , data["FOREST"].v(1990)
//		             , wprice["regprice"]
//		             , wprice["regprice0"].v(2000)
//                     , Rotation
//                     , harvWood );
//                abBiomassO = cohort.getBm();
//	            double pDefIncome = abBiomassO * 
//                         (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	            //Immediate Pay if deforested (Slash and Burn)
//                double sDefIncome = abBiomassO *
//		                 (decision.priceTimber() * data["FTIMBER"].v() 
//                       * (1. -coeff.HarvLoos.v()));
//                defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		                  + sDefIncome * data["SLASHBURN"].v();
         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              } else if (rotationForest.get(xi,yi) < rotMaxBmTh) {
                HarvestTmp = harvestGrid.get(xi,yi);
			    rotationForest.set(xi,yi, rotMaxBmTh);
                cohort.setRotPeriod(rotMaxBmTh);
                res = cohort.aging();
                sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v() ;
                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//                coeff.PriceC.insert(0,0.);
//                dima decision(1990
//	                 , data["NPP"]
//		             , data["SPOPDENS"]
//		             , data["SAGRSUIT"]
//		             , data["PRICEINDEX"]
//		             , coeff.PriceIndexE
//		             , data["R"]
//		             , coeff.PriceC
//		             , coeff.PlantingCostsR
//		             , coeff.PriceLandMinR
//		             , coeff.PriceLandMaxR
//		             , coeff.MaxRotInter
//		             , coeff.MinRotInter
//		             , coeff.decRateL
//		             , coeff.decRateS
//		             , data["FRACLONGPROD"]
//		             , coeff.baseline
//		             , data["FTIMBER"] 
//		             , coeff.PriceTimberMaxR
//		             , coeff.PriceTimberMinR
//		             , coeff.FCuptake
//		             , data["GDP"]
//		             , coeff.HarvLoos
//		             , data["FOREST"].v(1990)
//		             , wprice["regprice"]
//		             , wprice["regprice0"].v(2000)
//                     , Rotation
//                     , harvWood );
//                abBiomassO = cohort.getBm();
//	            double pDefIncome = abBiomassO * 
//                         (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	            //Immediate Pay if deforested (Slash and Burn)
//                double sDefIncome = abBiomassO *
//		                 (decision.priceTimber() * data["FTIMBER"].v() 
//                       * (1. -coeff.HarvLoos.v()));
//                defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		                  + sDefIncome * data["SLASHBURN"].v();
  			    newHarvestTmp = harvWood * forestArea0;
			    harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              }
            }
          }
          iter++;
        } while (((iter-1)->COUNTRY[2000]) == ((iter)->COUNTRY[2000]));   // Check are we in the same country  // end for Within current country
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\twoodHarvestStat= "<<woodHarvestStat[Country0-1] <<endl;}
      } else {iter++;}  // end for Country select
    } else {iter++;}  // end for IF unprotected
  } //end for WHILE
//************************End of Second Pass************************************

cout << "end of second pass" << endl;

//
///////////////////////////////////////////////////
////
////                   Third Pass
////
///////////////////////////////////////////////////
//
  iter = data_all.begin();

//cout << "Putting data for current cell into conteiner... "<< endl;
   while (iter != data_all.end())
   {
	if (iter->PROTECT[2000] == 0)
	{
      if (  
            (iter->COUNTRY[2000] == 10) 
          || (iter->COUNTRY[2000] == 11)
          || (iter->COUNTRY[2000] == 17)          
          || (iter->COUNTRY[2000] == 25)
          || (iter->COUNTRY[2000] == 27)
          || (iter->COUNTRY[2000] == 30)
          || (iter->COUNTRY[2000] == 33)        
          || (iter->COUNTRY[2000] == 43)  
          || (iter->COUNTRY[2000] == 46)    
          || (iter->COUNTRY[2000] == 47)  
          || (iter->COUNTRY[2000] == 56)                        
          || (iter->COUNTRY[2000] == 61)
          || (iter->COUNTRY[2000] == 62)
          || (iter->COUNTRY[2000] == 69) 
          || (iter->COUNTRY[2000] == 71) 
          || (iter->COUNTRY[2000] == 82)    
          || (iter->COUNTRY[2000] == 83)   
          || (iter->COUNTRY[2000] == 89)  
          || (iter->COUNTRY[2000] == 92)    
          || (iter->COUNTRY[2000] == 96)   
          || (iter->COUNTRY[2000] == 107) 
          || (iter->COUNTRY[2000] == 112)    
          || (iter->COUNTRY[2000] == 113)
          || (iter->COUNTRY[2000] == 114)          
          || (iter->COUNTRY[2000] == 128)          
          || (iter->COUNTRY[2000] == 135)           
          || (iter->COUNTRY[2000] == 137)            
          || (iter->COUNTRY[2000] == 142)    
          || (iter->COUNTRY[2000] == 150)   
          || (iter->COUNTRY[2000] == 151) 
          || (iter->COUNTRY[2000] == 155)                                                                                                             
          || (iter->COUNTRY[2000] == 156)             
          || (iter->COUNTRY[2000] == 165) 
          || (iter->COUNTRY[2000] == 166)       
          || (iter->COUNTRY[2000] == 170)                         
          || (iter->COUNTRY[2000] == 179) 
          || (iter->COUNTRY[2000] == 180) 
          || (iter->COUNTRY[2000] == 190)        
          || (iter->COUNTRY[2000] == 194)    
          || (iter->COUNTRY[2000] == 196)                                               
          || (iter->COUNTRY[2000] == 197) ) { // Test only some countries
        map<string, interpol> data;
        fillContainer(*iter,data);     
    	   //       xi = int((data["X"].v() - 0.25 + 180)*2);
    	   //       yi = int((data["Y"].v() - 0.25 + 90)*2);    
    	            int xi = (iter->x);
    	            int yi = (iter->y);
    	            double X = (iter->x)*0.5+0.25-180;
    	            double Y = (iter->y)*0.5+0.25-90;
    	          
//cout << "Xi = "<< xi <<"   Yi = "<< yi << endl;       
       int Country0 = (int)data["COUNTRY"].v();
                       

       double LandAreaHa = data["LANDAREA"].v()*100;
       double forestArea0 = LandAreaHa * data["FOREST"].v();
              
    int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
    int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)    
    double harvWood=0; //MG: harvestable wood, m3
    double abBiomassO = 0.;
    double HarvestTmp = 0;
    double newHarvestTmp = 0;




//*******************************************************************************
//********************************************************************************
//**********************************************************************************
//Setting up forest with initial biomass in the grid, "potential" rotation time
// for getting this biomass
// using Georg's FM tool
//------------------
//MG: From Georg's "Optimal rotation time"
    //Get optimal rotation time
    //0 .. Rotation time Unmanaged forests       // Very long RotTime - untill the forest breaks
    //1 .. Highest average increment             // Short RotTime - approx = age of max MAI
    //2 .. Maximum avarage Biomass               // Very long RotTime - approx = age when forest accumulates max biomass
    //3 .. Maximum average Biomass with thinning // Very long (longer than 2) RotTime - approx = age when forest accumulates max biomass
    //4 .. Maximum harvest at final cut          // Long RotTime (approx 2 times longer than at 1 but almost 2 times shorter than at 2)
    //5 .. Maximum average harvest with final cut // Very short RotTime (shorter than at 1)

//    int optRot = 0;  //MG: optimal rotation time (see above)

//#include "FM_initialization.cpp"
/////////////////////////////////// "FM_initialization
//MG: From Georg's "Optimal rotation time"
    //Get optimal rotation time
    //0 .. Rotation time Unmanaged forests       // Very long RotTime - untill the forest breaks
    //1 .. Highest average increment             // Short RotTime - approx = age of max MAI
    //2 .. Maximum avarage Biomass               // Very long RotTime - approx = age when forest accumulates max biomass
    //3 .. Maximum average Biomass with thinning // Very long (longer than 2) RotTime - approx = age when forest accumulates max biomass
    //4 .. Maximum harvest at final cut          // Long RotTime (approx 2 times longer than at 1 but almost 2 times shorter than at 2)
    //5 .. Maximum average harvest with final cut // Very short RotTime (shorter than at 1)

        int optimUnmanaged = 0;
        int optimMAI = 1;
        int optimMaxBm = 2;
        int optimMaxBmTh = 3;
        int optimHarvFin = 4;
        int optimHarvAve = 5;
        int rotUnmanaged = 0;
        int rotMAI = 0;
        int rotMaxBm = 0;
        int rotMaxBmTh = 0;
        int rotHarvFin = 0;
        int rotHarvAve = 0;

        g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
        g4m::ipol<double,double> hlv;  //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
        g4m::ipol<double,double> hle;  //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
        g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
        g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)

        sws.insert(10, .0);
        sws.insert(30, .6);
        hlv.insert(0, .0);
        hlv.insert(25, .7);
        hle.insert(0, .0);
        hle.insert(25, .7);
        dbv.insert(0, 2);
        dbe.insert(0, 3);
 
        g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
// End of FM initialisation

        double MAI = maiForest.get(xi,yi); //MG: mean annual increment in tC/ha/year
        double harvMAI = MAI*data["FTIMBER"].v()*(1-coeff.HarvLoos.v());



//cout << "xi= "<<xi<<"\t yi= "<<yi<<"\t MAI= "<< MAI<<"\t mai= "<< maiForest.get(xi,yi)<<"\t NPP= "<<iter->NPP[2000]<<"\t dNPP= "<<data["NPP"].v()<<endl; 
  

  double forFlag = 0.; //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
//  if (data["FOREST"].v(1990) >0 && data["BIOMASS"].v() > 0 && MAI > 0)
  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && MAI > 0)  
  {
//          biomasRot = fi.gU(data["BIOMASS"].v(), MAI, 1);       // rotation time to get current biomass (without thinning)
//          biomasRotTh = fi.gUt(data["BIOMASS"].v(), MAI, 1);     // rotation time to get current biomass (with thinning)     
          biomasRot = fi.gU(data["CABOVEHA"].v(), MAI, 1);       // rotation time to get current biomass (without thinning)
          biomasRotTh = fi.gUt(data["CABOVEHA"].v(), MAI, 1);     // rotation time to get current biomass (with thinning)     
          
          
//          rotUnmanaged = fi.gTopt(MAI, optimUnmanaged);
          rotMAI = fi.gTopt(MAI, optimMAI);
          rotMaxBm = fi.gTopt(MAI, optimMaxBm);                        
          rotMaxBmTh = fi.gTopt(MAI, optimMaxBmTh);
//          rotHarvAve = fi.gTopt(MAI, optimHarvAve);
          forFlag = 1.0;
  }     


// -- Initialise DIMA 
// #include "DIMA_initialization.cpp"    
        data["PRICEC"].insert(0,0.);
        dima decision((int)data["BYEAR"].v()
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
		       , data["FOREST"].v((int)data["BYEAR"].v())
		       , wprice["regprice"]
		       , wprice["regprice0"].v()
               , rotMAI
               , harvMAI);
// End of DIMA initialisation




//  g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
 double Thinning = -1.;
 double ThinningInit = -1.;
 int Rotation = 0;
 int RotationInit = 0;

//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}

if (woodHarvest[Country0-1] <= (0.9 * woodHarvestStat[Country0-1])) 
 {
//if (Country0 == 61){cout<<"Inside IF!!???<"<<endl;};
  if (thinningForest.get(xi,yi) < 0)
   {
                ThinningInit = 1.;
                thinningForest.set(xi,yi,1.);
                RotationInit = biomasRotTh+1;
                rotationType.set(xi,yi,11);
                g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, MAI, 
                                           RotationInit, ThinningInit,forFlag, 0.75);  
//          g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, MAI, RotationInit, ThinningInit,forFlag, 0.75);
     if (MAI > MAI_CountryUprotect[Country0-1])
     {   
         if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome))
                  {
                          managedForest.set(xi,yi,3.);
                          Rotation = rotMAI;
                          rotationType.set(xi,yi,1);
                    }else 
                    {     managedForest.set(xi,yi,2.);
                          Rotation = rotMaxBmTh;
                          rotationType.set(xi,yi,3);
                     }
       }else
       {                  managedForest.set(xi,yi,2.);
                          Rotation = rotMaxBmTh;
                          rotationType.set(xi,yi,3);
        }
       rotationForest.set(xi,yi,Rotation);	
       cohort.setRotPeriod(Rotation);
//       cohort.setStockingdegree(thinningForest.get(xi,yi));
//---grow forest for one year
//#include "growForest.cpp"
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
//          dima decision(1990
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
//		       , data["FTIMBER"] 
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v(1990)
//		       , wprice["regprice"]
//		       , wprice["regprice0"].v(2000)
//               , Rotation
//               , harvWood );
//          abBiomassO = cohort.getBm();
//	      double pDefIncome = abBiomassO * 
//                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	      //Immediate Pay if deforested (Slash and Burn)
//	      double sDefIncome = abBiomassO *
//		           (decision.priceTimber() * data["FTIMBER"].v()
//		         * (1. -coeff.HarvLoos.v()));
//	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		            + sDefIncome * data["SLASHBURN"].v();
// End of grow forest    
   		newHarvestTmp = harvWood * forestArea0;        
        harvestGrid.set(xi,yi,newHarvestTmp);   
        woodHarvest[Country0-1] += newHarvestTmp;
        woodLost[Country0-1] -= newHarvestTmp; 
 }   
 }else if (woodHarvest[Country0-1] >= (1.1 * woodHarvestStat[Country0-1])) 
 {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}
//if (Country0 == 61){cout<<"Inside IF!!??? >"<<endl;};
  if (thinningForest.get(xi,yi) > 0)
   {
                ThinningInit = -1.;
                thinningForest.set(xi,yi,-1.);
                RotationInit = biomasRot+1;
                rotationType.set(xi,yi,10);
                g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, MAI, 
                                           RotationInit, ThinningInit,forFlag, 0.75);  

      if (MAI > MAI_CountryUprotect[Country0-1])
       {   if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome))
           {   
                managedForest.set(xi,yi,0.);
//                 thinningForest.set(xi,yi,-1.);
                 Rotation = biomasRot+1;
                 rotationType.set(xi,yi,1);
             
            }else
            {     
                 managedForest.set(xi,yi,-1.);
                  Rotation = biomasRot+1;
                  rotationType.set(xi,yi,10);
              
            }
       }else 
       {  if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome))
            {          managedForest.set(xi,yi,-1.);
//                   thinningForest.set(xi,yi,-1.);
                   Rotation = biomasRot+1;
                   rotationType.set(xi,yi,10);
               
             }else
             {      
                    managedForest.set(xi,yi,-2.);
                    Rotation = biomasRot+1;
                    rotationType.set(xi,yi,10);
               
              }
        }
       rotationForest.set(xi,yi,Rotation);	
       cohort.setRotPeriod(Rotation);
//       cohort.setStockingdegree(thinningForest.get(xi,yi));
//---grow forest for one year
//#include "growForest.cpp"
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
//          dima decision(1990
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
//		       , data["FTIMBER"] 
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v(1990)
//		       , wprice["regprice"]
//		       , wprice["regprice0"].v(2000)
//               , Rotation
//               , harvWood );
//          abBiomassO = cohort.getBm();
//	      double pDefIncome = abBiomassO * 
//                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	      //Immediate Pay if deforested (Slash and Burn)
//	      double sDefIncome = abBiomassO *
//		           (decision.priceTimber() * data["FTIMBER"].v()
//		         * (1. -coeff.HarvLoos.v()));
//	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		            + sDefIncome * data["SLASHBURN"].v();
// End of grow forest     
   		newHarvestTmp = harvWood * forestArea0;        
        woodHarvest[Country0-1] -= (newHarvestTmp);
        harvestGrid.set(xi,yi,0);  
        woodLost[Country0-1] += newHarvestTmp; 
 }
}
//cout<< managedForest.get(xi,yi)<<endl;


/*  
cout <<"country=\t"<<Country0<<"\t woodHarvest=\t"<< woodHarvest[Country0-1];
cout <<"\t rotation=\t"<< Rotation<<"\t thinning= \t"<<thinningForest.get(xi,yi);
cout <<"\t sawnThW=\t"<<sawnThW*forestArea0*4<<"\t restThW=\t"<<restThW*forestArea0*4<<"\t";
cout <<"\t sawnW= \t"<<sawnW*forestArea0*4<<"\t restW= \t"<<restW*forestArea0*4;
cout << "\t sawnThWha=\t" << sawnThW << "\t restThWha=\t"<<restThW;
cout <<"\t sawnWha= \t"<<sawnW<<"\t restWha= \t"<<restW<<"\t abbiomassO= \t"<< abBiomassO<<"\t abbiomass= \t"<<data["BIOMASS"].v()<<"\t";
cout <<endl;
*/
//cout<< "managedForestSeting....\n";
  
    }//END for COUNTRY test
   } //End for PROTECT == 0
iter++;
  } // End for WHILE (cell loop) 

//cout << "First pass is finished"<< endl;
//    int optimUnmanaged = 0;
//    int optimMAI = 1;
//    int optimMaxBm = 2;
//    int optimMaxBmTh = 3;
//        biomasRot -> 10
//        biomasRotTh -> 11



//******************************************************************************
//**************************Forth Pass********************
//******************************************************************************
cout << "Start forth pass" << endl;

   iter = data_all.begin();

   while (iter != data_all.end())
   {

  int Country0 = 1;
  double HarvestTmp = 0;
  double newHarvestTmp = 0;
//  double forestArea0 = 0.;
  double abBiomassO = 0.;  	 
  int xi = 0;
  int yi = 0;
  double X = 0.;
  double Y = 0.;

	   if (iter->PROTECT[2000]==0)    // We consider only unprotected land
	   {

//cout<<Country0<<endl;       
//#include "countrySelect.cpp"
      if (  
            (iter->COUNTRY[2000] == 10) 
          || (iter->COUNTRY[2000] == 11)
          || (iter->COUNTRY[2000] == 17)          
          || (iter->COUNTRY[2000] == 25)
          || (iter->COUNTRY[2000] == 27)
          || (iter->COUNTRY[2000] == 30)
          || (iter->COUNTRY[2000] == 33)        
          || (iter->COUNTRY[2000] == 43)  
          || (iter->COUNTRY[2000] == 46)    
          || (iter->COUNTRY[2000] == 47)  
          || (iter->COUNTRY[2000] == 56)                        
          || (iter->COUNTRY[2000] == 61)
          || (iter->COUNTRY[2000] == 62)
          || (iter->COUNTRY[2000] == 69) 
          || (iter->COUNTRY[2000] == 71) 
          || (iter->COUNTRY[2000] == 82)    
          || (iter->COUNTRY[2000] == 83)   
          || (iter->COUNTRY[2000] == 89)  
          || (iter->COUNTRY[2000] == 92)    
          || (iter->COUNTRY[2000] == 96)   
          || (iter->COUNTRY[2000] == 107) 
          || (iter->COUNTRY[2000] == 112)    
          || (iter->COUNTRY[2000] == 113)
          || (iter->COUNTRY[2000] == 114)          
          || (iter->COUNTRY[2000] == 128)          
          || (iter->COUNTRY[2000] == 135)           
          || (iter->COUNTRY[2000] == 137)            
          || (iter->COUNTRY[2000] == 142)    
          || (iter->COUNTRY[2000] == 150)   
          || (iter->COUNTRY[2000] == 151) 
          || (iter->COUNTRY[2000] == 155)                                                                                                             
          || (iter->COUNTRY[2000] == 156)             
          || (iter->COUNTRY[2000] == 165) 
          || (iter->COUNTRY[2000] == 166)       
          || (iter->COUNTRY[2000] == 170)                         
          || (iter->COUNTRY[2000] == 179) 
          || (iter->COUNTRY[2000] == 180) 
          || (iter->COUNTRY[2000] == 190)        
          || (iter->COUNTRY[2000] == 194)    
          || (iter->COUNTRY[2000] == 196)                                               
          || (iter->COUNTRY[2000] == 197) ) { // Test only some countries

 do 
    {
        map<string, interpol> data;
        fillContainer(*iter,data);        

       int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
       int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
       double harvWood=0; //MG: harvestable wood, m3
  
    	            xi = (iter->x);
    	            yi = (iter->y);
//    	            X = (iter->x)*0.5+0.25-180;
//    	            Y = (iter->y)*0.5+0.25-90;
    	          
    	   //cout << "Xi = "<< xi <<"   Yi = "<< yi << endl;       
if (managedForest.get(xi,yi) > 0)
{ 
                       
       Country0 = (int)data["COUNTRY"].v();
       double LandAreaHa = data["LANDAREA"].v()*100;
       double forestArea0 = LandAreaHa * data["FOREST"].v();

//#include "DIMA_initialization.cpp"    
//#include "FM_initialization.cpp"
/////////////////////////////////// "FM_initialization
//MG: From Georg's "Optimal rotation time"
    //Get optimal rotation time
    //0 .. Rotation time Unmanaged forests       // Very long RotTime - untill the forest breaks
    //1 .. Highest average increment             // Short RotTime - approx = age of max MAI
    //2 .. Maximum avarage Biomass               // Very long RotTime - approx = age when forest accumulates max biomass
    //3 .. Maximum average Biomass with thinning // Very long (longer than 2) RotTime - approx = age when forest accumulates max biomass
    //4 .. Maximum harvest at final cut          // Long RotTime (approx 2 times longer than at 1 but almost 2 times shorter than at 2)
    //5 .. Maximum average harvest with final cut // Very short RotTime (shorter than at 1)

        int optimUnmanaged = 0;
        int optimMAI = 1;
        int optimMaxBm = 2;
        int optimMaxBmTh = 3;
        int optimHarvFin = 4;
        int optimHarvAve = 5;
        int rotUnmanaged = 0;
        int rotMAI = 0;
        int rotMaxBm = 0;
        int rotMaxBmTh = 0;
        int rotHarvFin = 0;
        int rotHarvAve = 0;

        g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
        g4m::ipol<double,double> hlv;  //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
        g4m::ipol<double,double> hle;  //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
        g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
        g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)

        sws.insert(10, .0);
        sws.insert(30, .6);
        hlv.insert(0, .0);
        hlv.insert(25, .7);
        hle.insert(0, .0);
        hle.insert(25, .7);
        dbv.insert(0, 2);
        dbe.insert(0, 3);
 
        g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
// End of FM initialisation



// double Thinning = 1.;
 int Rotation = 0;
 
  double forFlag = 0.; //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
//  if (data["FOREST"].v(1990) >0 && data["BIOMASS"].v() > 0 && maiForest.get(xi,yi) > 0)
  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  
  {
//          biomasRot = fi.gU(data["BIOMASS"].v(), maiForest.get(xi,yi), 1);       // rotation time to get current biomass (without thinning)
//          biomasRotTh = fi.gUt(data["BIOMASS"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
          biomasRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);       // rotation time to get current biomass (without thinning)
          biomasRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  

          if (thinningForest.get(xi,yi) == 1)
          {Rotation = biomasRotTh+1; 
          }else
          {Rotation = biomasRot+1;
          }
//          rotUnmanaged = fi.gTopt(MAI, optimUnmanaged);
          rotMAI = fi.gTopt(maiForest.get(xi,yi), optimMAI);
          rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), optimMaxBmTh);         
//          rotHarvAve = fi.gTopt(maiForest.get(xi,yi), optimHarvAve);
          forFlag = 1.0;
  }     

//cout << "\t Rotation1 = \t"<<Rotation<<endl;

			g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), 
			                            Rotation, thinningForest.get(xi,yi),forFlag, 0.75); 

cohort.aging();
            cohort.setRotPeriod(rotationForest.get(xi,yi));

//---grow forest for one year
//#include "growForest.cpp"
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
//          dima decision(1990
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
//		       , data["FTIMBER"] 
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v(1990)
//		       , wprice["regprice"]
//		       , wprice["regprice0"].v(2000)
//               , Rotation
//               , harvWood );
//          abBiomassO = cohort.getBm();
//	      double pDefIncome = abBiomassO * 
//                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	      //Immediate Pay if deforested (Slash and Burn)
//	      double sDefIncome = abBiomassO *
//		           (decision.priceTimber() * data["FTIMBER"].v()
//		         * (1. -coeff.HarvLoos.v()));
//	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		            + sDefIncome * data["SLASHBURN"].v();
// End of grow forest    

//harvestGrid.set(xi,yi,harvWood * forestArea0);   
//woodHarvest[Country0-1] += harvWood * forestArea0;


//#include "DIMA_initialization.cpp"    
//#include "FM_initialization.cpp"
//double decisionEstimate = (decission.forVal() * Hurdle_opt[Country0-1]) - 
//                            (decission.agrVal()+ defIncome);
//decisionGrid.set(xi,yi,decisionEstimate);

//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"countryB=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1];
//cout<<"\t 0.8harvestStat\t"<< 0.8 * woodHarvestStat[Country0-1]<<"\t 1.2harvestStat\t"<< 1.2 * woodHarvestStat[Country0-1];
//cout<<"\t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotationInit \t"<<Rotation<<endl;      

//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}
 if (woodHarvest[Country0-1] <= (0.9 * woodHarvestStat[Country0-1])) 
{
///if (Country0 == 61){cout<<"Inside IF (4th path) < !!???"<<endl;};
//
//		   if ((Rotation > rotMAI) && (rotationForest.get(xi,yi) > Rotation))
//		    { 
//	    		HarvestTmp = harvestGrid.get(xi,yi);
//                rotationForest.set(xi,yi,Rotation);		   
//                cohort.setRotPeriod(Rotation);
//#include "growForest.cpp"
//         		newHarvestTmp = harvWood * forestArea0;
//		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
////cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
////cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 0.8harvestStat\t"<< 0.8 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<Rotation<<endl;               
//                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
//
//            } else if (rotationForest.get(xi,yi) > rotMAI)
		   if (rotationForest.get(xi,yi) > rotMAI)
            { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,rotMAI);		   
                cohort.setRotPeriod(rotMAI);
//---grow forest for one year
//#include "growForest.cpp"
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
//          dima decision(1990
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
//		       , data["FTIMBER"] 
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v(1990)
//		       , wprice["regprice"]
//		       , wprice["regprice0"].v(2000)
//               , Rotation
//               , harvWood );
//          abBiomassO = cohort.getBm();
//	      double pDefIncome = abBiomassO * 
//                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	      //Immediate Pay if deforested (Slash and Burn)
//	      double sDefIncome = abBiomassO *
//		           (decision.priceTimber() * data["FTIMBER"].v()
//		         * (1. -coeff.HarvLoos.v()));
//	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		            + sDefIncome * data["SLASHBURN"].v();
// End of grow forest    
         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,newHarvestTmp);
//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 0.8harvestStat\t"<< 0.8 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<Rotation<<endl;               
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);

              }            
  

		  
 } else if (woodHarvest[Country0-1] >= 1.1 * woodHarvestStat[Country0-1]) 
 {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}
//if (Country0 == 61){cout<<"Inside IF (4th path) >!!???"<<endl;};        
//		   if (rotationForest.get(xi,yi) < Rotation)
//		    { 
//	    		HarvestTmp = harvestGrid.get(xi,yi);
//                rotationForest.set(xi,yi,Rotation);		   
//                cohort.setRotPeriod(Rotation);
//#include "growForest.cpp"
//         		newHarvestTmp = harvWood * forestArea0;
//		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
////cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";		     	
////cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 1.1harvestStat=\t"<< 1.1 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<Rotation<<endl;               
//                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
//             } else if (rotationForest.get(xi,yi) < rotMaxBmTh)
           if (rotationForest.get(xi,yi) < rotMaxBmTh)
			 {    HarvestTmp = harvestGrid.get(xi,yi);
			      rotationForest.set(xi,yi, rotMaxBmTh);
                  cohort.setRotPeriod(rotMaxBmTh);
//---grow forest for one year
//#include "growForest.cpp"
          res = cohort.aging();
          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
//          dima decision(1990
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
//		       , data["FTIMBER"] 
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v(1990)
//		       , wprice["regprice"]
//		       , wprice["regprice0"].v(2000)
//               , Rotation
//               , harvWood );
//          abBiomassO = cohort.getBm();
//	      double pDefIncome = abBiomassO * 
//                   (decision.priceTimber() * data["FTIMBER"].v() * (1. -coeff.HarvLoos.v()));
//	      //Immediate Pay if deforested (Slash and Burn)
//	      double sDefIncome = abBiomassO *
//		           (decision.priceTimber() * data["FTIMBER"].v()
//		         * (1. -coeff.HarvLoos.v()));
//	      defIncome = pDefIncome * (1. - data["SLASHBURN"].v())
//		            + sDefIncome * data["SLASHBURN"].v();
//// End of grow forest    
  			     newHarvestTmp = harvWood * forestArea0;
			     harvestGrid.set(xi,yi,newHarvestTmp);
//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";			     
//cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 1.1harvestStat=\t"<< 1.1 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<rotationForest.get(xi,yi)<<endl;               
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              }
 }
}

//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"countryE=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t harvestStat=\t"<< woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotationE \t"<<rotationForest.get(xi,yi);               
//cout <<"\t"<<sawnThW*forestArea0*4<<"\t"<<restThW*forestArea0*4<<"\t";
//cout <<sawnW*forestArea0*4<<"\t"<<restW*forestArea0*4;
//cout<<endl;

//cout <<rotationForest.get(xi,yi)<<"\t";
//cout <<sawnThW*forestArea0*4<<"\t"<<restThW*forestArea0*4<<"\t";
//cout <<sawnW*forestArea0*4<<"\t"<<restW*forestArea0*4;// <<"\t"<<abBiomassO*forestArea0<<"\t";
//cout << endl;

iter++;
//cout<<"iterC="<<((iter-1)->COUNTRY[2000])<<"\t iterC+1="<<((iter)->COUNTRY[2000]) <<endl;
} while (((iter-1)->COUNTRY[2000]) == ((iter)->COUNTRY[2000]));   // Check are we in the same country  // end for Within current country
//cout<<"iterS="<<((iter-1)->COUNTRY[2000])<<"\t iterCS+1="<<((iter)->COUNTRY[2000]) <<endl;

//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"countryE=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t harvestStat=\t"<< woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< data["BIOMASS"].v()<<"\t biomassO\t"<< abBiomassO<<"\t rotationE \t"<<rotationForest.get(xi,yi);               
//cout<<endl;

//cout <<X<<"\t"<<Y<<"\t"<<Country0<<"\t"<<year<<"\t"
//if (Country0 == 61) 
{cout<<"Country= "<<Country0<<"\t woodHarvestFinal= "<<woodHarvest[Country0-1]<<"\t woodHarvestStat= "<<woodHarvestStat[Country0-1] <<endl;}
} else{iter++;}  // end for Country select

} else{iter++;}  // end for IF unprotected
   

} //end for WHILE

	
//************************End of Forth Pass************************************
cout << "End of 4th pass"<<endl;



////************* TEST of INIT FM
//for (int i=1;i<=209;i++)
//{
//cout<<"Country= "<<i<<"\t managedCountInitFM = "<< managedCount[i-1]<<"\t woodHarvest = "<< woodHarvest[i-1]<<endl;
//}
//   iter = data_all.begin();
//
//   while (iter != data_all.end())
//   {
//
//  int Country0 = 1;
//  double HarvestTmp = 0;
//  double newHarvestTmp = 0;
////  double forestArea0 = 0.;
//  double abBiomassO = 0.;  	 
//  int xi = 0;
//  int yi = 0;
//  double X = 0.;
//  double Y = 0.;
//
//	   if (iter->PROTECT[2000]==0)    // We consider only unprotected land
//	   {
//
////cout<<Country0<<endl;       
////#include "countrySelect.cpp"
//      if (   (iter->COUNTRY[2000] == 11)
//          || (iter->COUNTRY[2000] == 25)
//          || (iter->COUNTRY[2000] == 33)        
//          || (iter->COUNTRY[2000] == 33)  
//          || (iter->COUNTRY[2000] == 61)
//          || (iter->COUNTRY[2000] == 62)
//          || (iter->COUNTRY[2000] == 69)             
//          || (iter->COUNTRY[2000] == 156)             
//          || (iter->COUNTRY[2000] == 165)             
//          || (iter->COUNTRY[2000] == 179)                     
//          || (iter->COUNTRY[2000] == 197) ) 
//        { // Test only some countries
//
//         map<string, interpol> data;
//        fillContainer(*iter,data);        
//
//       int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
//       int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
//       double harvWood=0; //MG: harvestable wood, m3
//
//
//    	            xi = (iter->x);
//    	            yi = (iter->y);
//    	            X = (iter->x)*0.5+0.25-180;
//    	            Y = (iter->y)*0.5+0.25-90;
//    	          
//    	   //cout << "Xi = "<< xi <<"   Yi = "<< yi << endl;       
//if (managedForest.get(xi,yi) > 0)
//{ 
//                       
//       Country0 = (int)data["COUNTRY"].v();
//       double LandAreaHa = data["LANDAREA"].v()*100;
//       double forestArea0 = LandAreaHa * data["FOREST"].v();
//
//
//        g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
//        g4m::ipol<double,double> hlv;  //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
//        g4m::ipol<double,double> hle;  //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
//        g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
//        g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)
//
//        sws.insert(10, .0);
//        sws.insert(30, .6);
//        hlv.insert(0, .0);
//        hlv.insert(25, .7);
//        hle.insert(0, .0);
//        hle.insert(25, .7);
//        dbv.insert(0, 2);
//        dbe.insert(0, 3);
// 
//        g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
//// End of FM initialisation
//
//
//
//// double Thinning = 1.;
// int Rotation = 0;
// 
//  double forFlag = 0.; //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
////  if (data["FOREST"].v(1990) >0 && data["BIOMASS"].v() > 0 && maiForest.get(xi,yi) > 0)
//  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
//  
//  {
//          biomasRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);       // rotation time to get current biomass (without thinning)
//          biomasRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
//
//          if (thinningForest.get(xi,yi) == 1)
//          {Rotation = biomasRotTh+1; 
//          }else
//          {Rotation = biomasRot+1;
//          }
//
//          forFlag = 1.0;
//  }     
//
////cout << "\t Rotation1 = \t"<<Rotation<<endl;
//
//			g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), 
//			                            Rotation, thinningForest.get(xi,yi),forFlag, 0.75); 
//
//cohort.aging();
//            cohort.setRotPeriod(rotationForest.get(xi,yi));
//
////---grow forest for one year
////#include "growForest.cpp"
//          res = cohort.aging();
//          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
//          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
//          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
//          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
//          harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v();
//          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//     CountriesWoodHarvestM3Year.inc(Country0,data["BYEAR"].v(),harvWood*forestArea0);
////     CountriesWoodHarvestM3Year.inc(Country0,data["BYEAR"].v(),harvWood); 
//}// end for if managed
//iter++;
//} else{iter++;}  // end for Country select
//
//} else{iter++;}  // end for IF unprotected
//	   
////int Nc=0;
//
//
//} //end for WHILE
///////////////////******   END of Test of init FM ************************

 } //END of void
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
 
void initLoop(int i, dataDetStruct &data_all, g4m::incrementTab &fi, ageStructVector &cohort_all, 
              ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
              griddata2<char> &thinningForest, griddata &rotationForest) 
 {
 
  double woodHarvest[209];
  double woodLost[209];
  double woodPot[209];

  for (int i=0; i<=208; i++){
    woodHarvest[i]=0.; 
    woodLost[i]=0.;
    woodPot[i]=0.;
  }
 
  int asID = 0;  // index in the ageStruct vector
  int kk=0;
  dataDetStruct::iterator iter = data_all.begin();
  while (iter != data_all.end()) {
    if (iter->PROTECT[2000]==0) {
      if (  
            (iter->COUNTRY[2000] == 10) 
          || (iter->COUNTRY[2000] == 11)
          || (iter->COUNTRY[2000] == 17)          
          || (iter->COUNTRY[2000] == 25)
          || (iter->COUNTRY[2000] == 27)
          || (iter->COUNTRY[2000] == 30)
          || (iter->COUNTRY[2000] == 33)        
          || (iter->COUNTRY[2000] == 43)  
          || (iter->COUNTRY[2000] == 46)    
          || (iter->COUNTRY[2000] == 47)  
          || (iter->COUNTRY[2000] == 56)                        
          || (iter->COUNTRY[2000] == 61)
          || (iter->COUNTRY[2000] == 62)
          || (iter->COUNTRY[2000] == 69) 
          || (iter->COUNTRY[2000] == 71) 
          || (iter->COUNTRY[2000] == 82)    
          || (iter->COUNTRY[2000] == 83)   
          || (iter->COUNTRY[2000] == 89)  
          || (iter->COUNTRY[2000] == 92)    
          || (iter->COUNTRY[2000] == 96)   
          || (iter->COUNTRY[2000] == 107) 
          || (iter->COUNTRY[2000] == 112)    
          || (iter->COUNTRY[2000] == 113)
          || (iter->COUNTRY[2000] == 114)          
          || (iter->COUNTRY[2000] == 128)          
          || (iter->COUNTRY[2000] == 135)           
          || (iter->COUNTRY[2000] == 137)            
          || (iter->COUNTRY[2000] == 142)    
          || (iter->COUNTRY[2000] == 150)   
          || (iter->COUNTRY[2000] == 151) 
          || (iter->COUNTRY[2000] == 155)                                                                                                             
          || (iter->COUNTRY[2000] == 156)             
          || (iter->COUNTRY[2000] == 165) 
          || (iter->COUNTRY[2000] == 166)       
          || (iter->COUNTRY[2000] == 170)                         
          || (iter->COUNTRY[2000] == 179) 
          || (iter->COUNTRY[2000] == 180) 
          || (iter->COUNTRY[2000] == 190)        
          || (iter->COUNTRY[2000] == 194)    
          || (iter->COUNTRY[2000] == 196)                                               
          || (iter->COUNTRY[2000] == 197) ) { // Test only some countries
        map<string, interpol> data;
        fillContainer(*iter,data);
        int xi = (iter->x);
        int yi = (iter->y);
//        double X = (iter->x)*0.5+0.25-180;
//        double Y = (iter->y)*0.5+0.25-90;
        int Country = (int)data["COUNTRY"].v();
        double LandAreaHa = data["LANDAREA"].v()*100;
        double forestArea0 = LandAreaHa * data["FOREST"].v();
        int aveRot = 0;                  // Average rotation time from statistics
        int PriceCi = PriceCiS[i];
        int iprice=i+1;
        coeff.PriceC.clear();
        coeff.PriceC.insert(0, PriceCi * data["CORRUPTION"].v());
        int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
        int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
        double harvWood=0.; //MG: harvestable wood, m3
        double harvWoodLost=0.; // wood lost in unmanaged forests, tC/ha
        double harvWoodNew=0.; //MG: harvestable wood, m3 (new forest)
        double harvWoodLostNew=0.; // wood lost in unmanaged forests, tC/ha (new forest)
//        dima decision((int)data["BYEAR"].v()
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
//		       , data["FTIMBER"]
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , data["GDP"]
//		       , coeff.HarvLoos
//		       , data["FOREST"].v((int)data["BYEAR"].v())
//		       , wprice["regprice"]
//		       , wprice["regprice0"].v()
//               , biomasRot+1
//               , harvWood );
        int beyear = (int)data["BYEAR"].v();
        double maxForInit = 1-data["BUILTUP"].v(beyear)-data["CROP"].v(beyear);
        if (maxForInit < 0) maxForInit = 0;
        double forestShare = data["FOREST"].v(); //Actual forest share
        if (forestShare > maxForInit) forestShare = maxForInit;
        double OforestShare = forestShare;
        double AforestShare = 0.;              //MG: Actual forest share
        double refForShare = forestShare;      //forest share of ref. year	
        double OfsNoPay = forestShare;         //MG: Forest share (deforested) without payment
        double AfsNoPay = 0.;                  //MG: Forest share (afforested) without payment
        double fsNoPay = OfsNoPay + AfsNoPay;  //MG: Forest share without payment
        g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
        g4m::ipol<double,double> hlv;   //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
        g4m::ipol<double,double> hle;   //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
        g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
        g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)
        sws.insert(10, .0);
        sws.insert(30, .6);
        hlv.insert(0, .0);
        hlv.insert(25, .7);
        hle.insert(0, .0);
        hle.insert(25, .7);
        dbv.insert(0, 2);
        dbe.insert(0, 3);
        double forFlag = 0.; //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
        double sawnWpot = 0.;
        double restWpot = 0.;
        double sawnThWpot = 0.;
        double restThWpot = 0.;
        int optimumType = 3;
        int MAIRot = 0;  //MG: optimal rotation time (see above)
        int optRotUnmanaged = 0;
        int rotationTimeCurr = 0; 
        int Rotation = 0;
        int rotMaxBm = 0;
        int rotMaxBmTh = 0;
        if (refForShare >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi)> 0) {
          biomasRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);
          biomasRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1); // with thinning
          MAIRot = fi.gTopt(maiForest.get(xi,yi), 1);
          forFlag = 1.0;
          rotMaxBm = fi.gTopt(maiForest.get(xi,yi), 2);                        
          rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), 3);
        }
        rotationTimeCurr = rotationForest.get(xi,yi);
        if (biomasRot < 0) biomasRot = 0;
        if (MAIRot < 0) MAIRot = 0;
        if (optRotUnmanaged < 0) optRotUnmanaged = 0;
        if (aveRot < 0) aveRot = 0;
        if (thinningForest.get(xi,yi) == 1) {
          Rotation = biomasRotTh+1; 
        } else {
          Rotation = biomasRot+1;
        }
// Existing (old forest)
        g4m::ageStruct *cohort = new g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
        g4m::ageStruct::v res;    // MG: results vector for the set (old) forest
        res = cohort->aging();
//        res = cohort.aging();        
        sawnWpot = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
        restWpot = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
        sawnThWpot = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
        restThWpot = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        double potHarvest = (sawnWpot + restWpot + sawnThWpot + restThWpot) * data["FTIMBER"].v(); // Potential sustainable harvest, m3/year
//        cohort->setRotPeriod(rotationTimeCurr);

////********Testing harvest
//        res = cohort->aging();
//  double sawnW = 0.;
//  double restW = 0.;
//  double sawnThW = 0.;
//  double restThW = 0.;
//  double sawnWlost = 0.;
//  double restWlost = 0.;
//  double sawnThWlost = 0.;
//  double restThWlost = 0.;
//
//    if (thinningForest.get(xi,yi)>0) {  
//    sawnW = res.enSw ;           // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
//    restW = res.enRw ;           // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
//    sawnThW = res.vnSw ;         // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
//    restThW = res.vnRw ;         // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
//    woodPot[Country-1]+=potHarvest * forestArea0;
//  } else {
//    sawnWlost = res.enSw ;           // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
//    restWlost = res.enRw ;           // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
//    sawnThWlost = res.vnSw ;         // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
//    restThWlost = res.vnRw ;         // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
//  }
//  harvWood = (sawnW + restW + sawnThW + restThW) * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
//  harvWoodLost = (sawnWlost + restWlost + sawnThWlost + restThWlost); // Total current "lost" wood in the cell, tC (in remote forests)
//  woodHarvest[Country-1]+=harvWood * forestArea0;
//  woodLost[Country-1]+=harvWoodLost * forestArea0;
//  woodPot[Country-1]+=potHarvest * forestArea0;
////********************************************************

// New forest (planted/afforested)
        g4m::ageStruct *newCohort = new g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),0, 0.75); 
        g4m::ageStruct::v newRes; // MG: results vector for the new (planted/afforested) forest
        dat singleCell;
        singleCell.Rotation = Rotation;
        singleCell.LandAreaHa = LandAreaHa;
        singleCell.potHarvest = potHarvest;
        singleCell.forestShare = forestShare;
        singleCell.OforestShare = OforestShare;
        singleCell.AforestShare = AforestShare;
        singleCell.prevOForShare = refForShare; //MG: Old Forest share in the previous reporting year
        singleCell.prevOForShareRP = OforestShare; //MG: Old Forest share in the previous reporting year
        singleCell.prevAForShareRP = 0.; //MG: New (afforested) Forest share in the previous reporting year
        singleCell.AforestSharePrev = 0.;
	    singleCell.savedCarbonPrev = 0.;
        singleCell.gainedCarbonPrev = 0.;
        singleCell.EmissionsTotPrev = 0.;
        singleCell.EmissionsAfforPrev = 0.;
        singleCell.prevPlantPhytHaBmGr = 0.;
        singleCell.prevPlantPhytHaBlGr = 0.;
        singleCell.deforestHaTot=0.;
        singleCell.afforestHaTot=0.;
        singleCell.EmissionsProduct= 0.;  
        singleCell.EmissionsLitter = 0.;  
        singleCell.EmissionsSOC = 0.;      
        singleCell.EmissionsSlashBurn = 0.;
        singleCell.EmissionsDeadBurn = 0.;
        singleCell.EmissionsCRootBurn = 0.;    
        singleCell.EmissionsTot = 0.;     
        singleCell.EmLitterAffor =0.;
        singleCell.EmSOCAffor = 0.; 
        singleCell.EmissionsAffor = 0.;
        for (int i=0; i<110; ++i) {
          singleCell.forestAgeShare[i] = 0.;
          singleCell.BDeadA[i]=0.;
          singleCell.LitterA[i]=0.;
          singleCell.SOCA[i]=0;
          singleCell.ProdLongA[i]=0.;
          singleCell.ProdShortA[i]=0.;
          singleCell.deforestA[i]=0;
          singleCell.FineRootA[i]=0;
          singleCell.LitterAffor[i]=0.;
          singleCell.SOCaffor[i]=0.;
        }
        singleCell.prevReportYear = data["BYEAR"].v();
        singleCell.ireportYear=0;
// saving results to initial vectors
        cohort_all.push_back(*cohort);
        newCohort_all.push_back(*newCohort);
        dat_all.push_back(singleCell);
        iter->asID = asID;
        asID++;

//if (thinningForest.get(xi,yi)>0){
//     CountriesManagedForHa.inc(Country,data["BYEAR"].v(),OforestShare*singleCell.LandAreaHa);
//     CountriesManagedCount.inc(Country,data["BYEAR"].v(),1);
//     }
//          
////     CountriesWoodHarvestM3Year.inc(Country,data["BYEAR"].v(),(harvWood*OforestShare+harvWoodNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa);
//     CountriesWoodLoosCYear.inc(Country,data["BYEAR"].v(),(harvWoodLost*OforestShare+harvWoodLostNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa);

if (asID % 1000 == 0) cout << "asID = " << asID << endl;
      }   // End for Country select
    }     // End for IF unprotected 
    iter++;
kk++;
//if (kk % 100 == 0)
//    cout << "kk = " << kk << endl;
  }       // End of WHILE cells loop
//for (int i=1;i<=209;i++){
//cout<<"Country= "<<i<<"\t harvWoodM3 = "<< woodHarvest[i-1]<<"\t woodPotM3 = "<< woodPot[i-1]<<"\t harvLostC = "<< woodLost[i-1]<<endl;
//}
 }
