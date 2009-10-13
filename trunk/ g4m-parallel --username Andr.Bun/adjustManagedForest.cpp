 void adjustManagedForest(dataDetStruct &data_all, g4m::incrementTab &fi, ageStructVector &cohort_all, 
              ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
              griddata2<char> &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
              griddata &rotationForestNew, griddata2<char> &thinningForestNew,
              griddata2<char> &rotationType, griddata &harvestGrid, int year) 
 
 
 {


  double woodHarvest[27];
//  double woodLost[27];
//  int managedCount[27];

  for (int i=0; i<=26; i++){
    woodHarvest[i]=CountriesWoodHarvestM3Year.getRegionSum(i, year); 
//    woodLost[i]=0.;
//    woodHarvestStat[i]=0.;
//    managedCount[i]=0;
  }
 
//------------------------------------------------------------------------------
// 
// -------Zero Adjust thinning if population density changed --------------------
//
//------------------------------------------------------------------------------

if ((year > 2000) & (year/10 % 10 == 0))
{
  dataDetStruct::iterator iter = data_all.begin();
  iter = data_all.begin();

//cout << "Putting data for current cell into conteiner... "<< endl;
   while (iter != data_all.end())
   {
	if (iter->PROTECT[2000] == 0)
	{

    int asID = iter->asID;
    int region = iter->POLESREG[2000];
    char regionch[2];
    int2str(region,regionch);
    int xi = (iter->x);
    int yi = (iter->y);
    
  if (woodHarvest[region-1] <= 0.9 * wprod[regionch].v(year))
    {if (managedForest.get(xi,yi)<=0)
      {map<string, interpol> data;
      fillContainer(*iter,data);
      if ((data["POPDENS"].v(1990) >0) && (data["GDP"].v(1990) > 0))
      {
      double sawnW = 0.;
      double restW = 0.;
      double sawnThW = 0.;
      double restThW = 0.;
      double sawnWnew = 0.;
      double restWnew = 0.;
      double sawnThWnew = 0.;
      double restThWnew = 0.;
      double harvestTmp = 0.;
      double newHarvestTmp = 0.;
      int biomassRot = 0;
      int rotMAI = 0;
      int rotMaxBmTh = 0;
      int Rotation = 0;


    g4m::ageStruct cohortTmp = cohort_all[asID];
    cout << "addr1 " << &cohortTmp << "   addr2 " << &(cohort_all[asID]) << endl;
    system("pause");
   	harvestTmp = harvestGrid.get(xi,yi);              

//cout << "xi= "<<xi<<"\t yi= "<<yi<<"\t MAI= "<< MAI<<"\t mai= "<< maiForest.get(xi,yi)<<"\t NPP= "<<iter->NPP[2000]<<"\t dNPP= "<<data["NPP"].v()<<endl; 

  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  {
          biomassRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
          rotMAI = fi.gTopt(maiForest.get(xi,yi), 1);
          rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), 3);         
  }     

//---------------------------------------------------
           int newForAge = newCohort_all[asID].getActiveAge();
           if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
            {
              g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
                         
             if (managedForest.get(xi,yi) == 0)
                  {
                         managedForest.set(xi,yi,3);
                         Rotation = rotMAI;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,2);
                         Rotation = rotMaxBmTh;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                
                rotationForestNew.set(xi,yi,Rotation);	
                thinningForestNew.set(xi,yi,1.);	
                cohortTmpNew.setRotPeriod(rotMaxBmTh);
                cohortTmpNew.setStockingdegree(1.);        
                newCohort_all[asID].setRotPeriod(rotMaxBmTh);
                newCohort_all[asID].setStockingdegree(1.);  
                g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
                sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
                newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                             (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                              data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3


           }else // New forest age <= rotation -> don't change FM for the new forest
           {
                      if (managedForest.get(xi,yi) == 0)
                  {
                         managedForest.set(xi,yi,3);
                         Rotation = rotMAI;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,2);
                         Rotation = rotMaxBmTh;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                 
                 newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                                   data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
            } // End    else // New forest age < rotation -> don't change FM for the new forest
              	harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        }  //End  if ((data["POPDENS"].v(1990) >0) && (data["GDP"].v(1990) > 0))
      } // End  if (managedForest(xi,yi)<=0)

 }else if (woodHarvest[region-1] >= 1.1 * wprod[regionch].v(year)) 
 {
    if (managedForest.get(xi,yi)>0) 
        {map<string, interpol> data;
        fillContainer(*iter,data);
     if ((data["POPDENS"].v(1990) == 0) && (data["GDP"].v(1990) == 0))
     {

      double sawnW = 0.;
      double restW = 0.;
      double sawnThW = 0.;
      double restThW = 0.;
      double sawnWnew = 0.;
      double restWnew = 0.;
      double sawnThWnew = 0.;
      double restThWnew = 0.;
//      double harvWoodTmp = 0.;
      double harvestTmp = 0.;
      double newHarvestTmp = 0.;
      int biomassRotTh = 0;
      int rotMAI = 0;
      int rotMaxBm = 0;
      int Rotation = 0;


    g4m::ageStruct cohortTmp = cohort_all[asID];
   	harvestTmp = harvestGrid.get(xi,yi);              

//cout << "xi= "<<xi<<"\t yi= "<<yi<<"\t MAI= "<< MAI<<"\t mai= "<< maiForest.get(xi,yi)<<"\t NPP= "<<iter->NPP[2000]<<"\t dNPP= "<<data["NPP"].v()<<endl; 

  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  {
          biomassRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (without thinning)  
          rotMAI = fi.gTopt(maiForest.get(xi,yi), 1);
          rotMaxBm = fi.gTopt(maiForest.get(xi,yi), 2);  
  }     


//---------------------------------------------------
           int newForAge = newCohort_all[asID].getActiveAge();
           if (newForAge > biomassRotTh)  // New forest age > rotation -> change FM for the new forest
            {
              g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
                         
             if (managedForest.get(xi,yi) == 2)
                  {
                         managedForest.set(xi,yi,-1);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,-2);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,-1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(-1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(-1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                
                rotationForestNew.set(xi,yi,Rotation);	
                thinningForestNew.set(xi,yi,-1.);	
                cohortTmpNew.setRotPeriod(rotMaxBm);
                cohortTmpNew.setStockingdegree(-1.);        
                newCohort_all[asID].setRotPeriod(rotMaxBm);
                newCohort_all[asID].setStockingdegree(-1.);  
                g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
                sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
                newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                             (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                              data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3


           }else // New forest age <= rotation -> don't change FM for the new forest
           {
             if (managedForest.get(xi,yi) == 2)
                  {
                         managedForest.set(xi,yi,-1);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,-2);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,-1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(-1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(-1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                 
                 newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                                   data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
            } // End    else // New forest age < rotation -> don't change FM for the new forest
              	harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
          } // End  if ((data["POPDENS"].v(1990) == 0) && (data["GDP"].v(1990) == 0))
        } // End  if ((managedForest(xi,yi)>0))
 } // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].v(year))   


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
   } //End for PROTECT == 0
iter++;
  } // End for WHILE (cell loop) 
} // End if year

//cout << "Zero pass is finished"<< endl; 
 
 
 
  
//----First pass = adjust rotation time -------  
dataDetStruct::iterator iter = data_all.begin();
while (iter != data_all.end()) {
 if (iter->PROTECT[2000]==0) {
  do 
  {
    int asID = iter->asID;
    int region = iter->POLESREG[2000];
    char regionch[2];
    int2str(region,regionch);
    int xi = (iter->x);
    int yi = (iter->y);
    
  double sawnW = 0.;
  double restW = 0.;
  double sawnThW = 0.;
  double restThW = 0.;
  double sawnWnew = 0.;
  double restWnew = 0.;
  double sawnThWnew = 0.;
  double restThWnew = 0.;
//  double harvWoodTmp = 0.;
  int biomassRotTh = 0;
  int rotMAI = 0;
  int rotMaxBmTh = 0;
  double harvestTmp = 0;
  double newHarvestTmp = 0;

    map<string, interpol> data;
    fillContainer(*iter,data);
  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  {
          biomassRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
          rotMAI = fi.gTopt(maiForest.get(xi,yi), 1);
          rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), 3);         
  }     

  if (woodHarvest[region-1] <= 0.9 * wprod[regionch].v(year))
    {if (managedForest.get(xi,yi)>=2)
     {

//cohortTmp = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
//cohortTmpNew = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
//      
//       for (int age=0;age<=200;age++){
//        cohortTmp.setArea(age,cohort_all[asID].getArea(age));
//        cohortTmp.setBm(age,cohort_all[asID].getBm(age));
//        cohortTmpNew.setArea(age,newCohort_all[asID].getArea(age));
//        cohortTmpNew.setBm(age,newCohort_all[asID].getBm(age));
//        }                
          g4m::ageStruct cohortTmp = cohort_all[asID];
          int newForAge = newCohort_all[asID].getActiveAge();         
    	harvestTmp = harvestGrid.get(xi,yi);
       if (rotMAI < rotationForest.get(xi,yi))
       {
        rotationForest.set(xi,yi,rotMAI);		   
        cohortTmp.setRotPeriod(rotMAI);
        cohort_all[asID].setRotPeriod(rotMAI);
        g4m::ageStruct::v resTmp = cohortTmp.aging();
        sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
        restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
        sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
        restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
           if ((newForAge > biomassRotTh) & (rotMAI < biomassRotTh)) {
             g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
             cohortTmpNew.setRotPeriod(rotMAI);
             newCohort_all[asID].setRotPeriod(rotMAI);    
             rotationForestNew.set(xi,yi,rotMAI);
             g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
             sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
             restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
             sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
             restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
             newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                         (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                          data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
           }else{
             newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                               data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
           }
      	harvestGrid.set(xi,yi,newHarvestTmp);
        woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        }

       }  // end for if (managedForest(xi,yi)>=2)          
    
    } //end for if (woodHarvest[region] <= 0.9 * wprod[region].v(year))
    else if (woodHarvest[region-1] >= 1.1 * wprod[regionch].v(year)) 
    {if (managedForest.get(xi,yi)<=2)
     {
    map<string, interpol> data;
    fillContainer(*iter,data);
//cohortTmp = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
//cohortTmpNew = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
      
//       for (int age=0;age<=200;age++){
//        cohortTmp.setArea(age,cohort_all[asID].getArea(age));
//        cohortTmp.setBm(age,cohort_all[asID].getBm(age));
//        cohortTmpNew.setArea(age,newCohort_all[asID].getArea(age));
//        cohortTmpNew.setBm(age,newCohort_all[asID].getBm(age));
//        }                
          g4m::ageStruct cohortTmp = cohort_all[asID];
          g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
          int newForAge = newCohort_all[asID].getActiveAge();  

    	harvestTmp = harvestGrid.get(xi,yi);
       if (rotMaxBmTh > rotationForest.get(xi,yi))
       {
        rotationForest.set(xi,yi,rotMaxBmTh);		   
        cohortTmp.setRotPeriod(rotMaxBmTh);
        cohort_all[asID].setRotPeriod(rotMaxBmTh);       
        g4m::ageStruct::v resTmp = cohortTmp.aging();
        sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
        restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
        sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
        restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
           if ((newForAge > biomassRotTh) & (rotMaxBmTh > biomassRotTh)) {
             g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
             cohortTmpNew.setRotPeriod(rotMaxBmTh);
             newCohort_all[asID].setRotPeriod(rotMaxBmTh);    
             rotationForestNew.set(xi,yi,rotMAI);
             g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
             sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
             restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
             sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
             restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
             newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                         (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                          data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
           }else{
             newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                               data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
           }
      	harvestGrid.set(xi,yi,newHarvestTmp);
        woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        }
        
       }  // end for if (managedForest(xi,yi)<=2)       
      }   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].v(year))
        iter++;
     } while (((iter-1)->POLESREG[2000]) == ((iter)->POLESREG[2000]));   // Check are we in the same region  // end for Within current country
        iter++;
  } //End Protect
} // End While
// ----- End of First pass
//------------------------------------------------------------------------------
// 
// -------Second Adjust thinning -----------------------------------------------
//
//------------------------------------------------------------------------------
  iter = data_all.begin();

//cout << "Putting data for current cell into conteiner... "<< endl;
   while (iter != data_all.end())
   {
	if (iter->PROTECT[2000] == 0)
	{

    int asID = iter->asID;
    int region = iter->POLESREG[2000];
    char regionch[2];
    int2str(region,regionch);
    int xi = (iter->x);
    int yi = (iter->y);
    
  if (woodHarvest[region-1] <= 0.9 * wprod[regionch].v(year))
    {if ((managedForest.get(xi,yi)<=0) & (managedForest.get(xi,yi)>-2))
     {
    map<string, interpol> data;
    fillContainer(*iter,data);

      double sawnW = 0.;
      double restW = 0.;
      double sawnThW = 0.;
      double restThW = 0.;
      double sawnWnew = 0.;
      double restWnew = 0.;
      double sawnThWnew = 0.;
      double restThWnew = 0.;
//      double harvWoodTmp = 0.;
      double harvestTmp = 0.;
      double newHarvestTmp = 0.;
      int biomassRot = 0;
      int rotMAI = 0;
      int rotMaxBmTh = 0;
      int Rotation = 0;
      


    g4m::ageStruct cohortTmp = cohort_all[asID];
   	harvestTmp = harvestGrid.get(xi,yi);              

//cout << "xi= "<<xi<<"\t yi= "<<yi<<"\t MAI= "<< MAI<<"\t mai= "<< maiForest.get(xi,yi)<<"\t NPP= "<<iter->NPP[2000]<<"\t dNPP= "<<data["NPP"].v()<<endl; 

  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  {
          biomassRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
          rotMAI = fi.gTopt(maiForest.get(xi,yi), 1);
          rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), 3);         
  }     


//---------------------------------------------------
           int newForAge = newCohort_all[asID].getActiveAge();  
           if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
            {
              g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
                         
             if (managedForest.get(xi,yi) == 0)
                  {
                         managedForest.set(xi,yi,3);
                         Rotation = rotMAI;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,2);
                         Rotation = rotMaxBmTh;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                
                rotationForestNew.set(xi,yi,Rotation);	
                thinningForestNew.set(xi,yi,1.);	
                cohortTmpNew.setRotPeriod(rotMaxBmTh);
                cohortTmpNew.setStockingdegree(1.);        
                newCohort_all[asID].setRotPeriod(rotMaxBmTh);
                newCohort_all[asID].setStockingdegree(1.);  
                g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
                sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
                newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                             (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                              data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3


           }else // New forest age <= rotation -> don't change FM for the new forest
           {
                      if (managedForest.get(xi,yi) == 0)
                  {
                         managedForest.set(xi,yi,3);
                         Rotation = rotMAI;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,2);
                         Rotation = rotMaxBmTh;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                 
                 newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                                   data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
            } // End    else // New forest age < rotation -> don't change FM for the new forest
              	harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        } // End  if ((managedForest(xi,yi)<=0) & (managedForest(xi,yi)>-2)
   
 }else if (woodHarvest[region-1] >= 1.1 * wprod[regionch].v(year)) 
 {
    if ((managedForest.get(xi,yi)>0) & (managedForest.get(xi,yi)<3))
     {
    map<string, interpol> data;
    fillContainer(*iter,data);

      double sawnW = 0.;
      double restW = 0.;
      double sawnThW = 0.;
      double restThW = 0.;
      double sawnWnew = 0.;
      double restWnew = 0.;
      double sawnThWnew = 0.;
      double restThWnew = 0.;
//      double harvWoodTmp = 0.;
      double harvestTmp = 0.;      
      double newHarvestTmp = 0;
      int biomassRot = 0;
      int rotMAI = 0;
      int rotMaxBm = 0;
      int Rotation = 0;


    g4m::ageStruct cohortTmp = cohort_all[asID];
   	harvestTmp = harvestGrid.get(xi,yi);              

//cout << "xi= "<<xi<<"\t yi= "<<yi<<"\t MAI= "<< MAI<<"\t mai= "<< maiForest.get(xi,yi)<<"\t NPP= "<<iter->NPP[2000]<<"\t dNPP= "<<data["NPP"].v()<<endl; 

  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  {
          biomassRot = fi.gU(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (without thinning)  
          rotMaxBm = fi.gTopt(maiForest.get(xi,yi), 2);  
  }     


//---------------------------------------------------
           int newForAge = newCohort_all[asID].getActiveAge();  
           if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
            {
              g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
                         
             if (managedForest.get(xi,yi) == 2)
                  {
                         managedForest.set(xi,yi,-1);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,-2);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,-1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(-1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(-1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                
                rotationForestNew.set(xi,yi,Rotation);	
                thinningForestNew.set(xi,yi,-1.);	
                cohortTmpNew.setRotPeriod(Rotation);
                cohortTmpNew.setStockingdegree(-1.);        
                newCohort_all[asID].setRotPeriod(Rotation);
                newCohort_all[asID].setStockingdegree(-1.);  
                g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
                sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
                newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                             (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                              data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3


           }else // New forest age <= rotation -> don't change FM for the new forest
           {
             if (managedForest.get(xi,yi) == 2)
                  {
                         managedForest.set(xi,yi,-1);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,1);
                   }else 
                   {     managedForest.set(xi,yi,-2);
                         Rotation = rotMaxBm;
                         rotationType.set(xi,yi,3);
                    }

                rotationForest.set(xi,yi,Rotation);	
                thinningForest.set(xi,yi,-1.);	
                cohortTmp.setRotPeriod(Rotation);
                cohortTmp.setStockingdegree(-1.);
                cohort_all[asID].setRotPeriod(Rotation);
                cohort_all[asID].setStockingdegree(-1.);
                g4m::ageStruct::v resTmp = cohortTmp.aging();
                sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        
                 
                 newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                                   data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
            } // End    else // New forest age < rotation -> don't change FM for the new forest
              	harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        } // End  if ((managedForest(xi,yi)>0) & (managedForest(xi,yi)<3)
 } // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].v(year))   


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
   } //End for PROTECT == 0
iter++;
  } // End for WHILE (cell loop) 

//cout << "Second pass is finished"<< endl;

//******************************************************************************
//**************************Third Pass********************
//******************************************************************************
cout << "Start third pass" << endl;
iter = data_all.begin();
while (iter != data_all.end()) {
 if (iter->PROTECT[2000]==0) {
  do 
  {
    int asID = iter->asID;
    int region = iter->POLESREG[2000];
    char regionch[2];
    int2str(region,regionch);
    int xi = (iter->x);
    int yi = (iter->y);
    
   
  double sawnW = 0.;
  double restW = 0.;
  double sawnThW = 0.;
  double restThW = 0.;
  double sawnWnew = 0.;
  double restWnew = 0.;
  double sawnThWnew = 0.;
  double restThWnew = 0.;
  double harvWoodTmp = 0.;
  int biomassRotTh = 0;
  int rotMAI = 0;
  int rotMaxBmTh = 0;
  double harvestTmp = 0;
  double newHarvestTmp = 0;

    map<string, interpol> data;
    fillContainer(*iter,data);

  if (data["FOREST"].v(1990) >0 && data["CABOVEHA"].v() > 0 && maiForest.get(xi,yi) > 0)  
  {
          biomassRotTh = fi.gUt(data["CABOVEHA"].v(), maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
          rotMAI = fi.gTopt(maiForest.get(xi,yi), 1);
          rotMaxBmTh = fi.gTopt(maiForest.get(xi,yi), 3);         
  }     

  if (woodHarvest[region-1] <= 0.9 * wprod[regionch].v(year))
    {if (managedForest.get(xi,yi)>=1)
     {
//cohortTmp = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
//cohortTmpNew = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
//      
//       for (int age=0;age<=200;age++){
//        cohortTmp.setArea(age,cohort_all[asID].getArea(age));
//        cohortTmp.setBm(age,cohort_all[asID].getBm(age));
//        cohortTmpNew.setArea(age,newCohort_all[asID].getArea(age));
//        cohortTmpNew.setBm(age,newCohort_all[asID].getBm(age));
//        }                
          g4m::ageStruct cohortTmp = cohort_all[asID];
          int newForAge = newCohort_all[asID].getActiveAge();          
    	harvestTmp = harvestGrid.get(xi,yi);
       if (rotMAI < rotationForest.get(xi,yi))
       {
        rotationForest.set(xi,yi,rotMAI);		   
        cohortTmp.setRotPeriod(rotMAI);
        cohort_all[asID].setRotPeriod(rotMAI);
        g4m::ageStruct::v resTmp = cohortTmp.aging();
        sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
        restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
        sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
        restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
           if ((newForAge > biomassRotTh) & (rotMAI < biomassRotTh)) {
             g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
             cohortTmpNew.setRotPeriod(rotMAI);
             newCohort_all[asID].setRotPeriod(rotMAI);    
             rotationForestNew.set(xi,yi,rotMAI);
             g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
             sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
             restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
             sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
             restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
             newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                         (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                          data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
           }else{
             newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                               data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
           }
      	harvestGrid.set(xi,yi,newHarvestTmp);
        woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        }

       }  // end for if (managedForest(xi,yi)>=2)          
    
    } //end for if (woodHarvest[region] <= 0.9 * wprod[region].v(year))
    else if (woodHarvest[region-1] >= 1.1 * wprod[regionch].v(year)) 
    {if ((managedForest.get(xi,yi)>0) & (managedForest.get(xi,yi)<=3))
     {
    map<string, interpol> data;
    fillContainer(*iter,data);
//cohortTmp = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
//cohortTmpNew = g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);
      
//       for (int age=0;age<=200;age++){
//        cohortTmp.setArea(age,cohort_all[asID].getArea(age));
//        cohortTmp.setBm(age,cohort_all[asID].getBm(age));
//        cohortTmpNew.setArea(age,newCohort_all[asID].getArea(age));
//        cohortTmpNew.setBm(age,newCohort_all[asID].getBm(age));
//        }                
          g4m::ageStruct cohortTmp = cohort_all[asID];
          g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
          int newForAge = newCohort_all[asID].getActiveAge();  
    	harvestTmp = harvestGrid.get(xi,yi);
       if (rotMaxBmTh > rotationForest.get(xi,yi))
       {
        rotationForest.set(xi,yi,rotMaxBmTh);		   
        cohortTmp.setRotPeriod(rotMaxBmTh);
        cohort_all[asID].setRotPeriod(rotMaxBmTh);       
        g4m::ageStruct::v resTmp = cohortTmp.aging();
        sawnW = resTmp.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
        restW = resTmp.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
        sawnThW = resTmp.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
        restThW = resTmp.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
           if ((newForAge > biomassRotTh) & (rotMaxBmTh > biomassRotTh)) {
                 g4m::ageStruct cohortTmpNew = newCohort_all[asID]; 
             cohortTmpNew.setRotPeriod(rotMaxBmTh);
             newCohort_all[asID].setRotPeriod(rotMaxBmTh);    
             rotationForestNew.set(xi,yi,rotMAI);
             g4m::ageStruct::v resTmpNew = cohortTmpNew.aging();
             sawnWnew = resTmpNew.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
             restWnew = resTmpNew.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
             sawnThWnew = resTmpNew.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
             restThWnew = resTmpNew.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
             newHarvestTmp = ((sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare +
                         (sawnWnew + restWnew + sawnThWnew + restThWnew) * dat_all[asID].AforestShare) * 
                          data["LANDAREA"].v() * 100 * data["FTIMBER"].v(); // Total current harvested wood in the cell, m3
           }else{
             newHarvestTmp = (sawnW + restW + sawnThW + restThW) * dat_all[asID].OforestShare * 
                               data["LANDAREA"].v() * 100 * data["FTIMBER"].v();
           }
      	harvestGrid.set(xi,yi,newHarvestTmp);
        woodHarvest[region-1] += (newHarvestTmp-harvestTmp);
        }
        
       }  // end for if (managedForest(xi,yi)<=2)       
      }   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].v(year))
        iter++;
     } while (((iter-1)->POLESREG[2000]) == ((iter)->POLESREG[2000]));   // Check are we in the same region  // end for Within current country
        iter++;
  } // End protect
} // End While
// 
//************************End of Third Pass************************************
cout << "End of Third pass"<<endl;
}
