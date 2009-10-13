#ifndef FILLCONTAINER_CPP
#define FILLCONTAINER_CPP

void fillContainer(g4m::dataStruct &, map<string, interpol> &);

void fillContainer(g4m::dataStruct &plot, map<string, interpol> &data)
 {
  data.erase("FOREST");
  if ((plot.FOREST[2000])+(plot.CROP[2000])+(plot.BUILTUP[2000])>1)
    {plot.FOREST.insert(2000,(1-((plot.CROP[2000])+(plot.BUILTUP[2000]))));}
  data.erase("COUNTRY");
  data["COUNTRY"].insert(2000, plot.COUNTRY[2000]);
  data.erase("POTVEG");
  data["POTVEG"].insert(2000,plot.POTVEG[2000]);
  data.erase("PROTECT");            
  data["PROTECT"].insert(2000,plot.PROTECT[2000]);
  data.erase("USED");            
  data["USED"].insert(2000,plot.USED[2000]);
  data.erase("LANDAREA");            
  data["LANDAREA"].insert(2000,plot.LANDAREA[2000]);
  data.erase("NPP");            
  data["NPP"].insert(2000,plot.NPP[2000]);
  data.erase("BIOMASS");
  data["BIOMASS"].insert(2000,plot.BIOMASS[2000]);
  data.erase("BIOMASSBL");            
  data["BIOMASSBL"].insert(2000,plot.BIOMASSBL[2000]);     
  data.erase("CABOVEHA");                   
  data["CABOVEHA"].insert(2000,plot.CABOVEHA[2000]);  
  data.erase("CBELOWHA");            
  data["CBELOWHA"].insert(2000,plot.CBELOWHA[2000]);
  data.erase("CDEADHA");            
  data["CDEADHA"].insert(2000,plot.CDEADHA[2000]);  
  data.erase("CLITTERHA");            
  data["CLITTERHA"].insert(2000,plot.CLITTERHA[2000]);     
  data.erase("SOCHA");                   
  data["SOCHA"].insert(2000,plot.SOCHA[2000]);              
  data.erase("FOREST");
  data["FOREST"].insert(2000,plot.FOREST[2000]);  
  data.erase("R");            
  data["R"].insert(2000,plot.R[2000]);
  data.erase("FRACLONGPROD");            
  data["FRACLONGPROD"].insert(2000,plot.FRACLONGPROD[2000]);   
  data.erase("CORRUPTION");                     
  data["CORRUPTION"].insert(2000,plot.CORRUPTION[2000]);  
  data.erase("SLASHBURN");            
  data["SLASHBURN"].insert(2000,plot.SLASHBURN[2000]); 
  data.erase("IIASA_REGION");                                   
  data["IIASA_REGION"].insert(2000,plot.IIASA_REGION[2000]);      
  data.erase("SAGRSUIT");                                   
  data["SAGRSUIT"].insert(2000,plot.SAGRSUIT[2000]);  
  data.erase("AGRSUIT");                                               
  data["AGRSUIT"].insert(2000,plot.AGRSUIT[2000]);
  data.erase("PRICEINDEX");                                               
  data["PRICEINDEX"].insert(2000,plot.PRICEINDEX[2000]);       
  data.erase("DECHERB");                                                    
  data["DECHERB"].insert(2000,plot.DECHERB[2000]);  
  data.erase("DECWOOD");                                               
  data["DECWOOD"].insert(2000,plot.DECWOOD[2000]);            
  data.erase("DECSOC");                                               
  data["DECSOC"].insert(2000,plot.DECSOC[2000]); 
  data.erase("POLESREG");
  data["POLESREG"].insert(2000,plot.POLESREG[2000]); 
  data.erase("FTIMBER");     
  data["FTIMBER"].insert(2000,plot.FTIMBER[2000]);   
  data.erase("POPDENS"); 
  data.erase("SPOPDENS"); 
  data.erase("GDP");             
  data.erase("BUILTUP");             
  data.erase("CROP");
  for (int year=2000;year<=2100;year+=10) {  
    data["POPDENS"].insert(year,plot.POPDENS[year]);
    data["SPOPDENS"].insert(year,plot.SPOPDENS[year]);
    data["GDP"].insert(year,plot.GDP[year]);
    data["BUILTUP"].insert(year,plot.BUILTUP[year]);
    data["CROP"].insert(year,plot.CROP[year]);                           
  }
//cout<<plot.FTIMBER[2000]<<endl;
//cout << "ftimber=\t" <<data["FTIMBER"].v()<<"\t polesreg=\t"<<data["POLESREG"].v()<<"\t CABOVEHA=\t"<<data["CABOVEHA"].v()<<endl;
 }



#endif
