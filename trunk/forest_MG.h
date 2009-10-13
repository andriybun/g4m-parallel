using namespace std;

struct dat {
  int Rotation;
  double LandAreaHa;
  double potHarvest;
  double forestShare;
  double OforestShare;
  double AforestShare;
  double prevOForShare;
  double prevOForShareRP;
  double prevAForShareRP;
  double AforestSharePrev;
  double savedCarbonPrev;
  double gainedCarbonPrev;
  double EmissionsTotPrev;
  double EmissionsAfforPrev;
  double prevPlantPhytHaBmGr;
  double prevPlantPhytHaBlGr;
  double deforestHaTot;
  double afforestHaTot;
  double EmissionsProduct;  
  double EmissionsLitter;  
  double EmissionsSOC;      
  double EmissionsSlashBurn;
  double EmissionsDeadBurn;
  double EmissionsCRootBurn;    
  double EmissionsTot;     
  double EmLitterAffor;
  double EmSOCAffor; 
  double EmissionsAffor;
  double forestAgeShare[110];
  double BDeadA[110];
  double LitterA[110];
  double SOCA[110];
  double ProdLongA[110];
  double ProdShortA[110];
  double deforestA[110];  
  double FineRootA[110];
  double LitterAffor[110];
  double SOCaffor[110];
  double prevReportYear;
  int ireportYear;
  };
//******************************************************************************
// types
//******************************************************************************
typedef vector< g4m::ipol<double,double> > cellCol;
typedef vector<cellCol> dataArray;
typedef vector<g4m::dataStruct> dataDetStruct;
typedef vector<g4m::ageStruct> ageStructVector;
typedef vector<dat> datGlobal;
//******************************************************************************
// containers of data
//******************************************************************************
//dataArray data05x_BUILTUP;
//dataArray data05x_CROP;
//dataArray data05x_AREA;
//dataArray data05x_PRICEINDEX;
//dataArray data05x_POPDENS;
//dataArray data05x_DISCRATE;
//dataArray data05x_GDP;
//dataArray data05x_COUNTRY;
//dataArray data05x_NOTUSABLE;
vector<g4m::dataStruct> dataDet;
g4m::coeffStruct coeff;
map<string, interpol> lprice; //datamap for land price corrections for current price  scenario (GLOBIOM)
map<string, interpol> wprice; //datamap for wood price corrections for current price  scenario (GLOBIOM)
map<string, interpol> wprod; //datamap for wood production corrections for POLES regions	(GLOBIOM)
//******************************************************************************
// constants and variables
//******************************************************************************
#ifdef unix
string homeDir = "./data/";
#else
string homeDir = "data\\";
#endif
int ResLatitude, ResLongitude;    // resolutions of model
int eyear, byear;
const double GridStepLat = 0.5;   // step by latitude
const double GridStepLon = 0.5;   // step by longitude
const int nYears = 110;
const int NumberOfCountries = 209;
double MAI_CountryUprotect[209];
double MAI_CountryAll[209];
//double woodHarvest[209]; 
double woodHarvestStat[209];
double Hurdle_opt[209];
double afforRate_opt[209];
double deforRate_opt[209];
double EmissionsCurCountry[210];
double EmissionsCurAfforCountry[210];
double LinPrice2050[51];
float CubPrice2050[51];
float MixPrice2050[51];
float LinPrice2030[51];
float CubPrice2030[51] ;
float MixPrice2030[51];
float LinPrice2020[51] ;
float CubPrice2020[51] ;
float MixPrice2020[51];
float LinPrice2015[51] ;
float CubPrice2015[51] ;
float MixPrice2015[51];
float LinPrice2010[51] ;
float CubPrice2010[51] ;
float MixPrice2010[51];
int PriceCiS[12] = {0,10,20,30,50,70,100,200,300,500,1000,0};
double deflator = 0.8807;

//*****************************************************************************
// country outputs
//*****************************************************************************
  countryData CountriesNforCover = countryData();
  countryData CountriesNforTotC = countryData();
  countryData CountriesAfforHaYear = countryData();  
  countryData CountriesAfforCYear = countryData();  
  
  countryData CountriesOforCover = countryData();
  countryData CountriesOforBiomassC = countryData();
  countryData CountriesOforTotC = countryData();  
  countryData CountriesDeforHaYear = countryData();  
  countryData CountriesDeforCYear = countryData();   
  
  countryData CountriesWoodHarvestM3Year = countryData();     
  countryData CountriesWoodLoosCYear = countryData();   

  countryData CountriesManagedForHa = countryData();     
  countryData CountriesManagedCount = countryData();   

//******************************************************************************
// functions
//******************************************************************************
void readInput(string, dataArray &);
int readInputDet(dataDetStruct &);
//void formInputDet(dataDetStruct &,int [],double []);

vector<double> procPlots(g4m::dataStruct &,double,int,double, double,double, double, double);
vector<double> calcPlots(double [],int);
vector<double> plotsData(g4m::dataStruct &,int);
void initManagedForest(dataDetStruct &, g4m::incrementTab &, datGlobal &,
                       griddata &, griddata2<char> &, griddata2<char> &, 
                       griddata2<char> &, griddata &,griddata &);
 void adjustManagedForest(dataDetStruct &data_all, g4m::incrementTab &fi, ageStructVector &cohort_all, 
              ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
              griddata2<char> &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
              griddata &rotationForestNew, griddata2<char> &thinningForestNew,
              griddata2<char> &rotationType, griddata &harvestGrid, int year); 
void initLoop(int, dataDetStruct &, g4m::incrementTab &, ageStructVector &, 
              ageStructVector &, datGlobal &, griddata &, griddata2<char> &, griddata &);
//void calc(g4m::dataStruct &, g4m::incrementTab &, g4m::ageStruct &, g4m::ageStruct &, 
//          dat &, griddata2<char> &, griddata &, griddata &, int, int, int);
void calc(g4m::dataStruct &, g4m::incrementTab &, g4m::ageStruct &, g4m::ageStruct &,
          dat &, griddata2<char> &, griddata &, griddata &,
          griddata &, griddata2<char> &, griddata2<char> &, 
          griddata &, int , int , int );          
void MAI_country(void);
void woodHarvestStatCountry(void);
void hurdle_aff_deff(void);
void cPrices(void);
void int2str(int i, char tmp[]);   // defined in misc.h

