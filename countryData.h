//  Name: CountryData class with statistics by regions
//  Author: Andriy Bun
//  Date: 12.10.09
//  MG: added printToFile
//  Description: 


#ifndef COUNTRYDATA_H_
#define COUNTRYDATA_H_

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "misc.h"

const int countriesNum = 242;

using namespace std;

class countryData {
private:
    vector< g4m::ipol<double,double> > countries;
    vector< unsigned char > regions;
    void setRegions();
public:
    countryData();
    ~countryData();
    void reset(void);
    void set(int, int, double);
    void inc(int, int, double);
    double get(int, int);
    double getRegionSum(unsigned char regID, int year);
    void PrintToFile(string, int, int,int);      // print array to file    
};

void countryData::reset(void)
{
  int ind = 0;
  do {
    countries[ind].clear();
    ind++;
  } while (ind <countriesNum);
}


void countryData::set(int ind, int year, double val)
{
  countries[ind].insert(year,val);
}

void countryData::inc(int ind, int year, double val)
{
  countries[ind].inc(year,val);
}

double countryData::get(int ind, int year)
 {
  return(countries[ind][year]);
 }

double countryData::getRegionSum(unsigned char regID, int year)
 {
  double sum = 0;
  int ind = 0;
  do {
    if (regions[ind] == regID) {
      sum += get(ind, year);
    }
    ind++;
  } while (ind <countriesNum);  
  return sum;
 }

void countryData::PrintToFile(string fileName, int firstYear, int lastYear, int step)
 {
  ofstream f;
  fileName = fileName + ".txt";
  f.open(fileName.c_str(),ios::out);
  if (f.is_open()) {

    f << "CountryN/Year"; 
    for (int j = firstYear; j <= lastYear; j += step)
    {    
    f <<"\t"<<j ;}
    f << endl;

//-----------------------------        
    for (int i = 0; i < countriesNum; i++)
     {
    f << i ;
    for (int j = firstYear; j <= lastYear; j += step)
    {
        f<<"\t"<<countries[i][j];
        }
      f << endl;  
     }
    f.close(); 
  } else {
    cout << "Unable to save to file!" << endl;
  }
 }

 
// default constructor
countryData::countryData()
{
  countries.resize(countriesNum);
  regions.resize(countriesNum);
  setRegions();
}

// destructor
countryData::~countryData()
{

}

void countryData::setRegions()
 {
  regions.assign(countriesNum, 0);
    regions[1] = 20;
    regions[2] = 17;
    regions[3] = 14;
    regions[5] = 25;
    regions[8] = 19;
    regions[9] = 10;
    regions[10] = 1;
    regions[11] = 7;
    regions[12] = 10;
    regions[13] = 16;
    regions[15] = 20;
    regions[16] = 16;
    regions[17] = 7;
    regions[18] = 16;
    regions[19] = 25;
    regions[21] = 20;
    regions[22] = 19;
    regions[23] = 17;
    regions[24] = 25;
    regions[25] = 2;
    regions[26] = 21;
    regions[27] = 6;
    regions[28] = 25;
    regions[29] = 25;
    regions[30] = 10;
    regions[31] = 22;
    regions[32] = 25;
    regions[33] = 3;
    regions[35] = 25;
    regions[36] = 25;
    regions[37] = 19;
    regions[38] = 4;
    regions[39] = 19;
    regions[40] = 25;
    regions[41] = 25;
    regions[42] = 16;
    regions[43] = 17;
    regions[44] = 16;
    regions[45] = 9;
    regions[46] = 6;
    regions[47] = 8;
    regions[49] = 16;
    regions[50] = 16;
    regions[51] = 19;
    regions[53] = 16;
    regions[53] = 25;
    regions[55] = 25;
    regions[56] = 5;
    regions[57] = 25;
    regions[60] = 15;
    regions[61] = 8;
    regions[62] = 7;
    regions[65] = 25;
    regions[65] = 25;
    regions[68] = 10;
    regions[69] = 7;
    regions[70] = 25;
    regions[71] = 9;
    regions[73] = 16;
    regions[75] = 16;
    regions[76] = 25;
    regions[76] = 25;
    regions[78] = 19;
    regions[79] = 16;
    regions[80] = 16;
    regions[82] = 6;
    regions[84] = 11;
    regions[85] = 21;
    regions[86] = 14;
    regions[87] = 14;
    regions[89] = 8;
    regions[91] = 14;
    regions[92] = 9;
    regions[93] = 25;
    regions[94] = 16;
    regions[96] = 12;
    regions[98] = 10;
    regions[99] = 25;
    regions[103] = 22;
    regions[103] = 24;
    regions[104] = 14;
    regions[105] = 10;
    regions[106] = 22;
    regions[107] = 5;
    regions[108] = 14;
    regions[109] = 25;
    regions[110] = 25;
    regions[111] = 14;
    regions[113] = 5;
    regions[114] = 7;
    regions[116] = 17;
    regions[117] = 25;
    regions[118] = 25;
    regions[119] = 21;
    regions[121] = 25;
    regions[123] = 25;
    regions[123] = 25;
    regions[126] = 13;
    regions[127] = 10;
    regions[129] = 22;
    regions[130] = 14;
    regions[131] = 21;
    regions[131] = 25;
    regions[133] = 25;
    regions[134] = 20;
    regions[135] = 7;
    regions[137] = 1;
    regions[138] = 16;
    regions[140] = 25;
    regions[142] = 18;
    regions[143] = 14;
    regions[144] = 20;
    regions[145] = 16;
    regions[146] = 15;
    regions[147] = 19;
    regions[148] = 19;
    regions[149] = 21;
    regions[150] = 6;
    regions[151] = 9;
    regions[155] = 6;
    regions[156] = 10;
    regions[157] = 15;
    regions[157] = 25;
    regions[159] = 25;
    regions[160] = 14;
    regions[161] = 25;
    regions[162] = 25;
    regions[163] = 25;
    regions[164] = 21;
    regions[165] = 6;
    regions[166] = 6;
    regions[167] = 15;
    regions[168] = 25;
    regions[169] = 23;
    regions[170] = 9;
    regions[171] = 20;
    regions[173] = 16;
    regions[174] = 16;
    regions[175] = 25;
    regions[176] = 19;
    regions[178] = 25;
    regions[179] = 8;
    regions[180] = 18;
    regions[181] = 14;
    regions[183] = 10;
    regions[183] = 25;
    regions[185] = 21;
    regions[186] = 15;
    regions[186] = 25;
    regions[188] = 16;
    regions[189] = 14;
    regions[190] = 26;
    regions[193] = 25;
    regions[194] = 10;
    regions[195] = 14;
    regions[196] = 8;
    regions[197] = 27;
    regions[198] = 19;
    regions[199] = 10;
    regions[200] = 15;
    regions[201] = 19;
    regions[202] = 22;
    regions[205] = 14;
    regions[206] = 17;
    regions[207] = 25;
    regions[208] = 25;
    regions[209] = 25;
 }
#endif
