void cPrices(void)
 {
/////////////////////// MAC2050 Price functions from Nat (DEFRA)
//
////MG: Definition of a table of linear price-time function for MACs 2050 scenario


	LinPrice2050[0]=0.;
	LinPrice2050[1]=0.;		LinPrice2050[2]=0.;		LinPrice2050[3]=0.;		LinPrice2050[4]=0.;	LinPrice2050[5]=1.21847E-14;LinPrice2050[6]=0.02;
	LinPrice2050[7]=0.04;	LinPrice2050[8]=0.07;	LinPrice2050[9]=0.09;	LinPrice2050[10]=0.11;	LinPrice2050[11]=0.13;	LinPrice2050[12]=0.16;
	LinPrice2050[13]=0.18;	LinPrice2050[14]=0.20;	LinPrice2050[15]=0.22;	LinPrice2050[16]=0.24;	LinPrice2050[17]=0.27;	LinPrice2050[18]=0.29;
	LinPrice2050[19]=0.31;	LinPrice2050[20]=0.33;	LinPrice2050[21]=0.36;	LinPrice2050[22]=0.38;	LinPrice2050[23]=0.40;	LinPrice2050[24]=0.42;
	LinPrice2050[25]=0.44;	LinPrice2050[26]=0.47;	LinPrice2050[27]=0.49;	LinPrice2050[28]=0.51;	LinPrice2050[29]=0.53;	LinPrice2050[30]=0.56;
	LinPrice2050[31]=0.58;	LinPrice2050[32]=0.60;	LinPrice2050[33]=0.62;	LinPrice2050[34]=0.64;	LinPrice2050[35]=0.67;	LinPrice2050[36]=0.69;
	LinPrice2050[37]=0.71;	LinPrice2050[38]=0.73;	LinPrice2050[39]=0.76;	LinPrice2050[40]=0.78;	LinPrice2050[41]=0.80;	LinPrice2050[42]=0.82;
	LinPrice2050[43]=0.84;	LinPrice2050[44]=0.87;	LinPrice2050[45]=0.89;	LinPrice2050[46]=0.91;	LinPrice2050[47]=0.93;	LinPrice2050[48]=0.96;
	LinPrice2050[49]=0.98;	LinPrice2050[50]=1.;

		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){CubPrice2050[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2050){CubPrice2050[i]=pow((i+2000.-2005.),3)/pow((2050.-2005.),3);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2050){CubPrice2050[i]=1.;};
			}

		for (int i=50;i>=0;i--){
			float currYear = float(i) + 2000.;
			float oneThird = 1./3.;
		if (i+2000 >= 2050){MixPrice2050[i]=1.;};
		if (i+2000 >= 2010 && i+2000 < 2050){MixPrice2050[i]=pow(float(i+2000-2005),oneThird)/3.55689330449006;};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 == 2010){MixPrice2050[i]=(pow((currYear-2005),oneThird)/2.46621207433047)*MixPrice2050[20];};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2006 && i+2000 < 2010){MixPrice2050[i]=((currYear-2005.)/5.)*MixPrice2050[10];};
		if (i+2000 < 2006){MixPrice2050[i]=0.;};
		}
///////////////////// MAC2030 Price functions
		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){LinPrice2030[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2030){LinPrice2030[i]=(i+2000.-2005.)/(2030.-2005.);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2030){LinPrice2030[i]=1.;};
			}
		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){CubPrice2030[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2030){CubPrice2030[i]=pow((i+2000.-2005.),3)/pow((2030.-2005.),3);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2030){CubPrice2030[i]=1.;};
			}

		for (int i=50;i>=0;i--){
		float currYear = float(i) + 2000.;
		float oneThird = 1./3.;
		if (i+2000 >= 2030){MixPrice2030[i]=1.;};
		if (i+2000 >= 2010 && i+2000 < 2030){MixPrice2030[i]=pow(float(i+2000-2005),oneThird)/2.9240177382;};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 == 2010){MixPrice2030[i]=(pow((currYear-2005),oneThird)/2.46621207433047)*MixPrice2030[20];};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2006 && i+2000 < 2010){MixPrice2030[i]=((currYear-2005.)/5.)*MixPrice2030[10];};
		if (i+2000 < 2006){MixPrice2030[i]=0.;};
		}
///////////////////// MAC2020 Price functions
		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){LinPrice2020[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2020){LinPrice2020[i]=(i+2000.-2005.)/(2020.-2005.);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2020){LinPrice2020[i]=1.;};
			}

		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){CubPrice2020[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2020){CubPrice2020[i]=pow((i+2000.-2005.),3)/pow((2020.-2005.),3);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2020){CubPrice2020[i]=1.;};
			}

		for (int i=50;i>=0;i--){
			float currYear = float(i) + 2000.;
			float oneThird = 1./3.;
		if (i+2000 >= 2020){MixPrice2020[i]=1.;};
		if (i+2000 >= 2010 && i+2000 < 2020){MixPrice2020[i]=pow(float(i+2000-2005),oneThird)/2.46621207433;};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 == 2010){MixPrice2020[i]=(pow((currYear-2005),oneThird)/2.46621207433047)*MixPrice2020[20];};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2006 && i+2000 < 2010){MixPrice2020[i]=((currYear-2005.)/5.)*MixPrice2020[10];};
		if (i+2000 < 2006){MixPrice2020[i]=0.;};
		}
				///////////////////// MAC2015 Price functions
		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){LinPrice2015[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2015){LinPrice2015[i]=(i+2000.-2005.)/(2015.-2005.);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2015){LinPrice2015[i]=1.;};
			}

		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){CubPrice2015[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2015){CubPrice2015[i]=pow((i+2000.-2005.),3)/pow((2015.-2005.),3);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2015){CubPrice2015[i]=1.;};
			}

		for (int i=50;i>=0;i--){
			float currYear = float(i) + 2000.;
			float oneThird = 1./3.;
		if (i+2000 >= 2015){MixPrice2015[i]=1.;};
		if (i+2000 >= 2010 && i+2000 < 2015){MixPrice2015[i]=pow(float(i+2000-2005),oneThird)/2.15443469;};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
//		if (i+2000 == 2010){MixPrice2015[i]=(pow((currYear-2005),oneThird)/2.15443469)*MixPrice2015[15];};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2006 && i+2000 < 2010){MixPrice2015[i]=((currYear-2005.)/5.)*MixPrice2015[10];};
		if (i+2000 < 2006){MixPrice2015[i]=0.;};
		}
				///////////////////// MAC2010 Price functions
		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){LinPrice2010[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2010){LinPrice2010[i]=(i+2000.-2005.)/(2010.-2005.);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2010){LinPrice2010[i]=1.;};
			}

		for (int i=0;i<=50;i++){
		if (i+2000 < 2006){CubPrice2010[i]=0;};
		if (i+2000 >= 2006 && i+2000 < 2010){CubPrice2010[i]=pow((i+2000.-2005.),3)/pow((2010.-2005.),3);};//{CubPrice=(year-2005)^3/(2050-2005)^3;};
		if (i+2000 >= 2010){CubPrice2010[i]=1.;};
			}

		for (int i=50;i>=0;i--){
			float currYear = float(i) + 2000.;
			float oneThird = 1./3.;
		if (i+2000 >= 2010){MixPrice2010[i]=1.;};
		if (i+2000 >= 2006 && i+2000 < 2010){MixPrice2010[i]=((currYear-2005.)/5.);};
		if (i+2000 < 2006){MixPrice2010[i]=0.;};
		}
 }
