#ifndef INTERPOL_H
#define INTERPOL_H

#include <map>

/////////////////
// INTERPOLATE //
/////////////////
class interpol {
//Interpolate: first..the observation year, second..the observed value
public:
  void insert(double, double);
  double v(double);  //returns the value
  friend interpol operator *(interpol a, double b);
  void clear(double, double);
private:
  typedef std::map<double, double> MapType;
  typedef MapType::value_type ValuePair;
  MapType aMap;
  double vinterpol(double x1, double y1, double x2, double y2, double x);
};


void interpol::clear(double i, double d) {
  aMap.insert(ValuePair(i, d));
}

interpol operator *(interpol a, double b) {
  interpol tmp;
  typedef std::map<double, double> MapType;
  typedef MapType::value_type ValuePair;
  MapType::const_iterator iter = a.aMap.begin();
  while(iter != a.aMap.end()) {
    tmp.insert(iter->first, iter->second * b);
    ++iter;
  }
  return(tmp);
}

double interpol::vinterpol(double x1, double y1, double x2, double y2, double x) {
  //interpolate/extrapolate the y value for value x
  double y;
  if(x1 == x2) {y = (y1 + y2)/2.;}
  else {
    y = y1 + (x-x1)/(x2-x1) * (y2-y1);
  }
  return(y);
}

double interpol::v(double i=0.) {
  MapType::const_iterator lo, up;
  double y = 0.;

  if(aMap.size() > 0) {
    up = aMap.lower_bound(i);
    if(up == aMap.end()) {--up; lo = up;}
    else {
      lo = up;
      --lo;
      if(up == aMap.begin()) {lo = aMap.begin();}
    }
    y = vinterpol(lo->first, lo->second, up->first, up->second, i);
  }
  return(y);
}

void interpol::insert(double i, double d) {
  aMap.insert(ValuePair(i, d));
}

#endif
