//  Run.cpp
//  sstcoremulti

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include "Scanner.hpp"

using namespace std;

bool Scanner::Run(){
  
  nfiles = int(datafiles.size());
  
  cerr << endl << endl;
  cerr << "  sstcoremulti v1.1 [2018-09]" << endl;
  cerr << "-----------------------------------------------------------------" << endl;
  cerr << "  processing " << nfiles << " binary OISST data files" << endl;
  cerr << "  from " << endl;
  cerr << "     " << datafiles.front() << endl;
  cerr << "  to " << endl;
  cerr << "     " << datafiles.back() << endl;
  
  if(!ExtractSST()) {
    return false;
  }
  
  return true;
}
