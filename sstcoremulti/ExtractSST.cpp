//  ExtractSST.cpp
//  sstcoremulti

#include <stdio.h>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include "Scanner.hpp"

using namespace std;

bool Scanner::ExtractSST() {
  ifstream infile;
  
  unsigned short val;
  unsigned int i, ii;
  int core;
  
  for (int f = 0; f < nfiles; f++) {
    i  = 0;
    ii = 0;
    infile.open(datafiles[f], ios::in|ios::binary);
    
    if (!infile) {
      cerr << "failed" << endl << "Could not open " << datafiles[f] << endl;
      return false;
    }
    
    core = cores[ii];
    while (infile.read((char *)&val, sizeof(unsigned short))) {
      if (i == core) {
        if (core != cores.front()) {
          cout << ',';
        }
        cout << val;
        if (core == cores.back()) break;
        ii++;
        core = cores[ii];
      }
      i++;
    }
    infile.close();
    
    cout << endl;
  }
  return true;
}
