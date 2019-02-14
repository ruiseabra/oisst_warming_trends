//  SSTfiles.cpp
//  sstcoremulti

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Scanner.hpp"

using namespace std;

bool Scanner::SSTfiles(const string& path) {
  ifstream infile(path);
  string str;
  
  if (!infile) {
    cerr << "failed" << endl << "Could not open the OISST paths' file " << path << endl;
    return false;
  }
  
  while (getline(infile, str)) {
    datafiles.push_back(str);
  }
  
  if(datafiles.size() == 0) return false;
  return true;
}
