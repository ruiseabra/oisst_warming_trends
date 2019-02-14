//  TargetCores.cpp
//  sstcoremulti

#include <stdio.h>
#include <iostream>
#include <string>
#include "Scanner.hpp"

using namespace std;

bool Scanner::TargetCores(const string& pattern){
  vector <string> CORES;
  int core;
  
  CORES = Str_split(pattern, ',');
  
  for (int i = 0; i < CORES.size(); i++) {
    core = stoi(CORES[i]) - 1;
    if (core < 0) return false;
    cores.push_back(core);
  }
  
  return true;
}
