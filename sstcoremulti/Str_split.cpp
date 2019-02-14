//  Str_split.cpp
//  sstcoremulti

#include <stdio.h>
#include <iostream>
#include "Scanner.hpp"

using namespace std;

vector<string> Scanner::Str_split(const string& s, char c) {
  auto end   = s.cend();
  auto start = end;
  vector<string> v;
  for(auto it = s.cbegin(); it != end; ++it) {
    if(*it != c) {
      if(start == end)
        start = it;
      continue;
    }
    if(start != end) {
      v.emplace_back(start, it);
      start = end;
    }
  }
  if(start != end)
    v.emplace_back(start, end);
  return v;
}
