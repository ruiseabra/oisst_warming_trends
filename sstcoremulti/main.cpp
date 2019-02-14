//  main.cpp
//  sstcoremulti
//
//  Created by Rui Seabra and António Múrias dos Santos on 03/09/2018.
//  Copyright © 2018 CIBIO. All rights reserved.
//
  
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include "Scanner.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  time_t    seconds {0};
  bool      status {true};
  Scanner   scanner;
  int       mins, secs;
  
  if(argc > 2){
    // get the list of sst files & target cores
    if(scanner.SSTfiles(argv[1]) && scanner.TargetCores(argv[2])) {
      seconds = time(NULL);
      status  = scanner.Run();
      seconds = time(NULL) - seconds;
      mins    = floor(seconds / 60);
      secs    = int(seconds) - (mins * 60);
      if(status) {
        cerr << "-----------------------------------------------------------------" << endl;
        cerr << "  sstcore runtime: " << mins << "m " << secs << "s" << endl;
        cerr << "---------------------------------------------------------------//" << endl;
      }
    }
  }else{
    cerr << endl << endl;
    cerr << "  sstcoremulti v1.0 [2018-08]" << endl;
    cerr << "------------------------------------------------------------------" << endl;
    cerr << "| reads files in binary format                                   |" << endl;
    cerr << "|    [integer, 2-bytes, ((x + 30) * 100) conversion]             |" << endl;
    cerr << "| -------------------------------------------------------------- |" << endl;
    cerr << "| run as: sstcoremulti <path> '<cores>'                          |" << endl;
    cerr << "| -------------------------------------------------------------- |" << endl;
    cerr << "| <path>  path for a txt file holding full paths for the         |" << endl;
    cerr << "|           target binary SST files to be processed              |" << endl;
    cerr << "| <core>  indexes for the cores to be extracted                  |" << endl;
    cerr << "|           - indexes start in 1, not 0                          |" << endl;
    cerr << "|           - provided as a string of cores connected by ','     |" << endl;
    cerr << "| -------------------------------------------------------------- |" << endl;
    cerr << "| -------------------------------------------------------------- |" << endl;
    cerr << "| Example:                                                       |" << endl;
    cerr << "| $ sstcoremulti ~/file_list.txt 1982 2017 '20560,20561,20563'   |" << endl;
    cerr << "------------------------------------------------------------------" << endl;
    cerr << endl;
  }
  return 0;
}
