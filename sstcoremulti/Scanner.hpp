//  Scanner.hpp
//  sstcoremulti

#ifndef Scanner_hpp
#define Scanner_hpp

#include <stdio.h>
#include <string>
#include <vector>


class Scanner{
  // paths for all binary sst files to process (obtained from argv[1])
  std::vector<std::string> datafiles;
  
  // number of sst files to be processed
  int nfiles {0};
  
  // indexes of the cores to be extracted (obtained from argv[2])
  // in R indexes start from 1, but here they start from 0, so they must be adjusted
  //   e.g. if in R we want the core[1000], here we will get core[999]
  std::vector<int> cores;
  
  // function to break strings using a pattern
  //   similar to str_split in R
  std::vector<std::string> Str_split(const std::string&, char);
  
  // function to extract sst data for the target cores from each OISST file in 'datafiles'
  bool ExtractSST();
  
  public:
  Scanner();
  ~Scanner();
  
  bool SSTfiles(const std::string&);
  bool TargetCores(const std::string&);
  bool Run();
};

#endif /* Scanner_hpp */
