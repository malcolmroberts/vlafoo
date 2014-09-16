#ifndef __clopts_hpp__
#define __clopts_hpp__ 1


#include <iostream>
#include <fstream>

#include <sys/types.h> // for making directories
#include <sys/stat.h>


#include <boost/program_options.hpp>
namespace po = boost::program_options;

void make_empty_file(std::string dir,std::string file);
void update_config_file(std::string config_file, po::variables_map &vm);

#endif
