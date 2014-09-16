#include "clopts.hpp"

void make_empty_file(std::string dir,std::string file)
{
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::ofstream outfile (file.c_str());
  outfile.close();
  
}

void update_config_file(std::string config_file, po::variables_map &vm) 
{ 
  //std::cout << "Add missing values to the config file" << std::endl;  

  std::ifstream input(config_file.c_str());
  std::vector<std::string> lines;
  { // read each line and add to lines:
    std::string line;
    while(getline(input,line))
      lines.push_back(line);
  }
  input.close();

  std::ofstream outfile;
  //std::cout << config_file << std::endl;
  outfile.open(config_file.c_str(), std::ios_base::app);

  for(po::variables_map::iterator vit = vm.begin(); 
      vit != vm.end(); 
      ++vit) {
    std::string v_name = vit->first;
    //std::cout << v_name  << std::endl;
	
    if(v_name != "config") {
      // TODO: instead of just checking config, check against
      // disallowed config-file variables.

      bool var_in_file=false;
	
      for(unsigned int i=0; i < lines.size(); ++i)
	if(lines[i].find(vit->first) == 0)
	  var_in_file=true;
	  
      if(!var_in_file) {
	// The variable is not in the config file, so add it.

	//std::cout << "val missing in config; adding" << std::endl;

	// NB: we do a lot of casting and catching exceptions so
	// that we can cout boost::any.
	outfile  << vit->first
		 << "=" ;
	try { outfile << vit->second.as<double>();
	} catch(...) {/* do nothing */ }
	try { outfile << vit->second.as<int>();
	} catch(...) {/* do nothing */ }
	try { outfile << vit->second.as<std::string>();
	} catch(...) {/* do nothing */ }
	try { outfile  << vit->second.as<bool>();
	} catch(...) {/* do nothing */ }
	outfile << std::endl;
      }
    }
  }
  outfile.close();
}

