#include "clopts.hpp"

// Show the variables_map used in boost program_options
void show_vm(po::variables_map vm)
{
  for(po::variables_map::iterator vit = vm.begin(); 
      vit != vm.end(); 
      ++vit) {
    std::cout << vit->first << "=";
    try { std::cout << vit->second.as<double>();
    } catch(...) {/* do nothing */ }
    try { std::cout << vit->second.as<int>();
    } catch(...) {/* do nothing */ }
    try { std::cout << vit->second.as<std::string>();
    } catch(...) {/* do nothing */ }
    try { std::cout  << vit->second.as<bool>();
    } catch(...) {/* do nothing */ }
    std::cout << std::endl;
  }
}

void make_empty_file(std::string dir, std::string file)
{
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string dir_file=dir+"/"+file;
  std::ofstream outfile(dir_file.c_str());
  outfile.close();
}

void fill_config_file(std::string config_dir_file, po::variables_map &vm) 
{ 
  //std::cout << "Add missing values to the config file" << std::endl;  

  std::ifstream input(config_dir_file.c_str());
  std::vector<std::string> lines;
  { // read each line and add to lines:
    std::string line;
    while(getline(input,line))
      lines.push_back(line);
  }
  input.close();

  std::ofstream outfile;
  //std::cout << config_dir_file << std::endl;
  outfile.open(config_dir_file.c_str(), std::ios_base::app);

  for(po::variables_map::iterator vit = vm.begin(); 
      vit != vm.end(); 
      ++vit) {
    std::string v_name = vit->first;
    //std::cout << v_name  << std::endl;
	
    // the "" is to get rid of positional arguments.
    if(v_name != "config" &&  v_name != "outdir" && v_name != "") {
      // TODO: instead of just checking config, check against
      // disallowed config-file variables.

      bool var_in_file=false;
	
      for(unsigned int i=0; i < lines.size(); ++i)
	if(lines[i].find(vit->first) == 0)
	  var_in_file=true;
	  
      if(!var_in_file) {
	// The variable is not in the config file, so add it.

	// std::cout << "val missing in config; adding" << std::endl;
	// std::cout << v_name << std::endl;
	// std::cout << vit->first << std::endl << std::endl;

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

void update_config_file(std::string config_dir_file, po::parsed_options &cl_opts,
			po::variables_map &vm)
{
  // Update the config file with the values from the command line.
  {
    //std::cout << "update config file with values from command-line" 
    //	<< std::endl;

    // Read the config file and store it, line-by-line, in lines.
    std::ifstream input(config_dir_file.c_str());
    std::vector<std::string> lines;
    { // read each line and add to lines:
      std::string line;
      while(getline(input,line))
	lines.push_back(line);
    }
    input.close();

    /*
      std::cout << "Lines from config file:" << std::endl;
      for(unsigned int i=0; i < lines.size(); ++i)
      std::cout << lines[i] << std::endl;
      std::cout << std::endl;
    */

    // Re-write the config file, line-by-line, replacing old
    // values with new values from the command-line.
    std::ofstream outfile;
    outfile.open(config_dir_file.c_str());
      
    for(unsigned int i=0; i < lines.size(); ++i) {
      // Search for updates from command-line:
      bool isfound=false;
      for(vec_opt::iterator ipo = cl_opts.options.begin();
	  ipo != cl_opts.options.end(); 
	  ++ipo) {
	po::basic_option<char>& l_option = *ipo;

	//std::cout << l_option.string_key << std::endl;

	// TODO improve the check if var name is allowed
	if(l_option.string_key != "config" 
	   && l_option.string_key != "outdir"
	   && l_option.string_key != "") {
	  if(!isfound) {

	    int found = lines[i].find(l_option.string_key);
	    isfound=(found == 0);

	    if(isfound) {
	      //std::cout <<  l_option.string_key << std::endl;
	      //std::cout << "\tadding to config file" << std::endl;
	      // Replace the line with the new value
	      outfile << l_option.string_key 
		      << "="
		      << l_option.value[0] << std::endl;
	      /*
		std::cout << "command line:\t"
		<< l_option.string_key 
		<< "=" 
		<< l_option.value[0] << std::endl;
		std::cout << "\tfound " << l_option.string_key << std::endl;
		std::cout << config_dir_file << ":\t" << lines[i] << std::endl;
	      */
	      //ipo=cl_opts.options.end(); // stop looping.
	    }
	  }
	}
      }
	
      if(!isfound) {
	// Keep the original line
	//std::cout << "keep old line:" << std::endl;
	//std::cout << lines[i] << std::endl;
	outfile << lines[i] << std::endl;
      }
	
    }
    outfile.close();
  }
}

void make_asy_header(po::variables_map vm)
{
  std::string asyfile="vlafoo.asy";

  std::ofstream outfile;
  outfile.open(asyfile.c_str());

  for(po::variables_map::iterator vit = vm.begin(); 
      vit != vm.end(); 
      ++vit) {
    try {  // doubles
      vit->second.as<double>();
      outfile << "real";
    } catch(...) {/* do nothing */ }
    try { 
      vit->second.as<int>();
      outfile << "int";
    } catch(...) {/* do nothing */ }
    try { 
      vit->second.as<std::string>();
      outfile << "string";
    } 
    catch(...) {/* do nothing */ }
    try { 
      vit->second.as<bool>();
      outfile << "bool";
    } 
    catch(...) {/* do nothing */ }
    outfile << " " <<  vit->first << ";" << std::endl;
  }

  outfile.close();
}
void make_asy_input(std::string dir, po::variables_map vm)
{
  std::string asyfile="vlafoo.asy";
  make_empty_file(dir, asyfile);

  std::string asy_dir_file=dir+"/"+asyfile;
  std::ofstream outfile;
  outfile.open(asy_dir_file.c_str());

  for(po::variables_map::iterator vit = vm.begin(); 
      vit != vm.end(); 
      ++vit) {
    outfile <<  vit->first << "=";
    try {  // doubles
      outfile << vit->second.as<double>();
    } catch(...) {/* do nothing */ }
    try { 
      outfile << vit->second.as<int>();
    } catch(...) {/* do nothing */ }
    try { 
      outfile << "\"" 
	      << vit->second.as<std::string>()
	      << "\"" ;
    } 
    catch(...) {/* do nothing */ }
    try { 
      if(vit->second.as<bool>())
	outfile << "true";
      else
	outfile << "false";
    } 
    catch(...) {/* do nothing */ }
    outfile << ";" << std::endl;
  }

  outfile.close();
}
