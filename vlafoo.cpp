#include "vlafoo.hpp"
#include <getopt.h>
#include <stdlib.h>     /* atoi */
#include <iomanip>      /* setw */
#include "xstream.h"
#include <fstream>

#include <iterator>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef std::vector< po::basic_option<char> > vec_opt; 

// Compute initial condition
void VlaFoo::initial_conditions(std::string & ic){
  real overs2pi=1.0/sqrt(2.0*PI);

  f.Load(0.0);
  
  // for(int ix=0;ix<nx;ix++) f[ix][1]=1.0;

  if(ic == "test0") {
    f[0][nv/4]=1.0;
    f[0][3*nv/4]=1.0;
    return;
  }
  if(ic == "test1") {
    for(int j=0; j < nv; ++j)
      f[0][j]=1.0;
    return;
  }

  if(ic == "landau") {
    std::cout << "Landau initial conditions." << std::endl;

    // landau damping
    // f0=exp(-vt*vt/2)*(1+eps*cos(kx*xt))/sqrt(2*pi)
    for(int i=0; i < nx; i++){
      real x=i*dx;
      real xterm=(1.0+eps*cos(kx*x))*overs2pi;
      for(int j=0; j < nv; j++){
	real v=index2v(j);
	f[i][j]=exp(-v*v*0.5)*xterm;
      }
    }
    return;
  }

  if(ic == "doublestream") {
    std::cout << "Double-stream instability initial conditions." 
	      << std::endl;
    //double-stream instability
    // f0=(exp(-(v-v0)*(v-v0)/2)+exp(-(v+v0)*(v+v0)/2))*
    // (1+eps*cos(kx*x))/sqrt(2*pi)/2
    real v0=3.0;
    real overs2pi=0.5/sqrt(2.0*PI);
    for(int i=0; i < nx; i++){
      real x=i*dx;
      real xterm=(1.0+eps*cos(kx*x))*overs2pi;
      for(int j=0; j < nv; j++){
	real v=index2v(j);
	real sum=v+v0;
	real diff=v-v0;
	f[i][j]=xterm*(exp(-0.5*sum*sum) + exp(-0.5*diff*diff));
      }
    }
    return;
  }
  
  std::cout << "Abort: no initial conditions found matching " 
	    << ic << std::endl;
  exit(1);
  
}

real VlaFoo::interpolate(real x, int nq, real* x0, real* y0)
{
  // Lagrange polynomial interpolation
  real val=0;
  if(nq > 0) {
    // NB: this is probably not the most numerically robust way to do
    // this.
    for(int j=0; j < nq; ++j) {
      real lj=1.0;
      real xj=x0[j];
      for(int m=0; m < j; ++m)
	lj *= (x-x0[m])/(xj-x0[m]);
      for(int m=j+1; m < nq; ++m)
	lj *= (x-x0[m])/(xj-x0[m]);
      val += y0[j]*lj;
      
    }
  } else {
    std::cout << "Error: spline of order "<< nq << " not implemented."
	      << std::endl;
    exit(1);
  }

  return val;
}

void VlaFoo::transport_x(real &dtx)
{
  // Advance f in time, taking accont the A\partial_x f term, with
  // A=diag[v].

  // The velocity for the particles is given by the distribution, ie
  // particles in the array at (i,j) have position x_i and velocity
  // v(j).

  real *x0=new real[nq];
  real *f0=new real[nq];

  real overdx=1.0/dx;

  // f(x,t+dt) = f(x-v*dt,t)
  for(int i=0; i < nx; i++) {
    for(int j=0; j < nv; j++) {
      real x=i*dx; // current position
      real v=index2v(j);
      if(x < 0.0) x += L;
      if(x > L) x -= L;

      real xold=x-v*dtx; // previous position
      if(xold < 0.0) xold += L;
      if(xold > L) xold -= L;
      
      // Find the neighbouring points:
      int iold=(int)floor(xold*overdx);
      int nleft=nq/2-1;
      int ileft=iold-nleft;
      if(ileft < 0) ileft += nx;
      x0[0]=dx*(ileft);
      for(int q=1; q < nq; q++)
	x0[q]=x0[0]+q*dx;

      for(int q=0; q < nq; q++) {
	if(abs(x0[q]-xold) > 0.5*L)
	  x0[q] -= L;
      }   

      // Find the f-values for the interpolation:
      for(int q=0; q < nq; q++) {
	//int qi=mod(nleft+i,nx);
	int qi=(ileft+q)%nx;

	f0[q]=f[qi][j]; // NB: periodic in x.
      }

      // Put the updated f in S:
      real val=interpolate(xold,nq,x0,f0);
      S[i][j]=val;
    }
  }
  delete[] x0;
  delete[] f0;
  

  // TODO: optimise into a swap?
  for(int i=0; i < nx; ++i) {
    array1<real>::opt fi=f[i];
    array1<real>::opt Si=S[i];
    for(int j=0; j < nv; ++j)
      fi[j]=Si[j];
  }
}

void VlaFoo::compute_rho(array2<real> &f)
{  
  // \rho = \int \approx  dv \sum_j f_j(x,t)
  for(int i=0; i < nx; i++) {
    real rhoi=0.0;
    array1<real>::opt fi=f[i];
    for(int j=0; j < nv; j++)
      rhoi += fi[j];
    rho[i]=rhoi*dv;
  }
}

void VlaFoo::compute_E(array2<real> &f)
{
  // Compute
  // E(x) = \int dx \int_{-vmax}^{vmax} dv f(x,v)
  // from
  // \del\cdot E = -1 + \rho
  // via a FFT.
  
  compute_rho(f);

  for(int i=0; i<nx; i++) rho[i] -= 1.0;

#if 1
  // using Fourier integration

  rc_x->fft(rho,rhok);

  // Set mean to zero:
  rhok[0]=0.0;
  real overnx=1.0/nx;
  Complex Iovernx=Complex(0,overnx);
  int stop=nxk-1;
  for(int i=1; i < stop; i++) {
    real k0xi=k0x*i;
    rhok[i] *= -Iovernx/k0xi;
  }
  rhok[nxk-1]=0.0; // kill the nyquist

  cr_x->fft(rhok,E);
#else
  // using crazy cotan stuff
  
  rc_x->fft(rho,rhok);

  // Set mean to zero:
  rhok[0]=0.0;
  real overnx=1.0/nx;
  for(int i=1; i < nxk-1; i++) {
    real k=k0x*i;
    Complex rhoki=rhok[i];
    real re=rhoki.real();
    real im=rhoki.imag();
    real ct=1.0/tan(PI*k*overnx);
    rhok[i]=0.5*k*dx*ct*Complex(im,-re);
  }
  rhok[nxk-1]=0.0; // kill the nyquist

#endif
}

void VlaFoo::compute_dfdv(array2<real> &f, array2<real> &S)
{
  // compute \partial_v f_i via Fourier differentiation.
  
  double overnv=1.0/nv; // FFT normalization term.
 
  Complex I=Complex(0,1);

  for(int i=0; i<nx; i++) {
    // copy velocity for each i into a temp buffer (NB real to complex
    // FFT overwrites input):
    array1<real>::opt fi=f[i];
    array1<real>::opt Si=S[i];

    for(int j=0; j < nv; ++j)
      vtemp[j]=fi[j];

    rc_v->fft(vtemp,vtempk);
    unsigned int stop=nvk-1;
    for(unsigned j=0; j < stop; j++) {
      real k=j*k0v;
      vtempk[j] *= I*k*overnv;
    }
    vtempk[stop]=0.0;  // Set Nyquist to zero.
    cr_v->fft(vtempk,Si); // transform back and put into source term
  }
}


// Compute -E_i \partial_v f_i
void VlaFoo::compute_negEdfdv(array2<real> &f, array2<real> &S)
{
  compute_E(f);
  compute_dfdv(f,S);

  for(int i = 0; i < nx; ++i) {
    array1<real>::opt Si = S[i];
    real Ei = E[i];
    for(int j = 0; j < nv; ++j)
      Si[j] *= -Ei;
  }
}

void VlaFoo::rk_source(real *f0, real *S0)
{
  // wrapper for compute_negEdfdv.

  array2<real> f(nx,nv,f0);
  array2<real> S(nx,nv,S0);
  compute_negEdfdv(f,S);
}

// Time step dtv of velocity resolution by Fourier
void VlaFoo::transport_v(array2<real> &f, real &dtv)
{
  // source term for f_i is -E_i \partial_v f_i and is called via
  // VlaFoo::rk_source.

  rk_step(f(),dtv);
}

void VlaFoo::time_step(real &dt)
{
  // Time-step by Strang splitting.

  // Advection term from rho (mean velocity) in space:
  transport_x(dt);

  // Compute the changes in velocity due to the electric field E:
  transport_v(f,dt);
}

void VlaFoo::solve(int itmax, real tmax, real tsave1, real tsave2)
{
  int it=0.0;
  int framenum=0.0;
  real tnow=0.0;

  //if(itmax < 0) itmax=INT_MAX;

  curves(tnow,true);   // erase contents of output files.
  curves(tnow);
  plot(framenum++);
  real dt0=dt;

  // FIXME: saving logic should deal more elegantly with multiple cases.
  bool savenow1=false;
  bool savenow2=false;
  real nextsave1=tnow+tsave1;
  real nextsave2=tnow+tsave2;

  std::cout << "t=" << tnow << std::endl;
  std::cout << "dt=" << dt << std::endl;
  std::cout << "time-stepping...\n" << std::endl;
  
  bool error=false;
  
  double wall0 = get_wall_time();
  double wallinterval=10.0; // update number of iterations every second
  int wallouts=0;
  bool go=true;

  int nout1=1;
  int nout2=1;

  double tsave=std::min(tsave1,tsave2);

  while(go) {
    it++;
    
    // Output 2D quantities:
    if(savenow2 && !error) {
      nout2++;
      plot(framenum++);
      nextsave2 = nout2 * tsave2;
      dt=dt0; // restore time-step
      savenow2=false;
    }

    // Time-step and increase time
    tnow += dt; // must be done before time_step, which may change dt.
    time_step(dt);
    dt=std::min(dt,tsave); // do not jump past tsave is we are dynamic.

    // Check for nans and infs.
    error=check_for_error();
    if(error) break;

    // Output scalar quantities:
    if(savenow1 && !error) {
      nout1++;
      curves(tnow);
      nextsave1 = nout1 * tsave1;
      dt = dt0; // restore time-step
      savenow1 = false;
    }

    if(tnow + dt > nextsave1) {
      dt0 = dt;
      dt = nextsave1 - tnow; // shorten dt so that we save at the right time
      savenow1 = true;
    }

    if(tnow + dt >= nextsave2) {
      if(! savenow1) dt0=dt;
      dt = nextsave2 - tnow; // shorten dt so that we save at the right time
      savenow2 = true;
    }

    // iteration limits:
    if(tnow >= tmax) {
      go = false;
      std::cout << "reached tmax=" << tmax << std::endl;
    }
    if(itmax > 0 && it >= itmax) {
      go=false;
      std::cout << "reached itmax="<<itmax << std::endl;
    }

    // diagnostic output:
    if(get_wall_time() - wall0 > wallinterval) {
      wall0=get_wall_time();
      wallouts++;
      std::cout << it << std::flush;
      if(wallouts % 10 == 0)
	std::cout << "\tt=" << tnow << "\tdt=" << dt << std::endl;
      else
	std::cout<<" ";
    }

  }
  if(!error) {
    curves(tnow);
    std::cout << "final time t=" << tnow << "\n" << std::endl;
  } else {
    std::cout << "Finished on error" << std::endl;
  }
}

void VlaFoo::curve(real tnow, real value, const char* fname,
		   bool clear_file)
{
  std::string outname;
  outname.append(outdir);
  outname.append("/");
  outname.append(fname);
  std::ofstream plot_file;
  if(clear_file) {
    plot_file.open(outname.c_str());
    plot_file.close();
  } else {
    plot_file.open(outname.c_str(),std::fstream::app);
    plot_file << tnow << "\t" << value << "\n";
    plot_file.close();
  }
}

void VlaFoo::curves(real tnow, bool clear_file)
{
  curve(tnow,compute_kin_energy(),"ekvt",clear_file);
  curve(tnow,compute_elec_energy(),"eEvt",clear_file);
  curve(tnow,fmin(),"fmin",clear_file);
  curve(tnow,fmax(),"fmax",clear_file);
  curve(tnow,fmax(),"ftot",clear_file);
}

void VlaFoo::plot(int framenum)
{
  // output the 2D field at a given time.  

  std::string outname;
  outname.append(outdir);
  outname.append("/");
  outname.append("plot");

  outname.append(string(framenum));

  xdr::oxstream plot_file;
  plot_file.open(outname.c_str());
  plot_file << nx << nv;
  for(int i=0; i < nx; i++)
    for(int j=0; j < nv; j++)
      plot_file << f[i][j];
  plot_file.close();
}

real VlaFoo::compute_kin_energy()
{
  real energy=0.0;
  for(int i=0; i < nx; i++) {
    for(int j=0; j < nv; j++) {
      real vj=index2v(j);
      energy += dv*dx*f[i][j]*vj*vj;
    }
  }
  return energy;
}

real VlaFoo::compute_elec_energy()
{
  compute_E(f);
  real energy=0.0;
  for(int i=0; i < nx; i++) {
    real Ei=E[i];
    energy += dx*Ei*Ei;
  }
  return sqrt(energy);
}


real VlaFoo::fmin()
{
  real val=f[0][0];
  
  for(int i=0; i < nx; ++i) {
    array1<real>::opt fi=f[i];
    for(int j=0; j < nv; ++j) {
      if(fi[j] < val) 
	val=fi[j];
    }
  }
  return val;
}

real VlaFoo::fmax()
{
  real val=f[0][0];
  
  for(int i=0; i < nx; ++i) {
    array1<real>::opt fi=f[i];
    for(int j=0; j < nv; ++j) {
      if(fi[j] > val) 
	val=fi[j];
    }
  }
  return val;
}

real VlaFoo::ftot()
{
  real val=0.0;
  
  for(int i=0; i < nx; ++i) {
    array1<real>::opt fi=f[i];
    for(int j=0; j < nv; ++j)
      val += fi[j];
  }
  return dx*dv*val;
}

bool VlaFoo::check_for_error()
{
  bool error=false;
  for(int i=0; i < nx; ++i) {
    real *fi=f[i];
    for(int j=0; j < nv; ++j) {
      real fij=fi[j];
      if(isnan(fij) || isinf(fij)) {
	error=true;
	break;
      }
    }
  }
  return error;
}


int main(int argc, char* argv[])
{
  std::cout << "command:" << std::endl;
  for(int i=0; i < argc; ++i) { 
    std::cout << argv[i];
    if(i < argc-1) 
      std::cout << " ";
    else
      std::cout << std::endl;
  }

  int nx;
  int nv;
  int itmax;
  real tmax;
  real cfl;
  real dt;
  real tsave1;
  real tsave2;
  real kx;
  real eps;
  real vmax;
  std::string ic;
  std::string outdir;
  std::string rk_name="euler";
  bool dynamic;
  double tolmin, tolmax;
  std::string config_file;

  // Set up and read the command line and config file.
  {
    // try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help,h", "produce help message")
      ("config,c",
       po::value<std::string>(&config_file)->default_value("test/vlafoo.cfg"),
       "name of a file of a configuration.")
      ;
    
    // Declare a group of options that will be allowed both on command
    // line and in config file
    po::options_description config("Configuration");
    config.add_options()
      ("nx", po::value<int>(&nx)->default_value(10),"nx")
      ("nv", po::value<int>(&nv)->default_value(10),"nv")
      ("tmax", po::value<double>(&tmax)->default_value(10.0),"tmax")
      ("vmax", po::value<double>(&vmax)->default_value(6.0),"vmax")
      ("tsave1", po::value<double>(&tsave1)->default_value(0.1),"tsave1")
      ("tsave2", po::value<double>(&tsave2)->default_value(1.0),"tsave2")
      ("dt", po::value<double>(&dt)->default_value(0.0),"dt")
      ("cfl", po::value<double>(&cfl)->default_value(0.4),"cfl")
      ("itmax", po::value<int>(&itmax)->default_value(1000000),"itmax")
      ("eps", po::value<double>(&eps)->default_value(5e-3),"eps")
      ("kx", po::value<double>(&kx)->default_value(0.2),"kx")
      ("ic", po::value<std::string>(&ic)->default_value("landau"),"ic")
      ("outdir", po::value<std::string>(&outdir)->default_value("test"),
       "outdir")
      ("rk_name", po::value<std::string>(&rk_name)->default_value("rk2"),
       "rk_name: euler or rk2")
      ("dynamic", po::value<bool>(&dynamic)->default_value(true),"dynamic")
      ("tolmin", po::value<double>(&tolmin)->default_value(0.0003),"tolmin")
      ("tolmax", po::value<double>(&tolmax)->default_value(0.0005),"tolmax")
      ;
    
    // Options for the command line:
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);
    
    // Options for the config file:
    po::options_description config_file_options;
    config_file_options.add(config);
    
    // The variables map
    po::variables_map vm;

    // Read the command-line options and store in the parsed version
    // in opts.
    po::parsed_options cl_opts = po::parse_command_line(argc,argv,
							cmdline_options);
    po::store(cl_opts,vm);
    po::notify(vm);    
    
    // Read the config file
    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
      // TODO: create a new file if the file doesn't already exist.
      std::cout << "can not open config file: " << config_file << std::endl;
      return 0;
    } else {
      po::parsed_options f_opts = parse_config_file(ifs, config_file_options);
      store(f_opts, vm);
      notify(vm);
      
      // Update the config file with the values from the command line.
      {
	// TODO: put this in a function.

	std::cout << "... finding command-line updates:" << std::endl;
	std::ifstream input(config_file.c_str());
	std::vector<std::string> lines;
	{
	  { // read each line and add to lines:
	    std::string line;
	    while(getline(input,line))
	      lines.push_back(line);
	  }
	  
	  // Search for updates from command-line:
	  for(vec_opt::iterator ipo = cl_opts.options.begin();
	      ipo != cl_opts.options.end(); 
	      ++ipo) {
	    po::basic_option<char>& l_option = *ipo;

	    // TODO check that the var name is allowed
	    //if(l_option.string_key != "config") {

	    std::cout << "command line:\t" 
		      << l_option.string_key 
		      << "=" 
		      << l_option.value[0] << std::endl;
	      
	    for(unsigned int i=0; i < lines.size(); ++i) {
	      std::size_t found = lines[i].find(l_option.string_key);
	      if(found == 0) {
		//std::cout << "\tfound " << l_option.string_key << std::endl;
		std::cout << config_file << ":\t" << lines[i] << std::endl;
	      }
	    }
	  }
	}
      }
      
      std::cout << "... finding missing variables in config file:" << std::endl;
      for(po::variables_map::iterator vit = vm.begin(); 
	  vit != vm.end(); 
	  ++vit) {
	std::string v_name = vit->first;
	//std::cout << v_name  << std::endl;
	
	bool found=false;
	
	for(vec_opt::iterator ipo = f_opts.options.begin();
	    ipo != f_opts.options.end();
	    ++ipo) {
	  po::basic_option<char>& l_option = *ipo;
	  //if(l_option.string_key != "config") {
	  //std::cout << l_option.string_key << std::endl;

	  if(v_name == l_option.string_key)
	    found=true;
	}
	if(!found) {
	  // NB: we do a lot of casting and catching exceptions so
	  // that we can cout boost::any.
	  std::cout << vit->first << "=";
	  try { std::cout << vit->second.as<double>() << std::endl;
	  } catch(...) {/* do nothing */ }
	  try { std::cout << vit->second.as<int>() << std::endl;
	  } catch(...) {/* do nothing */ }
	  try { std::cout << vit->second.as<std::string>() << std::endl;
	  } catch(...) {/* do nothing */ }
	  try { std::cout << vit->second.as<bool>() << std::endl;
	  } catch(...) {/* do nothing */ }
	}
      }
    }

    if (vm.count("help")) {
      std::cout << generic << std::endl;;
      return 0;
    }
  }
  
  VlaFoo vla(nx,nv,cfl,eps,kx,vmax,outdir,rk_name,dynamic,tolmin,tolmax);
  
  if(dt != 0.0) 
    vla.set_dt(dt);
    
  std::cout << "Setting up initial conditions..." << std::endl;
  vla.initial_conditions(ic);

  std::cout << "Solving..." << std::endl;
  vla.solve(itmax,tmax,tsave1,tsave2);
  
  std::cout << "Done." << std::endl;

  return 0;
}
