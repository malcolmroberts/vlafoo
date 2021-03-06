#include "vlafoo.hpp"
#include "xstream.h"
#include <fstream>
#include "clopts.hpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options;


// Compute initial condition
void VlaFoo::initial_conditions(std::string &ic, double v0, 
				double &tnow, double &dt,
				bool restart) {
  double overs2pi = 1.0 / sqrt(2.0 * PI);

  f.Load(0.0);
  framenum = 0;
  lastsave1 = 0.0;
  lastsave2 = 0.0;
  did_restart = false;

  // for(int ix=0;ix<nx;ix++) f[ix][1]=1.0;

  if(restart) {
    read_restart(tnow, dt);
    did_restart = true;
    return;
  }

  if(ic == "test0") {
    f[0][nv / 4] = 1.0;
    f[0][3 * nv / 4] = 1.0;
    return;
  }
  if(ic == "test1") {
    for(int j=0; j < nv; ++j)
      f[0][j] = 1.0;
    return;
  }

  if(ic == "landau") {
    std::cout << "Landau initial conditions." << std::endl;

    // landau damping
    // f0=exp(-vt*vt/2)*(1+eps*cos(kx*xt))/sqrt(2*pi)
    for(int i = 0; i < nx; i++){
      double x = i * dx;
      double xterm = (1.0 + eps * cos(kx * x)) * overs2pi;
      for(int j = 0; j < nv; j++){
	double v = index2v(j);
	f[i][j] = exp(-0.5 * v * v) * xterm;
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
    //double v0 = 3.0;
    double overs2pi = 0.5 / sqrt(2.0 * PI);
    for(int i = 0; i < nx; i++){
      double x = i * dx;
      double xterm = (1.0 + eps * cos(kx * x)) * overs2pi;
      for(int j = 0; j < nv; j++){
	double v = index2v(j);
	double sum = v + v0;
	double diff = v - v0;
	f[i][j] = xterm * (exp(-0.5 * sum * sum) + exp(-0.5 * diff * diff));
      }
    }
    return;
  }
  
  std::cout << "Abort: no initial conditions found matching " << ic 
	    << std::endl;
  exit(1);
}

double VlaFoo::interpolate(double x, int nq, double* x0, double* y0)
{
  // Lagrange polynomial interpolation
  double val = 0;
  if(nq > 0) {
    // NB: this is probably not the most numerically robust way to do
    // this.
    for(int j = 0; j < nq; ++j) {
      double lj = 1.0;
      double xj = x0[j];
      for(int m = 0; m < j; ++m)
	lj *= (x - x0[m]) / (xj - x0[m]);
      for(int m = j + 1; m < nq; ++m)
	lj *= (x - x0[m]) / (xj - x0[m]);
      val += y0[j] * lj;
      
    }
  } else {
    std::cout << "Error: spline of order "<< nq << " not implemented."
	      << std::endl;
    exit(1);
  }

  return val;
}

void VlaFoo::transport_x(double &dtx)
{
  // Advance f in time, taking accont the A\partial_x f term, with
  // A=diag[v].

  // The velocity for the particles is given by the distribution, ie
  // particles in the array at (i,j) have position x_i and velocity
  // v(j).

  double *x0 = new double[nq];
  double *f0 = new double[nq];

  double overdx = 1.0 / dx;

  // f(x,t+dt) = f(x-v*dt,t)
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < nv; j++) {
      double x = i * dx; // current position
      double v = index2v(j);
      if(x < 0.0) x += Lx;
      if(x > Lx) x -= Lx;

      double xold=x-v*dtx; // previous position
      if(xold < 0.0) xold += Lx;
      if(xold > Lx) xold -= Lx;
      
      // Find the neighbouring points:
      int iold = (int)floor(xold * overdx);
      int nleft = nq / 2 - 1;
      int ileft = iold - nleft;
      if(ileft < 0) ileft += nx;
      x0[0] = dx * ileft;
      for(int q = 1; q < nq; q++)
	x0[q] = x0[0] + q * dx;

      for(int q = 0; q < nq; q++) {
	if(abs(x0[q] - xold) > 0.5 * Lx)
	  x0[q] -= Lx;
      }

      // Find the f-values for the interpolation:
      for(int q = 0; q < nq; q++) {
	//int qi=mod(nleft+i,nx);
	int qi = (ileft + q) % nx;

	f0[q]=f[qi][j]; // NB: periodic in x.
      }

      // Put the updated f in S:
      double val = interpolate(xold, nq, x0, f0);
      S[i][j] = val;
    }
  }
  delete[] x0;
  delete[] f0;
  
  // TODO: optimise into a swap?
  for(int i = 0; i < nx; ++i) {
    array1<double>::opt fi = f[i];
    array1<double>::opt Si = S[i];
    for(int j = 0; j < nv; ++j)
      fi[j] = Si[j];
  }
}

void VlaFoo::compute_rho(array2<double> &f)
{  
  // \rho = \int \approx  dv \sum_j f_j(x,t)
  for(int i = 0; i < nx; i++) {
    double rhoi = 0.0;
    array1<double>::opt fi = f[i];
    for(int j = 0; j < nv; j++)
      rhoi += fi[j];
    rho[i] = rhoi * dv;
  }
}

void VlaFoo::compute_E(array2<double> &f)
{
  // Compute
  // E(x) = \int dx \int_{-vmax}^{vmax} dv f(x,v)
  // from
  // \del\cdot E = -1 + \rho
  // via a FFT.
  
  compute_rho(f);

  // for(int i = 0; i < nx; i++) 
  //   rho[i] -= 1.0;

#if 1
  // using Fourier integration

  rc_x->fft(rho,rhok);

  // Set mean to zero:
  rhok[0] = 0.0;
  double overnx = 1.0 / nx;
  Complex Iovernx = Complex(0, overnx);
  int stop = nxk - 1;
  for(int i = 1; i < stop; i++) {
    double k0xi = k0x * i;
    rhok[i] *= -Iovernx / k0xi;
  }
  rhok[nxk-1] = 0.0; // kill the nyquist

  cr_x->fft(rhok, E);
#else
  // using crazy cotan stuff
  
  rc_x->fft(rho,rhok);

  // Set mean to zero:
  rhok[0]=0.0;
  double overnx = 1.0 / nx;
  for(int i = 1; i < nxk - 1; i++) {
    double k = k0x * i;
    Complex rhoki = rhok[i];
    double re = rhoki.real();
    double im = rhoki.imag();
    double ct = 1.0 / tan(PI * k * overnx);
    rhok[i] = 0.5 * k * dx * ct * Complex(im, -re);
  }
  rhok[nxk - 1] = 0.0; // kill the nyquist

#endif
}

void VlaFoo::compute_dfdv(array2<double> &f, array2<double> &S)
{
  // compute \partial_v f_i via Fourier differentiation.
  
  double overnv=1.0/nv; // FFT normalization term.
 
  Complex I=Complex(0,1);

  for(int i=0; i<nx; i++) {
    // copy velocity for each i into a temp buffer (NB double to complex
    // FFT overwrites input):
    array1<double>::opt fi = f[i];
    array1<double>::opt Si = S[i];

    for(int j = 0; j < nv; ++j)
      vtemp[j] = fi[j];

    rc_v->fft(vtemp, vtempk);
    unsigned int stop = nvk - 1;
    for(unsigned j = 0; j < stop; j++) {
      double k = j * k0v;
      vtempk[j] *= I * k * overnv;
    }
    vtempk[stop] = 0.0;  // Set Nyquist to zero.
    cr_v->fft(vtempk, Si); // transform back and put into source term
  }
}


// Compute -E_i \partial_v f_i
void VlaFoo::compute_negEdfdv(array2<double> &f, array2<double> &S)
{
  compute_E(f);
  compute_dfdv(f, S);

  for(int i = 0; i < nx; ++i) {
    array1<double>::opt Si = S[i];
    double Ei = E[i];
    for(int j = 0; j < nv; ++j)
      Si[j] *= -Ei;
  }
}

void VlaFoo::rk_source(double *f0, double *S0)
{
  // wrapper for compute_negEdfdv.

  array2<double> f(nx, nv, f0);
  array2<double> S(nx, nv, S0);
  compute_negEdfdv(f,S);
}

// Time step dtv of velocity resolution by Fourier
void VlaFoo::transport_v(array2<double> &f, double &dtv)
{
  // source term for f_i is -E_i \partial_v f_i and is called via
  // VlaFoo::rk_source.

  rk_step(f(),dtv);
}

void VlaFoo::time_step(double &dt)
{
  // Time-step by Strang splitting.

  // NB: transport_v may redo time-steps (adjusting dt in the process)
  // and must therefore be called first.

  // Compute the changes in velocity due to the electric field E:
  transport_v(f,dt);
 
  // Advection term from rho (mean velocity) in space:
  transport_x(dt);
}

void VlaFoo::solve(double tnow, int itmax, double tmax, 
		   double tsave1, double tsave2)
{
  int it = 0;

  if(!did_restart) {
    // erase contents of output files and add initial value
    curves(tnow, true); 
    curves(tnow);
    plot(framenum++, tnow);
    write_stats(it, tnow, dt, true);
  }

  double dt0 = dt;

  // TODO: saving logic should deal more elegantly with multiple cases.
  double nextsave1 = tnow + tsave1;
  double nextsave2 = tnow + tsave2;

  double nextsave = std::min(nextsave1, nextsave2);
  bool savenow = false;

  bool error = false;
  
  double wall0 = get_wall_time();
  double wallinterval = 10;// 10  // update number of iterations every second
  unsigned int checkinterval = 10; //10
  int wallouts = 0;
  bool go = true;

  double tsave = std::min(tsave1, tsave2);
  dt = std::min(dt, tsave); // do not jump past tsave is we are dynamic.
  dt0 = dt;

  std::cout << "t=" << tnow << std::endl;
  std::cout << "tmax=" << tmax << std::endl;
  std::cout << "dt=" << dt << std::endl;
  std::cout << "time-stepping...\n" << std::endl;
  
  while(go) {
    it++;
    
    /*
    std::cout << "\nit" << it << std::endl; // FIXME: temp
    std::cout << "t=" << tnow << std::endl; // FIXME: temp
    std::cout << "dt=" << dt << std::endl; // FIXME: temp
    std::cout << "t+dt=" << tnow+dt << std::endl; // FIXME: temp
    std::cout << "nextsave=" << nextsave << std::endl; // FIXME: temp
    */ 

    // Do not overshoot tmax
    if(tnow >= tmax)
      dt = tmax-tnow;

    // Time-step and increase time
    tnow += dt; // must be done before time_step, which may change dt.
    time_step(dt);
    dt = std::min(dt, tsave); // do not jump past tsave is we are dynamic.

    // Check for nans and infs.
    error = check_for_error();
    if(error) 
      break;

    // Save the fields if it's time to do so
    if(savenow) {
      if(tnow >= nextsave2) {
	lastsave2 = tnow;
	nextsave2 += tsave2;
	plot(framenum++, tnow);
      }
      if(tnow >= nextsave1) {
	lastsave1 = tnow;
	nextsave1 += tsave1;
	curves(tnow);
      }
      nextsave = std::min(nextsave1, nextsave2);
      dt = dt0; // restore time-step
      savenow = false;
    }
        
    // Adjust time-step to reach next save-point:
    if(tnow + dt >= nextsave) {
      dt0 = dt;
      dt = nextsave - tnow; 
      //std::cout << "\t" << dt << std::endl;
      savenow = true;
    }

    // iteration limits:
    if(tnow >= tmax) {
      dt = tmax - tnow;
      go = false;
      std::cout << "reached tmax=" << tmax << std::endl;
    }
    if(itmax > 0 && it >= itmax) {
      go = false;
      std::cout << "reached itmax="<< itmax << std::endl;
    }

    // diagnostic output:
    if(get_wall_time() - wall0 > wallinterval) {
      wall0 = get_wall_time();
      wallouts++;
      std::cout << it << std::flush;
      write_stats(it, tnow, dt);
      if(wallouts % checkinterval == 0) {
	write_restart(tnow, dt);
	std::cout << "\tt=" << tnow << "\tdt=" << dt << std::endl;
      } else {
	std::cout<<" ";
      }
    }

  }
  if(!error) {
    curves(tnow);
    write_stats(it, tnow, dt0);
    write_restart(tnow, dt0);
    std::cout << "final time t=" << tnow << "\n" << std::endl;
  } else {
    std::cout << "Finished on error" << std::endl;
  }
}

void VlaFoo::curve(double tnow, double value, const char* fname,
		   bool clear_file)
{
  std::string outname;
  outname.append(outdir);
  outname.append("/");
  outname.append(fname);
  std::ofstream plot_file;
  if(clear_file) {
    plot_file.open(outname.c_str(), std::ofstream::out | std::ofstream::trunc);
    plot_file.close();
  } else {
    plot_file.open(outname.c_str(),std::fstream::app);
    plot_file << tnow << "\t" << value << std::endl;
    plot_file.close();
  }
}

void VlaFoo::curves(double tnow, bool clear_file)
{
  curve(tnow, compute_kin_energy(), "ekvt", clear_file);
  curve(tnow, compute_elec_energy(), "eEvt", clear_file);
  curve(tnow, fmin(), "fminvt", clear_file);
  curve(tnow, fmax(), "fmaxvt", clear_file);
  curve(tnow, ftot(), "ftotvt", clear_file);
  curve(tnow, fL1(), "fL1", clear_file);
  curve(tnow, fL2(), "fL2", clear_file);
}

void VlaFoo::plot(int framenum, double tnow)
{
  // output the 2D field at a given time.  

  std::string outname;
  outname.append(outdir);
  outname.append("/");
  outname.append("plot");

  outname.append(string(framenum));

  xdr::oxstream plot_file;
  plot_file.open(outname.c_str());
  plot_file << tnow;
  plot_file << nx << nv;
  for(int i = 0; i < nx; i++)
    for(int j = 0; j < nv; j++)
      plot_file << f[i][j];
  plot_file.close();
}

void VlaFoo::write_restart(const double tnow, const double dt)
{
  // Write to the restart file
  xdr::oxstream xout;
  std::string path_file = outdir + "/" + restart_filename;
  xout.open(path_file.c_str());
  
  xout << tnow << dt << framenum << lastsave1 << lastsave2;
  for(int i = 0; i < nx; i++)
    for(int j = 0; j < nv; j++)
      xout << f[i][j];
  xout.close();
}

void VlaFoo::read_restart(double &tnow, double & dt)
{
  xdr::ixstream xin;
  std::string path_file = outdir + "/" + restart_filename;
  xin.open(path_file.c_str());
  if(xin) {  
    xin >> tnow;
    xin >> dt;
    xin >> framenum;
    xin >> lastsave1;
    xin >> lastsave2;
    for(int i = 0; i < nx; i++)
      for(int j = 0; j < nv; j++)
	xin >> f[i][j];
    xin.close();
  } else {
    std::cout << "restart file "
	      << path_file
	      << " does not exist"
	      << std::endl;
    exit(1);
  }
}

void VlaFoo::write_stats(int it, double tnow, double dt, bool reset)
{
  std::string outname;
  outname.append(outdir);
  outname.append("/stats");
  
  std::ofstream plot_file;
  if(reset) {
    plot_file.open(outname.c_str(), std::ofstream::out | std::ofstream::trunc);
    plot_file << "#"
	      << "it" << "\t" 
	      << "tnow" << "\t"
	      << "dt" 
	      << std::endl;
  } else {
    plot_file.open(outname.c_str(),std::fstream::app);
  }
  plot_file << it << "\t" 
	    << tnow << "\t"
	    << dt 
	    << std::endl;
  plot_file.close();
  
}

double VlaFoo::compute_kin_energy()
{
  double energy=0.0;
  for(int i=0; i < nx; i++) {
    for(int j=0; j < nv; j++) {
      double vj=index2v(j);
      energy += dv*dx*f[i][j]*vj*vj;
    }
  }
  return energy;
}

double VlaFoo::compute_elec_energy()
{
  compute_E(f);
  double energy = 0.0;
  for(int i = 0; i < nx; i++) {
    double Ei = E[i];
    energy += dx * Ei * Ei;
  }
  return 0.5 * sqrt(energy);
}


double VlaFoo::fmin()
{
  double val=f[0][0];
  
  for(int i = 0; i < nx; ++i) {
    array1<double>::opt fi = f[i];
    for(int j = 0; j < nv; ++j) {
      if(fi[j] < val) 
	val = fi[j];
    }
  }
  return val;
}

double VlaFoo::fmax()
{
  double val = f[0][0];
  
  for(int i = 0; i < nx; ++i) {
    array1<double>::opt fi = f[i];
    for(int j = 0; j < nv; ++j) {
      if(fi[j] > val) 
	val = fi[j];
    }
  }
  return val;
}

double VlaFoo::ftot()
{
  double val = 0.0;
  
  for(int i = 0; i < nx; ++i) {
    array1<double>::opt fi = f[i];
    for(int j = 0; j < nv; ++j)
      val += fi[j];
  }
  return dx * dv * val;
}

double VlaFoo::fL1()
{
  double val = 0.0;
  
  for(int i = 0; i < nx; ++i) {
    array1<double>::opt fi = f[i];
    for(int j = 0; j < nv; ++j)
      val += abs(fi[j]);
  }
  return dx * dv * val;
}

double VlaFoo::fL2()
{
  double val = 0.0;
  
  for(int i = 0; i < nx; ++i) {
    array1<double>::opt fi = f[i];
    for(int j = 0; j < nv; ++j)
      val += fi[j] * fi[j];
  }
  return dx * dv * val;
}


bool VlaFoo::check_for_error()
{
  bool error = false;
  for(int i = 0; i < nx; ++i) {
    double *fi = f[i];
    for(int j = 0; j < nv; ++j) {
      double fij = fi[j];
      if(isnan(fij) || isinf(fij)) {
	error = true;
	break;
      }
    }
  }
  return error;
}


int main(int argc, char* argv[])
{
  std::cout << "Welcome to VlaFoo.\n"
	    <<"This is VlaFoo, Welcome to VlaFoo.\n" 
	    << "\tAnything is possible at VlaFoo...\n" 
	    <<"\t...the only limit is yourself.\n" 
	    <<"\tWelcome to VlaFoo." 
	    << std::endl;

  std::cout << "git branch: " << GITBRANCH
	    << "\tgit commit version: " << GITHASH 
	    << std::endl;

  std::cout << "command:" << std::endl;
  for(int i = 0; i < argc; ++i) { 
    std::cout << argv[i];
    if(i < argc-1) 
      std::cout << " ";
    else
      std::cout << std::endl;
  }

  int nx;
  int nv;
  int itmax;
  double tmax;
  double cfl;
  double dt, dtmax;
  double tsave1;
  double tsave2;
  double Lx;
  double kx;
  double eps;
  double vmax;
  std::string ic;
  double v0;
  std::string outdir, config_file;
  std::string rk_name;
  bool restart;
  double tolmin, tolmax;
  bool redo;
  std::string config_dir_file;

  // option maps
  po::options_description generic("Allowed options");
  po::options_description config("Configuration");
  po::options_description cmdline_options;
  po::options_description config_file_options;
  po::variables_map vm;

  // Set up and read the command line and config file.
  {
    generic.add_options()
      ("help,h", "produce help message")
      ("config,c",
       po::value<std::string>(&config_file)->default_value("vlafoo.cfg"),
       "name of a file of a configuration.")
      ("outdir", po::value<std::string>(&outdir)->default_value("test"),
       "outdir")
      ;
    
    // Declare a group of options that will be allowed both on command
    // line and in config file
    config.add_options()
      ("nx", po::value<int>(&nx)->default_value(10), "nx")
      ("nv", po::value<int>(&nv)->default_value(10), "nv")
      ("tmax", po::value<double>(&tmax)->default_value(10.0), "tmax")
      ("vmax", po::value<double>(&vmax)->default_value(10.0), "vmax")
      ("tsave1", po::value<double>(&tsave1)->default_value(0.1), "tsave1")
      ("tsave2", po::value<double>(&tsave2)->default_value(1.0), "tsave2")
      ("dt", po::value<double>(&dt)->default_value(0.0), "dt")
      ("dtmax", po::value<double>(&dtmax)->default_value(0.0), "dtmax")
      ("cfl", po::value<double>(&cfl)->default_value(0.4), "cfl")
      ("itmax", po::value<int>(&itmax)->default_value(1000000), "itmax")
      ("eps", po::value<double>(&eps)->default_value(5e-3), "eps")
      ("Lx", po::value<double>(&Lx)->default_value(31.4159265358979), "Lx")
      ("kx", po::value<double>(&kx)->default_value(0.2), "kx")
      ("ic", po::value<std::string>(&ic)->default_value("landau"), 
       "initial conditions: landau, doublestream, test0, test1")
      ("v0", po::value<double>(&v0)->default_value(3.0), "v0")
      ("rk_name", po::value<std::string>(&rk_name)->default_value("rk2"),
       "rk_name: euler, rk2, or rk2d (dynamic)")
      ("tolmin", po::value<double>(&tolmin)->default_value(0.0003), "tolmin")
      ("tolmax", po::value<double>(&tolmax)->default_value(0.0005), "tolmax")
      ("restart", po::value<bool>(&restart)->default_value(false), "restart")
      ("redo", po::value<bool>(&redo)->default_value(false), "redo")
      ;

    // Options for the command line:
    cmdline_options.add(generic).add(config);
    // Options for the config file:
    config_file_options.add(config);
  }

  {
    // Read the command-line options and store in the parsed version
    // in opts.

    // Keep a copy for to update the config file:
    po::parsed_options cl_opts 
      = po::parse_command_line(argc, argv, cmdline_options);

    // positional arguments (outdir):
    po::positional_options_description p;
    p.add("outdir", -1);
    store(po::command_line_parser(argc, argv).
	  options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << cmdline_options << std::endl;
      return 0;
    }
    std::cout << "Finished with command-line arguments." << std::endl;

    // Read the config file
    {
      config_dir_file = outdir + "/" + config_file;
      std::ifstream ifs(config_dir_file.c_str());
      if (!ifs) {
	std::cout << "Creating empty config file " 
		  << config_dir_file 
		  << std::endl;
	
	make_empty_file(outdir,config_file);
	
	ifs.open(config_dir_file.c_str()); 
      }
      
      po::parsed_options f_opts = parse_config_file(ifs, config_file_options);
      store(f_opts, vm);
      notify(vm);
      ifs.close();
    } 

    // Add the any missing entries to the config file and then update
    // with command-line arguments.
    fill_config_file(config_dir_file, vm);
    update_config_file(config_dir_file, cl_opts, vm);
    
    // Show all of the parameters in the variable map:
    std::cout << "Parameters used in this simulation:" << std::endl;
    show_vm(vm);

    make_asy_header(vm);
    make_asy_input(outdir, vm);
  }
  
  VlaFoo vla(nx, nv, cfl, eps, Lx, kx, vmax, outdir, rk_name, tolmin, 
	     tolmax, dtmax, redo);
  
  if(dt != 0.0) 
    vla.set_dt(dt);
    
  std::cout << "Setting up initial conditions..." << std::endl;
  double tnow = 0.0;
  vla.initial_conditions(ic, v0, tnow,dt, restart);

  std::cout << "Solving..." << std::endl;
  vla.solve(tnow, itmax, tmax, tsave1, tsave2);
  
  std::cout << "Done." << std::endl;

  return 0;
}
