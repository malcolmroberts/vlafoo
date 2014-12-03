#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "Array.h"
#include "timestepper.hpp"
#include "complex"
#include "fftw++.h"
#include <fstream>
#include <string.h> // strcpy
#include <sys/time.h> // wall time
#include <sys/stat.h> //mkdir

typedef std::complex<double> Complex;
using namespace Array;
//typedef array1<double>::opt vector;

// A C++ class for solving Vlasov equation
// with a velocity Fourier transform 

class VlaFoo: public timestepper<double>{
private:
  int nvk, nxk, nq;
  double PI, dx, dv, k0x, k0v;

  // The Fourier transforms:
  fftwpp::rcfft1d *rc_v, *rc_x;
  fftwpp::crfft1d *cr_v, *cr_x;

  array2<double> f;
  array2<double> S; // Source term
  array1<double> E; // Electric field in space
  array1<double> rho; // Mean velocity (ie velocity density)
  array1<Complex> rhok; // Fourier-transformed version (for derivative)
  array1<Complex> vtempk; // 1D array for Fourier transformed v per each x.
  array1<double> vtemp; // 1D array for Fourier transformed v per each x.

  int nx, nv;

  double cfl, dt;
  double eps, Lx, kx; // Landau damping parameters 
  double vmax;
  std::string outdir;
  std::string restart_filename;
  bool did_restart;

  int framenum;
  double lastsave1, lastsave2;
  
  double index2v(int j) {return j*dv - vmax;}
  bool check_for_error();

  double interpolate(double x, int nq, double* x0, double*f);
  void transport_x(double &dtx);

  void compute_rho(array2<double> &f);// Compute the mean velocity
  void compute_E(array2<double> &f); // Compute the electric field
  // Compute the source term due to the electric field.
  void compute_dfdv(array2<double> &f, array2<double> &S);
  // Compute the source term due to the electric field.
  void compute_negEdfdv(array2<double> &f, array2<double> &S);
  // Step velocity by Fourier for the term E \grad_v f
  void transport_v(array2<double> &f, double &dtv);

  // Output:

  // functions to compute output quantities:
  double compute_kin_energy();
  double compute_elec_energy();
  double fmin();
  double fmax();
  double ftot();
  // scalar output:
  void curves(double tnow, bool clear_file=false);
  void curve(double tnow, double value, const char *fname, bool clear_file);
  // 2D output:
  void plot(int framenum, double tnow);

  void write_restart(const double tnow, const double dt);
  void read_restart(double & tnow, double & dt);

public:
    
  // Constructor (no default constructor for you!  Ha!)
  VlaFoo(int nx, int nv, double cfl, double eps, double Lx, double kx, 
	 double vmax, std::string &outdir, std::string &rk_name, 
	 double tolmin, double tolmax, double dtmax, bool redo): 
    nx(nx), nv(nv), cfl(cfl), eps(eps), Lx(Lx), kx(kx), vmax(vmax), 
    outdir(outdir){
    // Set parameters:
    
    // Compute derived parameters:
    PI = 4.0 * atan((double)1.0);
    k0v = PI / vmax;
    k0x = 2 * PI / Lx;
    dx = Lx / nx;
    dv = 2 * vmax / nv;
    dt = dx * cfl / vmax;
    nvk = nv / 2 + 1;
    nxk = nx / 2 + 1;
    nq = 4; // number of quadrature points
    
    // Allocate the time-stepper from the base class
    rk_allocate(nx * nv, rk_name, tolmin, tolmax, dtmax, redo);

    std::cout << "Creating new 1D Vlasov solver:" << std::endl;
    std::cout << "nx\t" << nx << std::endl;
    std::cout << "nv\t" << nv << std::endl;
    std::cout << "nvk\t" << nvk << std::endl;
    std::cout << "Lx\t" << Lx << std::endl;
    std::cout << "vmax\t" << vmax << std::endl;
    std::cout << "eps\t" << eps << std::endl;
    std::cout << "kx\t" << kx << std::endl;
    std::cout << "dx\t" << dx << std::endl;
    std::cout << "k0v\t" << k0v << std::endl;
    std::cout << "k0x\t" << k0x << std::endl;
    std::cout << "cfl\t" << cfl << std::endl;
    std::cout << "Output directory\t" << outdir << std::endl;
    std::cout << "Runge-Kutta method:\t" << rk_name << std::endl;
    std::cout << "tolmin:\t" << tolmin << std::endl;
    std::cout << "tolmax:\t" << tolmax << std::endl;

    // Allocate data arrays (NB: FFTs need double data to be
    // complex-aligned):
    size_t align = sizeof(Complex);
    f.Allocate(nx, nv, align);
    S.Allocate(nx, nv, align);
    E.Allocate(nx, align);
    vtempk.Allocate(nvk, align);
    vtemp.Allocate(nv, align);
    rho.Allocate(nx, align);
    rhok.Allocate(nxk, align);

    rc_v = new fftwpp::rcfft1d(nv, f[0], vtempk());
    cr_v = new fftwpp::crfft1d(nv, vtempk(), f[0]);
    rc_x = new fftwpp::rcfft1d(nx, rho, rhok);
    cr_x = new fftwpp::crfft1d(nx, rhok, rho);

    restart_filename = "restart";

    int mkdirnotok = mkdir(outdir.c_str(),S_IRWXU);
    if(mkdirnotok) {
      // TODO: this claims to fail if directory is already present.
      // std::cout << "ERROR: mkdir failed!" << std::endl;
    }
    
  }

  // Destructor
  ~VlaFoo() {
    delete rc_x;
    delete cr_x;
    delete rc_v;
    delete cr_v;
 
    f.Deallocate();
    S.Deallocate();
    E.Deallocate();
    vtemp.Deallocate();
    vtempk.Deallocate();
    rho.Deallocate();
    rhok.Deallocate();
  }

  void set_dt(double new_dt) {dt = new_dt;}

 // Set initial conditions
  void initial_conditions(std::string &ic, double &tnow, double &dt, 
			  bool restart);

  void rk_source(double *f, double *S);

  // Full resolution by Strang splitting
  void time_step(double &dt);
  void solve(double tnow, int itmax, double tmax, double tsave1, double tsave2);

  void write_stats(int it, double tnow, double dt, bool reset=false);
};

#include <sstream>
template <typename T>
std::string string(T Number)
{
  std::ostringstream ss;
  ss << Number;
  return ss.str();
}

// return a mod b
int mod(int a, int b) {
  while(a < 0) a += b;
  return a % b;
}

#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
// double get_cpu_time(){
//     return (double)clock() / CLOCKS_PER_SEC;
// }

