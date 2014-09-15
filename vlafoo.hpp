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

typedef double real;
typedef std::complex<real> Complex;
using namespace Array;
//typedef array1<real>::opt vector;

// A C++ class for solving Vlasov equation
// with a velocity Fourier transform 

class VlaFoo: public timestepper<real>{
private:
  int nvk, nxk, nq;
  real PI, dx, dv, k0x, k0v;

  // The Fourier transforms:
  fftwpp::rcfft1d *rc_v, *rc_x;
  fftwpp::crfft1d *cr_v, *cr_x;

  array2<real> f;
  array2<real> S; // Source term
  array1<real> E; // Electric field in space
  array1<real> rho; // Mean velocity (ie velocity density)
  array1<Complex> rhok; // Fourier-transformed version (for derivative)
  array1<Complex> vtempk; // 1D array for Fourier transformed v per each x.
  array1<real> vtemp; // 1D array for Fourier transformed v per each x.

  int nx, nv;

  real cfl, dt;
  real eps, kx; // Landau damping parameters 
  real L, vmax;
  std::string outdir;
  std::string restart_filename;
  bool restart;

  int framenum;
  real lastsave1, lastsave2;
  
  real index2v(int j) {return j*dv - vmax;}
  bool check_for_error();

  real interpolate(real x, int nq, real* x0, real*f);
  void transport_x(real &dtx);

  void compute_rho(array2<real> &f);// Compute the mean velocity
  void compute_E(array2<real> &f); // Compute the electric field
  // Compute the source term due to the electric field.
  void compute_dfdv(array2<real> &f, array2<real> &S);
  // Compute the source term due to the electric field.
  void compute_negEdfdv(array2<real> &f, array2<real> &S);
  // Step velocity by Fourier for the term E \grad_v f
  void transport_v(array2<real> &f, real &dtv);

  // Output:

  // functions to compute output quantities:
  real compute_kin_energy();
  real compute_elec_energy();
  real fmin();
  real fmax();
  real ftot();
  // scalar output:
  void curves(real tnow, bool clear_file=false);
  void curve(real tnow, real value, const char *fname, bool clear_file);
  // 2D output:
  void plot(int framenum, double tnow);

  void write_restart(const double tnow, const double dt);
  void read_restart(double & tnow, double & dt);

public:
    
  // Constructor (no default constructor for you!  Ha!)
  VlaFoo(int nx, int nv, real cfl, real eps, real kx, real vmax, 
	 std::string &outdir, std::string &rk_name, bool dynamic, 
	 double tolmin, double tolmax): 
    nx(nx),nv(nv),cfl(cfl),eps(eps),kx(kx),vmax(vmax),outdir(outdir){
    // Set parameters:
    
    // Compute derived parameters:
    PI=4.0*atan((real)1.0);
    L=2.0*PI/kx;
    k0v=PI/vmax;
    k0x=2*PI/L;
    dx=L/nx;
    dv=2*vmax/nv;
    dt=dx*cfl/vmax;
    nvk=nv/2+1;
    nxk=nx/2+1;
    nq=4; // number of quadrature points
    
    // Allocate the time-stepper from the base class
    rk_allocate(nx*nv,rk_name,dynamic,tolmin,tolmax);

    std::cout << "Creating new 1D Vlasov solver:" << std::endl;
    std::cout << "nx\t" << nx << std::endl;
    std::cout << "nv\t" << nv << std::endl;
    std::cout << "nvk\t" << nvk << std::endl;
    std::cout << "L\t" << L << std::endl;
    std::cout << "vmax\t" << vmax << std::endl;
    std::cout << "eps\t" << eps << std::endl;
    std::cout << "kx\t" << kx << std::endl;
    std::cout << "dx\t" << dx << std::endl;
    std::cout << "k0v\t" << k0v << std::endl;
    std::cout << "k0x\t" << k0x << std::endl;
    std::cout << "cfl\t" << cfl << std::endl;
    std::cout << "Output directory\t" << outdir << std::endl;
    std::cout << "Runge-Kutta method:\t" << rk_name << std::endl;
    std::cout << "dynamic?:\t" << dynamic << std::endl;
    std::cout << "tolmin:\t" << tolmin << std::endl;
    std::cout << "tolmax:\t" << tolmax << std::endl;

    // Allocate data arrays (NB: FFTs need real data to be
    // complex-aligned):
    size_t align=sizeof(Complex);
    f.Allocate(nx,nv,align);
    S.Allocate(nx,nv,align);
    E.Allocate(nx,align);
    vtempk.Allocate(nvk,align);
    vtemp.Allocate(nv,align);
    rho.Allocate(nx,align);
    rhok.Allocate(nxk,align);
    nvk=nv/2+1;
    rc_v=new fftwpp::rcfft1d(nv,f[0],vtempk());
    cr_v=new fftwpp::crfft1d(nv,vtempk(),f[0]);
    rc_x=new fftwpp::rcfft1d(nx,rho,rhok);
    cr_x=new fftwpp::crfft1d(nx,rhok,rho);

    restart_filename="restart";

    int mkdirnotok=mkdir(outdir.c_str(),S_IRWXU);
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

  void set_dt(real new_dt) {dt = new_dt;}

 // Set initial conditions
  void initial_conditions(std::string &ic, double &tnow, double &dt);

  void rk_source(real *f, real *S);

  // Full resolution by Strang splitting
  void time_step(real &dt);
  void solve(double tnow, int itmax, real tmax, real tsave1, real tsave2);
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

