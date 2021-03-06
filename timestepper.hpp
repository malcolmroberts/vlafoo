#include <iostream>
#include <stdlib.h>
#include <assert.h>  

template<class T>
class timestepper
{
private:
  T **S; // array of source buffers
  T *f_save;
  int rk_stages;
  bool f_save_allocated;
  bool fold_allocated;
  enum RKTYPE{EULER,RK2,RK2D};
  RKTYPE rk;
  double tolmin, tolmax;
  double dtmax;
protected:
  unsigned int rk_n; // number of data points
  bool redo;
  
public:
  timestepper()
  {
    f_save_allocated=false;
  }

  T max(T a, T b) {
    if(a > b)
      return a;
    return b;
  }

  T abs(T a) {
    if(a < 0)
      return -a;
    return a;
  }

  void rk_allocate(int n, std::string &rk_name, 
		   double tolmin0, double tolmax0,
		   double dtmax0,
		   bool redo0) {
    redo = redo0;
    dtmax = dtmax0;
    rk_n = n;
    rk_stages = 0;
    tolmin = tolmin0;
    tolmax = tolmax0;

    f_save_allocated = false;

    assert(tolmin < tolmax);

    if(rk_name == "euler") {
      rk_stages = 1;
      rk = EULER;
    }

    if(rk_name == "rk2") {
      rk_stages = 2;
      rk = RK2;
    }

    if(rk_name == "rk2d") {
      rk_stages = 2;
      rk = RK2D;
      f_save_allocated = true;
      f_save = new T[rk_n];
    }

    if(rk_stages == 0) {
      std::cerr << "ERROR: unknown integrator " 
		<< rk_name << std::endl;
      exit(1);
    }

    S = new T*[rk_stages];
    for(int i = 0; i < rk_stages; ++i)
      S[i] = new T[rk_n];
  }

  ~timestepper() {
    for(int i = 0; i < rk_stages; i++)
      delete S[i];
    delete[] S;
    if(f_save_allocated) 
      delete[] f_save;
  }

  virtual void rk_source(T *f, T *S){}
  
  void rk_step(T* f, double &dt) {
    switch(rk) {
    case EULER:
      euler_step(f, dt);
      break;
    case RK2:
      rk2_step(f, dt);
      break;
    case RK2D:
      rk2d_step(f, dt);
      break;
    default:
      std::cerr << "ERROR: unknown integrator in rk_step." << std::endl;
      exit(1);
      break;
    }
  }

  T compute_error(T *s0, T *s1, T *f, T *f_save, double dt) {
    const T eps = 1e-15; // Avoids division by zero, implies that values
		       // below eps are not reliable.

    T error = 0.0;
    for(unsigned int i = 0; i < rk_n; ++i) {
      T diff 
	= dt * abs(s0[i] - s1[i]) / (0.5 * (abs(f[i]) + abs(f_save[i])) + eps);
      if(diff > error)
	error = diff;
    }
    return error;
  }

  bool all_finite(T *f) {
    bool isfinite = true;
    for(unsigned int i = 0; i < rk_n; i++) {
      if(isnan(f[i]) || isinf(f[i])) {
	isfinite = false;
	i = rk_n;
      }
    }
    return isfinite;
  }

  bool step_success(T error, bool finite) {
    if(!finite) 
      return false;
    
    if(!redo)
      return true;

    if(error > tolmax)  
      return false;
    
    return true;
  }

  void adjust_dt(T error, bool finite, double &dt) {
    if(!finite) {
      dt *= 0.7;
      return;
    }

    if(error > tolmax) {
      dt *= 0.7;
      return;
    }

    if(error < tolmin)
      dt *= 1.2;
    
    if(dtmax > 0 && dt > dtmax)
      dt=dtmax;

    return;
  }

  void euler_step(T* f, const double dt) {
    rk_source(f, S[0]);
    T *S0 = S[0];
    for(unsigned int i = 0; i < rk_n; ++i)
      f[i] += dt*S0[i];
  }

  void rk2_step(T* f, double &dt) {
    T* s;

    s = S[0];
    rk_source(f, s);
    double halfdt = 0.5 * dt;
    for(unsigned int i = 0; i < rk_n; i++)
      f[i] += halfdt * s[i];
    
    s = S[1];
    rk_source(f,s);
    for(unsigned int i = 0; i < rk_n; ++i)
      f[i] += dt*s[i];
  }

  void rk2d_step(T* f, double &dt)
  {
    T* s;

    T error = 0.0;

    bool reset=false;
    // backup the data in case we need to redo the step.
    for(unsigned int i = 0; i < rk_n; ++i) {
      f_save[i] = f[i]; 
    }

    bool done=false;
    while(!done) {

      if(reset) {
	for(unsigned int i = 0; i < rk_n; ++i) {
	  f[i] = f_save[i];
	}
      }
      reset = true;

      // First stage
      s = S[0];
      rk_source(f, s);
      double halfdt = 0.5 * dt;
      for(unsigned int i = 0; i < rk_n; ++i)
	f[i] += halfdt * s[i];
      
      // Second stage
      s = S[1];
      rk_source(f, s);
      for(unsigned int i = 0; i < rk_n; ++i)
	f[i] += dt * s[i];
     
      // check for error and react:
      bool finite = all_finite(f);
      if(finite)
	error = compute_error(S[0],S[1], f, f_save, dt);
      adjust_dt(error, finite, dt);
      
      done = step_success(error, finite);
    }
  }

};
