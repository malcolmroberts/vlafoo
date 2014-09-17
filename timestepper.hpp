#include <iostream>
#include <stdlib.h>
#include <assert.h>  

template<class T>
class timestepper
{
private:
  T** S; // array of source buffers
  T* f0;
  int rk_stages;
  bool rk_allocated, f0_allocated;
  bool dynamic;
  enum RKTYPE{EULER,RK2};
  RKTYPE rk;
  double tolmin, tolmax;
  double dtmax;
protected:
  unsigned int rk_n; // number of data points
  
public:
  timestepper()
  {
    rk_allocated=false;
    f0_allocated=false;
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
		   bool dynamic0, 
		   double tolmin0, double tolmax0,
		   double dtmax0)
  {
    dtmax=dtmax0;
    rk_n=n;
    rk_stages=0;
    dynamic=dynamic0;
    tolmin=tolmin0;
    tolmax=tolmax0;

    rk_allocated=false;
    f0_allocated=false;

    assert(tolmin < tolmax);

    if(rk_name == "euler") {
      rk_stages=1;
      rk=EULER;
      dynamic=false;
    }
    if(rk_name == "rk2") {
      rk_stages=2;
      rk=RK2;
      f0=new T[rk_n];
    }

    if(rk_stages == 0) {
      std::cerr << "ERROR: unknown integrator " 
		<< rk_name << std::endl;
      exit(1);
    }

    S=new T*[rk_stages];
    for(int i=0; i < rk_stages; i++)
      S[i] = new T[rk_n];
    rk_allocated=true;
    f0_allocated=true;
  }

  ~timestepper()
  {
    if(rk_allocated) {
      for(int i=0; i < rk_stages; i++)
	delete S[i];
      delete[] S;
    }
    if(f0_allocated)
      delete[] f0;

  }

  virtual void rk_source(T *f, T *S){}
  
  void rk_step(T* f, double &dt)
  {
    switch(rk) {
    case EULER:
      euler_step(f,dt);
      break;
    case RK2:
      rk2_step(f,dt);
      break;
    default:
      std::cerr << "ERROR: unknown integrator in rk_step." << std::endl;
      exit(1);
      break;
    }
  }

  void euler_step(T* f, const double dt)
  {
    rk_source(f,S[0]);
    T* S0=S[0];
    for(unsigned int i=0; i < rk_n; i++)
      f[i] += dt*S0[i];
  }

  void rk2_step(T* f, double &dt)
  {
    T* s;

    bool done=false;

    while(!done) {
      s = S[0];
      rk_source(f,s);
      double halfdt = 0.5 * dt;
      for(unsigned int i = 0; i < rk_n; i++)
	f0[i] = f[i] + halfdt * s[i];
      
      s = S[1];
      rk_source(f0,s);
      for(unsigned int i = 0; i < rk_n; i++)
	f0[i] += dt*s[i];
      
      if(dynamic) {
	T error=0.0;
	{
	  T* S0=S[0];
	  T* S1=S[1];
	  const T eps=1e-10;
	  for(unsigned int i=0; i < rk_n; i++) {
	    T S0i = S0[i];
	    T S1i = S1[i];
	    T fi = f[i];
	    //T f0i = f0[i];
	    //T diff = dt*abs(S0i-S1i)
	    //   / (max(abs(fi+dt*S0i),abs(fi+dt*S1i)) + eps);
	    T diff = abs(S0i-S1i) / (abs(fi) + eps);
	    if(diff > error)
	      error = diff;
	  }
	}
	
	bool isfinite=true;
	for(unsigned int i = 0; i < rk_n; i++) {
	  if(isnan(f0[i]) || isinf(f0[i]))
	    isfinite=false;
	}
 
	if(!isfinite)
	  std::cout << isfinite << std::endl;
	
	if(error < tolmin && isfinite) {
	  dt *= 1.2;
	  done=true;
	  for(unsigned int i = 0; i < rk_n; i++)
	    f[i] = f0[i];
	}

	if(!isfinite || error > tolmax)
	  dt *= 0.7;
	if(isfinite)      
	  done=true;

	if(dtmax > 0 && dt > dtmax)
	  dt=dtmax;

      }
    }
  }
};
