#include <iostream>
#include <stdlib.h>

template<class T>
class timestepper
{
private:
  T** S; // array of source buffers
  int rk_stages;
  bool rk_allocated;
  bool dynamic;
  enum RKTYPE{EULER,RK2};
  RKTYPE rk;
  double tolmin, tolmax;
protected:
  unsigned int rk_n; // number of data points
  
public:
  timestepper()
  {
    rk_allocated=false;

    // FIXME: set parameters via command-line
    dynamic=true;
    tolmin=0.03;
    tolmax=0.05;
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

  void rk_allocate(int n, std::string &rk_name)
  {
    rk_n=n;
    rk_stages=0;
    if(rk_name == "euler") {
      rk_stages=1;
      rk=EULER;
      dynamic=false;
    }
    if(rk_name == "rk2") {
      rk_stages=2;
      rk=RK2;
    }

    if(rk_stages == 0) {
      std::cerr << "ERROR: unknown integrator " 
		<< rk_name << std::endl;
      exit(1);
    }

    S=new T*[rk_stages];
    for(int i=0; i < rk_stages; i++)
      S[i] = new T[rk_n] ;
    rk_allocated=true;
  }

  ~timestepper()
  {
    if(rk_allocated) {
      for(int i=0; i < rk_stages; i++)
	delete S[i];
      delete[] S;
    }
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

  void euler_step(T* f, double &dt)
  {
    rk_source(f,S[0]);
    T* S0=S[0];
    for(unsigned int i=0; i < rk_n; i++)
      f[i] += dt*S0[i];
  }

  void rk2_step(T* f, double &dt)
  {
    T* s=S[0];
    rk_source(f,s);
    double halfdt=0.5*dt;
    for(unsigned int i=0; i < rk_n; i++)
      f[i] += halfdt*s[i];
    
    s=S[1];
    rk_source(f,s);
    for(unsigned int i=0; i < rk_n; i++)
      f[i] += dt*s[i];

    T error=0.0;
    {
      T* S0=S[0];
      T* S1=S[1];
      T eps=1e-10;
      for(unsigned int i=0; i < rk_n; i++) {
	T S0i = S0[i];
	T S1i = S1[i];
	T diff = abs(S0i-S1i) / (max(abs(S0i),abs(S1i)) + eps);
	if(diff > error)
	  error = diff;
      }

    }
    std::cout << error << std::endl;

    // TODO: adaptive version
  }
  
};
