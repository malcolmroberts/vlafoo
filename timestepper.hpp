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

    // FIXME: set via command-line
    dynamic=true;
    tolmin=0.03;
    tolmax=0.05;
  }

  void rk_allocate(int n, std::string &rk_name)
  {
    rk_n=n;
    rk_stages=0;
    if(rk_name == "euler") {
      rk_stages=1;
      rk=EULER;
    }
    if(rk_name == "rk2") {
      rk_stages=2;
      rk=RK2;
    }

    if(rk_stages == 0) {
      std::cerr << "ERROR: unkown integrator " 
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
      std::cerr << "ERROR: unkown integrator in rk_step." << std::endl;
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
    

    double error=0.0;
    // TODO: compute error

    // TODO: adaptive version
  }
  
};
