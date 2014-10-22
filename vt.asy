import graph;
import utils;
include "./vlafoo.asy";

size(300,200,IgnoreAspect);

// specify legend for run names:
// asy vt -u "runlegs=\"asdf,asdf1\""

string runs=getstring("run");
string file=getstring("file");

string xlabel="$t$";
string ylabel="$y$";
scale(Linear,Log);

if(file == "ekvt")
  ylabel="Kinetic energy";
if(file == "eEvt")
  ylabel="Electrical energy";
if(file == "fminvt") {
  ylabel="minimum density";
  scale(Linear,Linear);
}
if(file == "fmaxvt") {
  ylabel="maximum density";
  scale(Linear,Linear);
}
if(file == "ftotvt") {
  scale(Linear,Linear);
  ylabel="total number of particles";
}
  
bool docrop=false;
real tmax=infinity;


string runlegs="";
usersetting();

bool myleg=((runlegs== "") ? false: true);
string[] legends=set_legends(runlegs);



string filename;
int n=-1;
bool flag=true;
int lastpos;
while(flag) {
  ++n;
  int pos=find(runs,",",lastpos);
  if(lastpos == -1) {filename=""; flag=false;}
  filename=substr(runs,lastpos,pos-lastpos);
  if(flag) {

    // Include the asy information
    eval("include \""+filename+"/vlafoo.asy\";",true);

    filename += "/";
    filename += file;
    write(filename);
    
    lastpos=pos > 0 ? pos+1 : -1;
    
    file fin=input(filename).line();
    real[][] a=fin.dimension(0,0);
    a=transpose(a);
    real[] t=a[0];
    real[] Ek=a[1];

    write("max time="+string(max(t)));

    string legend= myleg ? legends[n] : texify(filename);
    draw(graph(t,Ek,t<=tmax),Pen(n),legend);
    if(docrop) ylimits(7e-2,1e-1,Crop);
  }
}


xaxis(xlabel,BottomTop,LeftTicks);
if(docrop)
  // TODO: see fftw++'s mpi/timing.asy for a better way to deal with
  // the case where we don't have enough decades of data.
  yaxis(ylabel,LeftRight,
	LeftTicks(DefaultFormat,new real[] {1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,7e-2,8e-2,9e-2,1e-1}));
else
yaxis(ylabel,LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);
