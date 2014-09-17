import graph;
import utils;

size(300,200,IgnoreAspect);

// specify legend for run names:
// asy stats -u "runlegs=\"asdf,asdf1\""

string runs=getstring("run");
string file="stats";
string val=getstring("val");

string xlabel="$t$";
string ylabel=val;

scale(Linear,Linear);

if(val == "dt") scale(Linear,Log);
if(file == "tnow") scale(Linear,Linear);
  
bool docrop=false;

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
    write(filename);
    filename += "/";
    filename += file;
    lastpos=pos > 0 ? pos+1 : -1;
    
    file fin=input(filename).line();
    real[][] a=fin.dimension(0,0);
    a=transpose(a);
    real[] it=a[0];
    real[] twow=a[1];
    real[] dt=a[1];

    string legend= myleg ? legends[n] : texify(filename);
    draw(graph(it,a[1]),Pen(n),legend);
  }
}

xaxis(xlabel,BottomTop,LeftTicks);
yaxis(ylabel,LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);
