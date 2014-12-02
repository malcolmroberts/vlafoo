import graph;
import utils;

size(300,200,IgnoreAspect);

// specify legend for run names:
// asy stats -u "runlegs=\"asdf,asdf1\""

string runs=getstring("run");
string file="stats";
string xval=getstring("xval");
string yval=getstring("yval");

string xlabel=xval;
string ylabel=yval;

int xcol=0;
int ycol=1;
string sscale="";

scale(Linear,Linear);

if(xval == "it") xcol=0;
if(xval == "t")  xcol=1;
if(xval == "dt") xcol=2;

if(yval == "it") ycol=0;
if(yval == "t")  ycol=1;
if(yval == "dt") ycol=2;

if(sscale == "")
  sscale=getstring("scale linlin or loglin or linlog or loglog");
if(sscale == "loglog") scale(Log,Log);
if(sscale == "loglin") scale(Log,Linear);
if(sscale == "linlog") scale(Linear,Log);
if(sscale == "linlin") scale(Linear,Linear);


bool docrop=false;

string runlegs="";
bool doleg=true;

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
    real[] dt=a[2];

    string legend= myleg ? legends[n] : texify(filename);
    draw(graph(a[xcol],a[ycol],a[ycol]>1e-10),Pen(n),legend);
  }
}

xaxis(xlabel,BottomTop,LeftTicks);
yaxis(ylabel,LeftRight,RightTicks);

if(doleg)
  attach(legend(),point(plain.E),20plain.E);
