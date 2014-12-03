import graph;
import palette;
import cooltowarm;
include "./vlafoo.asy";

import animate; 
//settings.render=2; 
//settings.tex="pdflatex"; 
//settings.prc=false;
//settings.thick=false;
settings.outformat="mp4";
 
size(10cm,10cm,IgnoreAspect); 

animation A;

string run=getstring("run");

eval("include \""+run+"/vlafoo.asy\";",true);

string basename=getstring("basename");
//string basename=getstring("basename");
int n=getint("n");
//real vmax=6.0;

real vcut=0;

usersetting();

if(vcut == 0)
  vcut=getreal("vcut");

pair a=(0,-vmax);
pair b=(Lx,vmax);

pen[] Palette=BWRainbow();

picture pic;

pair p0=(0,0);
pair p1=(1,1);

real[][][] ff;

// find the image pallette bounds globallly:
real fmin=infinity, fmax=-infinity;
real t[] = new real[n];

for(int i=0; i < n; ++i) { 

  string filename=run+"/"+basename+string(i);
  write(filename);
  //file fin=input(filename).line();
  //real[][] f=fin.dimension(0,0);
  
  file fin=input(filename,mode="xdr");

  t[i]=fin;
  
  real[][] f=fin.read(2);
 
  ff[i]=f;
  
  if(min(f) < fmin) fmin=min(f);
  if(max(f) > fmax) fmax=max(f);
} 

write("image bounds:",fmin,fmax);
Palette=cooltowarm(fmin,0,fmax);

// re-read the files and make the movie:
for(int i=0; i < n; ++i) { 
  add(pic); 

  //string filename=basename+string(i);
  //file fin=input(filename).line();
  real[][] f=ff[i];//fin.dimension(0,0);

  bounds range;
  
  range=image(f,Range(fmin,fmax),a,b,Palette); //global full colour bar
  //range=image(f,Full,a,b,Palette);  // Full colour bar

  if(vcut > 0 && vcut < vmax) {
    real f=vcut/vmax;
    limits((a.x,f*a.y),(b.x,f*b.y),Crop);
  }
  
  
  xaxis("$x$",BottomTop,LeftTicks,above=true);
  yaxis("$v$",LeftRight,RightTicks,above=true);

  // Add the palette bar:
  picture bar;
  string barlegend="";
  real paletteheight=6cm;
  palette(bar,barlegend,range,(0,0),(0.5cm,paletteheight),Right,Palette,
	  PaletteTicks(ptick=linewidth(0.5*linewidth())));
  add(bar.fit(),point(E),30E);

  draw(Label("\texttt{"+texify(run+"/"+basename)+"}, $t="+string(t[i])+"$"),point(S),10S);
       //  draw(Label("\texttt{"+texify(run+"/"+basename)+"} t="+string(t[i]))),point(S),10S);
  
  A.add(); 
  erase(); 
} 

// 100/delay FPS
//A.movie(delay=50,options="-density 288x288 -geometry 50%x");
//A.movie(delay=50,options="-density 288x288");
A.movie(delay=50,options="-density 144x144");
//A.movie(options="-density 288x288 -geometry 50%x");
