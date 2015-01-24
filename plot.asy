size(10cm,10cm,IgnoreAspect);

import graph;
import palette;
import cooltowarm;
include "./vlafoo.asy";

string run=getstring("run");

eval("include \""+run+"/vlafoo.asy\";",true);


string filename=getstring("filename");
//file fin=input(filename).line();
file fin=input(run+"/"+filename,mode="xdr");

real vcut=0;

usersetting();

if(vcut == 0)
  vcut=getreal("vcut");

real t=fin;

real[][] f;
f=fin.read(2);

write("time="+string(t));

int nx=10;
int nv=10;


pair a=(0,-vmax);
pair b=(Lx,vmax);

write("min=",min(f));
write("max=",max(f));

pen[] Palette;
//Palette=cooltowhite;
//Palette=whitetocool;
//Palette=whitetowarm;
//Palette=warmtowhite;
//Palette=cooltowarm;
Palette=cooltowarm(min(f),0,max(f));
//Palette=warmtocool(min(f),0,max(f));
//Palette=warmtocool;
//Palette=BWRainbow();


int nneg=0; // number of negative points
real negmean=0; // mean value over negative points
real negcut=0;//-0.001;
for(int i=0; i < f.length; ++i) {
  for(int j=0; j < f[i].length; ++j) {
    if(f[i][j] < negcut) {
      ++nneg;
      negmean += f[i][j];
      //f[i][j]=1;
    } else {
      //f[i][j]=0;
    }
  }
}

bounds range;
range=image(f,Full,a,b,Palette);  // Full colour bar

if(vcut > 0 && vcut < vmax) {
  real f=vcut/vmax;
  limits((a.x,f*a.y),(b.x,f*b.y),Crop);
}
  
if(nneg > 0) {
  //write(nneg/(f.length*f[0].length));
  negmean /= nneg;
  write("mean of negative points: "+string(negmean));
}

xaxis("$x$",BottomTop,LeftTicks,above=true);
yaxis("$v$",LeftRight,RightTicks,above=true);

// Add the palette bar:
picture bar;
string barlegend = "";
real paletteheight = 6cm;
palette(bar, barlegend, range, (0,0), (0.5cm, paletteheight), Right, Palette,
        PaletteTicks(ptick = linewidth(0.5 * linewidth())));
add(bar.fit(),point(E),30E);

draw(Label("\texttt{"+texify(run+"/"+filename)+"}, $t="+string(t)+"$"),
     point(S),10S);
