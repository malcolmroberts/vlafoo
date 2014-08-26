size(15cm,10cm,IgnoreAspect);

import graph;
import palette;
import cooltowarm;

//scale(Linear,Log);

string filename;
filename=getstring("external data:");
file fin=input(filename).line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
real[] x=a[0];
real[] v=a[1];
real[] f0=a[2];


real vmax=6.0;
real L=31.4159;
L=getreal("L");
vmax=getreal("vmax");
pair a=(0,-vmax);
pair b=(L,vmax);

int nx=256;
int nv=201;
write(f0.length);
write(nx*nv);

real[][] f=new real[nx][nv];

for(int x=0; x < nx; ++x) {
  for(int v=0; v < nv; ++v) {
    f[x][v]=f0[(x * nv) +v];
  }
}


pen[] Palette;
Palette=cooltowarm(min(f),0,max(f));

bounds range;
range=image(f,Full,a,b,Palette);  // Full colour bar


// Add the palette bar:
picture bar;
string barlegend="";
real paletteheight=6cm;
palette(bar,barlegend,range,(0,0),(0.5cm,paletteheight),Right,Palette,
        PaletteTicks(ptick=linewidth(0.5*linewidth())));
add(bar.fit(),point(E),30E);

