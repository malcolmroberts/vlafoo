/***********************************************************************

  Color maps based on the diverging color maps by Kenneth Moreland
  www.sandia.gov/~kmorel/documents/ColorMaps/

  Copyright 2014 Malcolm Roberts, www.malcolmiwroberts.com

  -------------------------------------------------------------------

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  
************************************************************************/

import palette;

pen[] cooltowhite_pens={rgb(59,76,192),
		       rgb(68,90,204),
		       rgb(77,104,215),
		       rgb(87,117,225),
		       rgb(98,130,234),
		       rgb(108,142,241),
		       rgb(119,154,247),
		       rgb(130,165,251),
		       rgb(141,176,254),
		       rgb(152,185,255),
		       rgb(163,194,255),
		       rgb(174,201,253),
		       rgb(184,208,249),
		       rgb(194,213,244),
		       rgb(204,217,238),
		       rgb(213,219,230)};

pen[] cooltowarm_white={rgb(221,221,221)};

pen[] whitetowarm_pens={rgb(229,216,209),
		       rgb(236,211,197),
		       rgb(241,204,185),
		       rgb(245,196,173),
		       rgb(247,187,160),
		       rgb(247,177,148),
		       rgb(247,166,135),
		       rgb(244,154,123),
		       rgb(241,141,111),
		       rgb(236,127,99),
		       rgb(229,112,88),
		       rgb(222,96,77),
		       rgb(213,80,66),
		       rgb(203,62,56),
		       rgb(192,40,47),
		       rgb(180,4,38)};
			
pen[] reverse_pens(pen[] P) {
  pen[] parray;
  for(int i=0; i < P.length; ++i) {
    parray[i]=P[P.length-i-1];
  }
  return parray;
}

// A palette varying linearly over the specified array of pens, using
// NColors in each interpolation interval.
pen[] Gradient(int NColors=256, pen[] p) 
{
  pen[] P;
  if(p.length < 2) abort("at least 2 colors must be specified");
  real step=NColors > 1 ? (1/(NColors-1)) : 1;
  for(int i=0; i < p.length-1; ++i) {
    pen begin=p[i];
    pen end=p[i+1];
    P.append(sequence(new pen(int j) {
          return interp(begin,end,j*step);
        },NColors));
  }
  return P;
}

pen[] cooltowarm=Gradient(concat(cooltowhite_pens,
				 cooltowarm_white,
				 whitetowarm_pens));

pen[] warmtocool=Gradient(reverse_pens(concat(cooltowhite_pens,
					     cooltowarm_white,
					     whitetowarm_pens)));

pen[] cooltowhite=Gradient(256,concat(cooltowhite_pens,
				      cooltowarm_white));

pen[] whitetocool=Gradient(256,reverse_pens(concat(cooltowhite_pens,
						  cooltowarm_white)));

pen[] whitetowarm=Gradient(concat(cooltowarm_white,
				  whitetowarm_pens));

pen[] warmtowhite=Gradient(reverse_pens(concat(cooltowarm_white,
					      whitetowarm_pens)));

// Passing min, mean, and max values sets the white point to represent
// the mean value with the colour map truncated if the range is
// greater in one direction than the other.
pen[] cooltowarm(real min, real mean, real max)
{
  real ratio=abs(max-mean)/abs(mean-min);
  if(ratio > 1) {
    // the data range is greater between max and mean
    int len=cooltowhite_pens.length;
    int start=len-floor(len/ratio);
    pen[] cooltowhite_subarray;
    for(int i=start; i < len; ++i)
      cooltowhite_subarray[i-start]=cooltowhite_pens[i];
    return Gradient(256,concat(cooltowhite_subarray,
			       cooltowarm_white,
			       whitetowarm_pens));
  }
  if(ratio < 1) {
    // the data range is greater between min and mean
    int len=whitetowarm_pens.length;
    int stop=floor(ratio*len);
    pen[] whitetowarm_subarray;
    for(int i=0; i < stop; ++i)
      whitetowarm_subarray[i]=whitetowarm_pens[i];
    return Gradient(256,concat(cooltowhite_pens,
			       cooltowarm_white,
			       whitetowarm_subarray));
  }
  return cooltowarm;
}

pen[] warmtocool(real min, real mean, real max)
{
  return(reverse_pens(cooltowarm(max,mean,min)));
}
