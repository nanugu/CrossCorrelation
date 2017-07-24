#include "libCorrelationAlgorithmsGenerationV4.i"
require, "yoco.i";

func fit_poly(x, a){
  return a(1)+a(2)*1/x;
}

//uses: Test_Systematic_errors_with_increasing_interpolation.i

Data=fits_read("Fits/Data_bias_vs_interpolation.fits");


yocoNmCreate, 19, 2, 1, square=1, dx=0.06,dy=0.045, fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;
fma;


Wid=1;
///Plotting
dS=0.5;

plsys,1;
index=1;
//deviation
plp, Data(,2), Data(,1),  symbol=9, size=dS, width=Wid;
plg, Data(,2), Data(,1);

plp, Data(,3), Data(,1),  symbol="^", size=dS, color="blue", width=Wid;
plg, Data(,3), Data(,1),   color="blue";

plp, Data(,4), Data(,1), color="red", symbol="*", size=dS, width=Wid;
plg, Data(,4), Data(,1),  color="red";


plp, Data(,5), Data(,1), color="green", symbol="x", size=dS, width=Wid;
plg, Data(,5), Data(,1),  color="green";

plp, Data(,6), Data(,1), color="magenta", symbol="#", size=dS, width=Wid;
plg, Data(,6), Data(,1),  color="magenta";

plp, Data(,7), Data(,1), color="black", symbol="v", size=dS, width=Wid;
plg, Data(,7), Data(,1),  color="black";

temp_x=span(1.0,10.0, 10);
a=[1, 1];
r= lmfit(fit_poly, temp_x, a, Data(,3), 1.);
plg, fit_poly(temp_x,a), temp_x, type=2, color="cyan", width=4*Wid;

  
//xytitles, "Improvemnt of sampling", "deviation (!d) [pixel]", [0.01, 0.01];
//pltitle, "Systematic error";


plp,   0.143,  4.0, symbol="^", size=dS, color="blue", width=Wid;
plp,   0.133,  4.0, symbol="*", size=dS, color="red", width=Wid;
plp,   0.123,  4.0, symbol="x", size=dS, width=Wid, color="green";
plp,   0.113,  4.0, symbol=9, size=dS, width=Wid, color="black";
plp,   0.103,  4.0, symbol="#", size=dS, width=Wid,color="magenta";
plp,   0.095,  4.0, symbol="v", size=dS, width=Wid;

plt, "--",  3.1, 0.14, tosys=1, height=12, color="blue";
plt, "--",  3.1, 0.13, tosys=1, height=12, color="red";
plt, "--",  3.1, 0.12, tosys=1, height=12, color="green";
plt, "--",  3.1, 0.11, tosys=1, height=12, color="black";
plt, "--",  3.1, 0.101, tosys=1, height=12, color="magenta";
plt, "--",  3.1, 0.091, tosys=1, height=12, color="black";




plt, "  -- LGS+corr. step2",  4.0, 0.14, tosys=1, height=14, color="blue", font="times";
plt, "  -- CF+corr. step2",  4.0, 0.13, tosys=1, height=14, color="red", font="times";
plt, "  -- Solar+corr. step2",  4.0, 0.12, tosys=1, height=14, color="green", font="times";
plt, "  -- Point+corr. step2",  4.0, 0.11, tosys=1, height=14, color="black", font="times";
plt, "  -- Point+TCoG step2",  4.0, 0.10, tosys=1, height=14, color="magenta", font="times";
plt, "  -- LGS+TCoG step2",  4.0, 0.09, tosys=1, height=14, color="black", font="times";
plt, swrite(format="---- %0.2f/K", a(2)),  3.0, 0.08, tosys=1, height=14, color="cyan", font="times";


limits, 0.0, 10.5, -0.02, 0.16;



write, "Improvement in percentage for point, LGS, GC and solar =", (Data(,2)(1)-Data(,2)(5))/Data(,2)(1)*100, (Data(,3)(1)-Data(,3)(5))/Data(,3)(1)*100, (Data(,4)(1)-Data(,4)(5))/Data(,4)(1)*100, (Data(,5)(1)-Data(,5)(5))/Data(,5)(1)*100;

write, "Improvement in factor for point, LGS, GC and solar =", (Data(,2)(1))/Data(,2)(5), (Data(,3)(1))/Data(,3)(5), (Data(,4)(1))/Data(,4)(5), (Data(,5)(1))/Data(,5)(5);


pdf, "Fig7_Interpolation_sub-pixel_resolutionV4_1";






