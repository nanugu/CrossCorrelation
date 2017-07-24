//uses: Test_Fig3_5_SNR_vs_bias_vs_sigma_for_all_objects.i


#include "yoco.i"
Data=fits_read("Fits/point_SNR.fits");
Data=abs(Data(,:));

dS=0.5;
Wid=1.5;

yocoNmCreate, 1, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys, 1;

P=where(Data(,1)<=40);

plsys, 1;
plp, Data(,8)(P), Data(,1)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,9)(P), Data(,1)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,10)(P), Data(,1)(P),  color="cyan", symbol="+", size=dS, width=Wid;
plp, Data(,11)(P), Data(,1)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,12)(P), Data(,1)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,13)(P), Data(,1)(P),  symbol="square", size=dS, color="magenta", width=Wid;

/*
plp,  0.28+0.05, 20, symbol="square", size=dS, color="magenta", width=Wid;
plp,  0.24+0.05, 20, symbol=9, size=dS, color="blue", width=Wid;
plp,  0.20+0.05, 20, symbol="x", size=dS, width=Wid, color="red";
plp,  0.16+0.05,  20, symbol="^", size=dS, color="black", width=Wid;
plp,  0.12+0.05, 20, symbol="*", size=dS, color="green", width=Wid;


plt, "  -- TCoG",  20, 0.272+0.05, tosys=1, height=12, color="magenta", font="times";
plt, "  -- QPF",  20, 0.232+0.05, tosys=1, height=12, color="blue", font="times";
plt, "  -- PF",  20, 0.192+0.05, tosys=1, height=12, color="red", font="times";
plt, "  -- PYF",  20, 0.152+0.05, tosys=1, height=12, color="black", font="times";
plt, "  -- GF",  20, 0.112+0.05, tosys=1, height=12, color="green", font="times";
*/
limits, 0, 42, -0.0, 0.5;
pdf, "Sigma_SNR_Point";



yocoNmCreate, 2, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys, 1;
plp, Data(,2)(P), Data(,1)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,3)(P), Data(,1)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,4)(P), Data(,1)(P),  color="cyan", symbol="+", size=dS, width=Wid;
plp, Data(,5)(P), Data(,1)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,6)(P), Data(,1)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,7)(P), Data(,1)(P),  symbol="square", size=dS, color="magenta", width=Wid;
limits, 0, 42, -0.02, 0.3; 
pdf, "Bias_SNR_Point";



X=span(0, 0.5, 100); 
yocoNmCreate, 3, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys, 1;
plp, Data(,2)(P), Data(,8)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,3)(P), Data(,9)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,4)(P), Data(,10)(P),  color="cyan", symbol="+", size=dS, width=Wid;
plp, Data(,5)(P), Data(,11)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,6)(P), Data(,12)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,7)(P), Data(,13)(P),  symbol="square", size=dS, color="magenta", width=Wid;

//plg, [0.05,0.05], [0, 1], type="dot";
//plg, [-0.3, 0.3], [0.05, 0.05], type= "dot";

plg, X, X, type="dot";


limits, -0.01, 0.5, -0.02, 0.3; 
pdf, "Bias_Sigma_Point";
