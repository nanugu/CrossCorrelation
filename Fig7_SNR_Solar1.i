//uses: Test_Fig3_5_SNR_vs_bias_vs_sigma_for_all_objects.i

#include "yoco.i"
Data=fits_read("Fits/solar_SNR.fits");

Data=abs(Data(,:));

dS=0.5;
Wid=1.5;

yocoNmCreate, 1, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

Data(,1)=Data(,1); //SNR per image

P=where(Data(,1)<=800);

plsys, 1;
plp, Data(,8)(P), Data(,1)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,9)(P), Data(,1)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,11)(P), Data(,1)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,12)(P), Data(,1)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,13)(P), Data(,1)(P),  symbol="square", size=dS, color="magenta", width=Wid;
//xytitles, "SNR", "RMS error [pixels]";

limits, 0, 800, -0.01, 0.5;

pdf, "Sigma_SNR_Solar";



yocoNmCreate, 2, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys, 1;
plp, Data(,2)(P), Data(,1)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,3)(P), Data(,1)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,5)(P), Data(,1)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,6)(P), Data(,1)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,7)(P), Data(,1)(P),  symbol="square", size=dS, color="magenta", width=Wid;
limits, -0.01,800, -0.03, 0.3; 

pdf, "Bias_SNR_Solar";


yocoNmCreate, 3, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

X=span(0, 0.5, 100); 
plsys, 1;
plp, Data(,2)(P), Data(,8)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,3)(P), Data(,9)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,5)(P), Data(,11)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,6)(P), Data(,12)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,7)(P), Data(,13)(P),  symbol="square", size=dS, color="magenta", width=Wid;
limits, -0.01, 0.5, -0.03, 0.3;
plg, X, X, type="dot";
pdf, "Bias_Sigma_Solar";

