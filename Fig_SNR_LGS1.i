//uses: Test_Fig3_5_SNR_vs_bias_vs_sigma_for_all_objects.i


#include "yoco.i"
Data=fits_read("Fits/Two_step_vs_SNR_LGS.fits");
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
plp, Data(,11)(P), Data(,1)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,12)(P), Data(,1)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,13)(P), Data(,1)(P),  symbol="square", size=dS, color="magenta", width=Wid;

limits, 0, 42, -0.01, 0.6;
pdf, "Sigma_SNR_LGS_TwoStep";



yocoNmCreate, 2, 2, 1, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys, 1;
plp, Data(,2)(P), Data(,1)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,3)(P), Data(,1)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,5)(P), Data(,1)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,6)(P), Data(,1)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,7)(P), Data(,1)(P),  symbol="square", size=dS, color="magenta", width=Wid;
limits, 0, 42, -0.02, 0.1; 
pdf, "Bias_SNR_LGS_TwoStep";




yocoNmCreate, 3, 2, 3, square=1, dx=0.15,  dy=0.15,  fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys, 1;
plp, Data(,2)(P), Data(,8)(P), color="blue", symbol=9, size=dS, width=Wid;
plp, Data(,3)(P), Data(,9)(P),  color="red", symbol="x", size=dS, width=Wid;
plp, Data(,5)(P), Data(,11)(P),  symbol="^", size=dS, width=Wid;
plp, Data(,6)(P), Data(,12)(P),  symbol="*", size=dS, color="green", width=Wid;
plp, Data(,7)(P), Data(,13)(P),  symbol="square", size=dS, color="magenta", width=Wid;
limits, -0.01, 0.5, -0.02, 0.2; 
pdf, "Bias_sigma_LGS_TwoStep";
