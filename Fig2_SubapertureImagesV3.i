require, "yoco.i";
require, "libImages.i";


Point_Sub=MakePoint(9.0, 9.0);
LGS_Sub=MakeLGS(9.0, 9.0);
GC_Sub=MakeGC(8.0, 8.0);
Solar_Sub=MakeSolar(8, 8);

yocoNmCreate, 20, 4, 1, square=1, dx=0.05, dy=0.05, fx=1, fy=1, ytitleSpace=0.1;
fma;
get_style,landscape,systems,legends,clegends;

plsys,1;
pli, Point_Sub;
//xytitles, "X (pixels)", "Y (pixels)", [0.02, 0.02];
//pltitle, "Point source";
//plt, "a",  2, 15, tosys=1, height=16, color="yellow", font="timesB", justify="LT";
limits, 1, 16, 1, 16;



plsys,2;
pli, LGS_Sub;
//xytitles, "X (pixels)", "",[0.0, 0.02];
//pltitle, "LGS elongated";
//plt, "b",  2, 15, tosys=1, height=16, color="yellow", font="timesB", justify="LT";
limits, 1, 16, 1, 16;


plsys,3;
pli,  GC_Sub;
//xytitles, "X (pixels)", "",[0.0, 0.02];
//pltitle, "Galactic Center";
//plt, "c",  2, 15, tosys=1, height=16, color="yellow", font="timesB", justify="LT";
limits, 1, 16, 1, 16;

plsys,4;
pli,  Solar_Sub;
//xytitles, "X (pixels)", "",[0.0, 0.02];
//pltitle, "Solar photosphere";
//plt, "d",  2, 15, tosys=1, height=16, color="yellow", font="timesB", justify="LT";
limits, 1, 16, 1, 16;

pdf, "Fig1_Subaperture_imagesV2";

