//uses Test_Fig2_SubPixel_Fitting_all_2017Jan29.i

require, "yoco.i";
Data1=fits_read("Fits/Data_performance_subpixel_pointV2_5_5.fits");
Data2=fits_read("Fits/Data_performance_subpixel_LGSV2_5_5.fits");
Data3=fits_read("Fits/Data_performance_subpixel_GCV2_5_5.fits");
Data4=fits_read("Fits/Data_performance_subpixel_solarV2_5_5.fits");

ytop= 0.33;
ybottom= -ytop;
xleft= -1;
xright= -xleft;

e=dimsof(Data1)(2);
dS= 0.5;
Wid= 1.5;

//define data colours
fitcode= [2, 3, 5, 6, 7, 4];
fitname= ["Corr + QPF", "Corr + PF", "Corr + PYF", "Corr + GF", "Corr + TCoG", "TCoG"];
colours= ["blue", "red", "black", "green", "magenta", "cyan"];
symbols= [9, 6, 3, 8, 1, 2];

winkill,1;
yocoNmCreate, 1, 2, 2, square=1, dx=0.04, fy=[1,0], fx= [0,1]; //style="boxed_times.gs";
yocoNmLimits, xleft, xright, ybottom, ytop;
fma;

for (i=1; i<=4; i++) {

  if (i == 1) {data= Data1; yocoNmPlsys, 1, 1;} //top left: point source
  else if (i == 2) {data= Data2; yocoNmPlsys, 2, 1;} //top right: laser guide star
  else if (i == 3) {data= Data3; yocoNmPlsys, 1, 2;} //bottom left: galactic centre
  else if (i == 4) {data= Data4; yocoNmPlsys, 2, 2;} //bottom right: solar
  
  plg, [0,0], [xleft, xright], type="dot";
  plg, [ybottom, ytop], [0.5, 0.5], type= "dot";
  plg, [ybottom, ytop], [-0.5, -0.5], type= "dot";
  
  
  for (f=1; f<=5; f++)
    plp, data(,fitcode(f)), data(,1), symbol=symbols(f), color= colours(f),
      size= dS, width= Wid;

  if(i<3){
    f=6;
    plp, data(,fitcode(f)), data(,1), symbol=symbols(f), color= colours(f), size= dS, width= Wid;
  }
 }


yocoNmPlsys, 1, 1;

dy=-0.035;
for (f=1; f<=6; f++) {
  x= -0.5;
  y= 0.33+f*dy;
  plp, y, x, symbol= symbols(f), color= colours(f), size= dS, width= Wid;
  plt, fitname(f), x+0.1, y-0.01, height= 12, color= colours(f), tosys=1;
 }



pdf, "Fig2_Subpixel_performance_20160921";
