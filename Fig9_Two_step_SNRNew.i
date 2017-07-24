//uses Test_Fig2_SubPixel_Fitting_all_2017Jan29.i

require, "yoco.i";
Data1=fits_read("Fits/Data_performance_subpixel_pointV2_5_5.fits");
Data2=fits_read("Fits/Data_performance_subpixel_LGSV2_5_5.fits");
Data3=fits_read("Fits/Data_performance_subpixel_GCV2_5_5.fits");
Data4=fits_read("Fits/Data_performance_subpixel_solarV2_5_5.fits");

Data1=fits_read("Fits/Two_step_vs_SNR_point.fits");
Data2=fits_read("Fits/Two_step_vs_SNR_LGS.fits");
Data3=fits_read("Fits/Two_step_vs_SNR_GC.fits");
Data4=fits_read("Fits/Two_step_vs_SNR_solar.fits");

ytop= 0.05;
ybottom= -0.01;
xleft= 0;
xright= 45;

//e=dimsof(Data1)(2);
dS= 0.5;
Wid= 1.5;

//define data colours
fitcode= [2, 3, 5, 6, 7, 4];
fitname= ["QPF", "PF", "PYF", "GF", "TCoG", "DTCoG"];
colours= ["blue", "red", "black", "green", "magenta", "cyan"];
symbols= [9, 6, 3, 8, 1, 2];

winkill,1;
yocoNmCreate, 1, 2, 2, square=1, dx=0.04, fy=[1,0], fx= [0,1]; //style="boxed_times.gs";
yocoNmLimits, xleft, xright, ybottom, ytop;
fma;

for (i=1; i<=2; i++) {

  if (i == 1) {data= Data1; yocoNmPlsys, 1, 2;} //top left: point source
  else if (i == 2) {data= Data2; yocoNmPlsys, 2, 2;} //top right: laser guide star
  else if (i == 3) {data= Data3; yocoNmPlsys, 1, 2;} //bottom left: galactic centre
  else if (i == 4) {data= Data4; yocoNmPlsys, 2, 2;} //bottom right: solar
  
  //plg, [-0.01,-0.01], [xleft, xright], type="dot";
  //plg, [0.06,0.06], [xleft, xright], type="dot";
  //plg, [-0.01, 0.06], [0.5, 0.5], type= "dot";
   
  Data=abs(data(,:));
  if (i<4){
      P=where(Data(,1)<=40);
      limits, 0, 42, -0.01, 0.05;
  }else{
      P=where(Data(,1)<=800);
      limits, 0, 805, -0.01, 0.05;
  }

  
  for (f=1; f<=6; f++)
      plp, Data(,fitcode(f))(P), Data(,1)(P), symbol=symbols(f), color= colours(f),
      size= dS, width= Wid;
 }

pdf, "Bias_error_SNR_Two_step_method_1";

//------------------------------------------------------------------




yocoNmCreate, 2, 2, 2, square=1, dx=0.04, fy=[1,0], fx= [0,1]; //style="boxed_times.gs";
//yocoNmLimits, xleft, xright, ybottom, ytop;
fma;

for (i=3; i<=3; i++) {

  if (i == 1) {data= Data1; yocoNmPlsys, 1, 1;} //top left: point source
  else if (i == 2) {data= Data2; yocoNmPlsys, 2, 1;} //top right: laser guide star
  else if (i == 3) {data= Data3; yocoNmPlsys, 1, 2;} //bottom left: galactic centre
  else if (i == 4) {data= Data4; yocoNmPlsys, 2, 2;} //bottom right: solar
  
  //plg, [-0.01,-0.01], [xleft, xright], type="dot";
  //plg, [0.06,0.06], [xleft, xright], type="dot";
  //plg, [-0.01, 0.06], [0.5, 0.5], type= "dot";
   
  Data=abs(data(,:));
  if (i<4){
      P=where(Data(,1)<=40);
      limits, 0, 42, -0.01, 0.05;
  }else{
      P=where(Data(,1)<=800);
      limits, 0, 805, -0.01, 0.05;
  }
  
  for (f=1; f<=5; f++)
      plp, Data(,fitcode(f))(P), Data(,1)(P), symbol=symbols(f), color= colours(f),
      size= dS, width= Wid;
 }


yocoNmPlsys, 2, 2;
data= Data4;
Data=abs(data(,:));
P=where(Data(,1)<=800);

for (f=1; f<=5; f++)
      plp, Data(,fitcode(f))(P), Data(,1)(P), symbol=symbols(f), color= colours(f),
      size= dS, width= Wid;

limits, 0, 850, -0.01, 0.05;

pdf, "Bias_error_SNR_Two_step_method_2";
