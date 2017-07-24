require, "yoco.i";
Data1=fits_read("Fits/Data_Two_step_traditional_point.fits");
Data2=fits_read("Fits/Data_Two_step_traditional_LGS.fits");
Data3=fits_read("Fits/Data_Two_step_traditional_GC.fits");
Data4=fits_read("Fits/Data_Two_step_traditional_solar.fits");

ytop= 0.22;
ybottom= -ytop;
xleft= -1;
xright= -xleft;

e=dimsof(Data1)(2);
dS= 0.5;
Wid= 1.5;

//define data colours
fitcode= [2, 3,  5, 6];
fitname= ["Corr. conventional", "Corr. two-step", "TCoG conventional","TCoG two-step"];
colours= ["blue", "red", "black", "green"];
symbols= [4, 5,  6, 8];

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
  
  if(i<3){
    for (f=1; f<=4; f++){
      plp, data(,fitcode(f)), data(,1), symbol=symbols(f), color= colours(f), size= dS, width= Wid;
    }
  }
  
  if(i>2){
    for (f=1; f<=2; f++){
      plp, data(,fitcode(f)), data(,1), symbol=symbols(f), color= colours(f), size= dS, width= Wid;
    }
  }
  
  
 }


yocoNmPlsys, 1, 1;

dy=-0.03;
for (f=1; f<=2; f++) {
  x= -0.5;
  y= 0.22+f*dy;
  plp, y, x, symbol= symbols(f), color= colours(f), size= dS, width= Wid;
  plt, fitname(f), x+0.1, y-0.01, height= 13, color= colours(f), tosys=1;
 }

for (f=3; f<=4; f++) {
  x= -0.5;
  y= -0.08+f*dy;
  plp, y, x, symbol= symbols(f), color= colours(f), size= dS, width= Wid;
  plt, fitname(f), x+0.1, y-0.01, height= 13, color= colours(f), tosys=1;
 }

//pdf, "Fig5_Comparision_two_step_traditional_1.pdf";


yocoNmPlsys, 2, 2;
/*
yocoNmCreate, 2, 2, 2, square=1, dx=0.04, fy=[1,0], fx= [0,1]; //style="boxed_times.gs";
yocoNmLimits, xleft, xright, ybottom, ytop;
fma;
*/

for (i=3; i<=4; i++) {

  if (i == 1) {data= Data1; yocoNmPlsys, 1, 1;} //top left: point source
  else if (i == 2) {data= Data2; yocoNmPlsys, 2, 1;} //top right: laser guide star
  else if (i == 3) {data= Data3; yocoNmPlsys, 1, 2;} //bottom left: galactic centre
  else if (i == 4) {data= Data4; yocoNmPlsys, 2, 2;} //bottom right: solar
  
  plg, [0,0], [xleft, xright], type="dot";
  plg, [ybottom, ytop], [0.5, 0.5], type= "dot";
  plg, [ybottom, ytop], [-0.5, -0.5], type= "dot";
  

  for (f=1; f<=2; f++){
     
      plp, data(,fitcode(f)), data(,1), symbol=symbols(f), color= colours(f),size= dS, width= Wid;
  }
}


pdf, "Fig6_Comparision_two_step_traditional_2.pdf";
