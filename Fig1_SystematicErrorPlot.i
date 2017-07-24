require, "yoco.i";
Data1=fits_read("Fits/Point_systematicError.fits");

ytop= 0.15;
ybottom= -ytop;


e=dimsof(Data1)(2);
dS= 0.5;
Wid= 1.5;

P=where(Data1(,1) >= -1 & Data1(,1)  <= 1.0);

/*
Q=array(1, 6);
Q(1)=where(Data1(,1) == 0);
Q(2)=where(Data1(,1) == 0.22);
Q(3)=where(Data1(,1) == 0.4);
Q(4)=where(Data1(,1) == 0.6);
Q(5)=where(Data1(,1) == 0.8);
Q(6)=where(Data1(,1) == 1);
*/

//define data colours
fitcode= [2, 3, 5, 6, 7];
fitname= ["QPF", "PF", "PYF", "GF", "TCoG"];
colours= ["blue", "red", "black", "green", "magenta"];
symbols= [9, 6, 3, 8, 1];

winkill,1;
yocoNmCreate, 1, 2, 2, square=1, dx=0.04, fy=[1,0], fx= 1; //style="boxed_times.gs";
yocoNmLimits, -1.05, 1.05, ybottom, ytop;
fma;

data= Data1;
plsys, 1; //top left: point source
for (i=1; i<=1; i++) {    
    
    plg, [ybottom, ytop], [0.5, 0.5], type= "dot";
    plg, [ybottom, ytop], [-0.5, -0.5], type= "dot";
    plg, [0,0], [-1, 1], type="dot";
    
    for (f=5; f<=5; f++){
        plg, data(,fitcode(f))(P), data(,1)(P);
        //plp, data(,fitcode(f))(Q), data(,1)(Q),symbol=9, width=3, size=0.5;
    }
}



pdf, "Fig1_SystematicErrors";
