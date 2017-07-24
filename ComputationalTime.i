require, "libImages.i";
require, "libCorrelationAlgorithmsGenerationV4.i";

/*
  Compute computational cost ratio

  steps:
  1. Compute convetional, Method1 and Method2 in a loop for 500 times.

  Measure time T1 for Convetuonal method
  Measure time T2 for Method1
  Measure time T3 for Method2

  print T2/T1, T3/T1

  And do it for 16x16 pixels
      and 32x32 pixels
 */


Ref=MakeLGS(8.,8.);
Temp=MakeLGS(8.45,8.45);

t0=t1=array(0.0, 3);
step=500;


/* 1. For conventional method */
timer, t0;//start time
for(i=1;i<=step; i++){
    a=IM_center(Ref, Temp, 1, 1, 1);
    a=ParabolaApp(PeakIM);
}
timer, t1; //end time
T1= (t1-t0);


/*2. Method 1 */
timer, t0; //start time
for(i=1;i<=step; i++){
    a=Method1(Ref, Temp, 1, 1, 3, 1);
}
timer, t1; //end time
T2= (t1-t0);

/* Method 2  */
timer, t0;
for(i=1;i<=step; i++){
    a=SPcenter(Ref, Temp, 1, 1, 3, 1);
}
timer, t1;
T3= (t1-t0);

write, "image (16x16 pixels), Method1/Conventional, Method2/Conventional = ", T2(1)/T1(1), T3(1)/T1(1);




//Creating a 32 pixel image by padding zeros
Ref=PadImage2x(Ref);
Temp=PadImage2x(Temp);

t0=t1=array(0.0, 3);

/* 1. For conventional method */
timer, t0;//start time
for(i=1;i<=step; i++){
    a=IM_center(Ref, Temp, 1, 1, 1);
    a=ParabolaApp(PeakIM);
}
timer, t1;//end time
T1= (t1-t0);


/*2. Method 1 */
timer, t0;//start time
for(i=1;i<=step; i++){
    a=Method1(Ref, Temp, 1, 1, 3, 1);
}
timer, t1;//end time
T2= (t1-t0);


/*2. Method 2 */
timer, t0;//start time
for(i=1;i<=step; i++){
    a=SPcenter(Ref, Temp, 1, 1, 3, 1);
}
timer, t1;//end time
T3= (t1-t0);


write, "image (32x32 pixels), Method1/Conventional, Method2/Conventional = ", T2(1)/T1(1), T3(1)/T1(1);
