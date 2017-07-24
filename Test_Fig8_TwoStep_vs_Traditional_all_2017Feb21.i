#include "libCorrelationAlgorithmsGenerationV4.i"
require, "libImages.i";

//Make plot using: Fig8_TwoMethods_vs_Traditional.i

func write_two_step_vs_traditional(object)
/* DOCUMENT func write_two_step_vs_traditional(object)

   Writes the bias data for the four object cases 
   object is one of "point", "LGS", "GC", "solar".

   bias comparison: two step method vs traditional

   flux solar 5e4; flux for oether objects 5e3;
   SNR solar= 5e4/(sqrt(5e4)+16*16*1.0^2)=104 per image;
   SNR other objects = 5e3/(sqrt(5e3)+16*16*1.0^2) = 15 per image;

   TCoG selected for peakfinding because of its best performance against RMS noise
 */
{
  //setup object specific variables
  datafile= "Fits/Data_Two_step_traditional_"+object+".fits";

  if (object=="point") reference= image= MakePoint(8.0, 8.0);
  if (object=="LGS") reference= image= MakeLGS(8.0, 8.0);
  if (object=="GC") reference= image= MakeGC(8.0, 8.0);
  if (object=="solar") reference= image= MakeSolar(8.0, 8.0);

  Nx= dimsof(reference)(2);

  
  flux_i= 5e3;//total flux in image e-/sub-aperture | SNR = 5e3/(sqrt(5e3+16*16*1.0)) = 15.0;
  flux_solar=5e4; //solar flux
  flux_r= 16*16*flux_i;//total flux in reference = number of subaperture times flux_i
  RON= 1.0;//read out noise

   if(object=="solar") {
      reference= poisson(16*16*flux_r*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
  }else {
      reference= poisson(flux_r*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
  }
  image*= flux_i; //fluxify image 

  N= 40;//number of shifts
  NR= 500;//number of random realizations
  
  IN= TwoStep2=TwoStep1=Traditional=array(0.0, 4, N);
  CoG_2 = CoG_1 = array(0.0, 4, N);
  RTwoStep1=RTwoStep2=RTraditional= array(0.0, 2, NR);
  RCoG_2 = RCoG_1 = array(0.0, 2, NR);

  
  for(i=1; i<=N; i++){//for each shift
    
    sh= 2*double(i)/N-1.0;//sh moves from -1-ish to +1-ish pixel
    shx= sh;
    shy= sh;

    xc= Nx/2.0 + shx;
    yc= Nx/2.0 + shy;
    
    if (object=="point") image_sh= MakePoint(xc, yc);
    if (object=="LGS") image_sh= MakeLGS(xc, yc);
    if (object=="GC") image_sh= MakeGC(xc, yc);
    if (object=="solar") image_sh= MakeSolar(xc, yc);
    
    //currently flux of image_sh = 1 and before noisfying multiply with flux_i 
    if(object=="solar") {
        image_sh*= flux_solar;
        CorrlationMethod=7;
    }else if (object == "LGS") {
       image_sh*= flux_i*10; //5e4
       CorrlationMethod=1;
    }else {
        image_sh*= flux_i; //5e3
        CorrlationMethod=1;
    }
    
    for(j=1; j<=NR; j++){//for reach random realization
      image_noise= poisson(image_sh) + RON*random(dimsof(image_sh)); //noisify image
     
      //call IM_center to output external variable PeakIM an image with the
      //correlation peak
      reference1=reference;
      tmp= IM_center(reference1, image_noise, CorrlationMethod, 1, 1);
      RTraditional(,j)=CoGApp(PeakIM)(1:2);
      
      RTwoStep1 (,j)= SPcenter(reference1, image_noise, CorrlationMethod, 1, 5, 5)(1:2); //TCoG

      RTwoStep2 (,j)= Method1(reference1, image_noise, CorrlationMethod, 1, 5, 5)(1:2); //TCoG
      
      RCoG_1(,j)= CoGApp(image_noise);
      RCoG_2(,j)= TCoG2(image_noise, 1,5, 5)(1:2);
    }//NR
    
    IN(,i)= [xc, yc, xc, yc];//save input shifts
    
    Traditional(,i)= [ RTraditional(1,avg), RTraditional(2,avg),
                       RTraditional(1,rms)/sqrt(NR), RTraditional(2,rms)/sqrt(NR)]; //standard error on the mean
    
    TwoStep1(,i)= [ RTwoStep1(1,avg), RTwoStep1(2,avg),
                    RTwoStep1(1,rms)/sqrt(NR), RTwoStep1(2,rms)/sqrt(NR)]; //standard error on the mean
    
    TwoStep2(,i)= [ RTwoStep2(1,avg), RTwoStep2(2,avg),
                    RTwoStep2(1,rms)/sqrt(NR), RTwoStep2(2,rms)/sqrt(NR)];

    CoG_1(,i)= [RCoG_1(1,avg),RCoG_1(2,avg), RCoG_1(1,rms),RCoG_1(1,rms)];
    CoG_2(,i)= [RCoG_2(1,avg),RCoG_2(2,avg), RCoG_2(1,rms),RCoG_2(1,rms)];
                 
    
    write, i, sh, TwoStep2(,i)(1)-Nx/2., TwoStep2(,i)(2)-Nx/2., TwoStep2(,i)(2)-Nx/2.0-sh;
  }
  
  //initialize output data
  Data= array(0.0, N, 13);
  
  Data(,1)= IN(1,)-8; //x value
  
  // x value
  Data(,2)= Traditional(1,)-IN(1,);  
  Data(,3)= TwoStep1(1,)-IN(1,);
  Data(,4)= TwoStep2(1,)-IN(1,); 
  Data(,5)= CoG_1(1,)-IN(1,);
  Data(,6)= CoG_2(1,)-IN(1,);
  
  //error
  Data(,7)= Traditional(3,);
  Data(,8)= TwoStep1(3,);
  Data(,9)= TwoStep2(3,);
  Data(,10)= CoG_1(3,);
  Data(,11)= CoG_2(3,);
    
  //this file is used in Fig. 8
  fits_write, datafile, Data, overwrite=1;

  
  return Data;
}

// possible objects are "point", "LGS", "GC", "solar"

objects= ["point", "LGS", "GC", "solar"];

for(i=1;i<=2; i++){
    tmp=  write_two_step_vs_traditional( objects(i) );
}


//tmp=  write_two_step_vs_traditional("solar" );
