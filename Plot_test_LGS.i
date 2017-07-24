#include "libCorrelationAlgorithmsGenerationV4.i"
#include "libImages.i"

//make fig plot: Fig3_SubPixel_fitting_alg_performance_20160921.i

func Correlation_Peak(object)
/* DOCUMENT func CorrelationPeak(object)

   Correplation peak
   object is one of "point", "LGS", "GC", "solar".

   Compare bias error for sub-pixel algorithms.

   Flux for point, LGS and GC = 5e3/image
   for solar 5e4/image  (line 30 and 31)

   flux_i= 5e3;//total flux in image | SNR = 5e3/(sqrt(5e3+16*16*RON=1)) = 15.0;
   flux_solar=5e4; //solar image flux 
 
 */
{
  //setup object specific variables
  if (object=="point") reference= image= MakePoint(8.0, 8.0);
  if (object=="LGS") reference= image= MakeLGS(8.0, 8.0);
  if (object=="GC") reference= image= MakeGC(8.0, 8.0);
  if (object=="solar") reference= image= MakeSolar(8.0, 8.0);

  Nx= dimsof(reference)(2);

  
  flux_i= 5e3;//total flux in image | SNR = 5e3/(sqrt(5e3+16*16*RON=1)) = 15.0;
  flux_solar=5e4; //solar image flux 
  flux_r= 16*16*flux_i;//total flux in reference = number of subaperture times flux_i
  RON= 1.0;//read out noise

  if(object=="solar") {
      reference= poisson(16*16*flux_solar*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
  }else {
      reference= poisson(16*16*flux_i*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
  }
  
  image*= flux_i; //fluxify image 

  N= 10;//number of shifts
  NR= 500;//number of random realizations
  
  DATA= array(0.0, Nx, Nx, N);
  
  for(i=1; i<=N; i++){//for each shift
    
    sh= double(i)/10.0; 
    shx= sh;
    shy= sh;

    xc= Nx/2.0 + shx;
    yc= Nx/2.0 ;
    
    if (object=="point") image_sh= MakePoint(xc, yc);
    if (object=="LGS") image_sh= MakeLGS(xc, yc);
    if (object=="GC") image_sh= MakeGC(xc, yc);
    if (object=="solar") image_sh= MakeSolar(xc, yc);
    
    if(object=="solar") {
        image_sh*= flux_solar; //5e4
    }else{
        image_sh*= flux_i; // 5e3
    }
    
    image_noise= poisson(image_sh) + RON*random(dimsof(image_sh)); //noisify image
      
      //call IM_center to output external variable PeakIM an image with the
      //correlation peak
    reference1=reference;
    if(object=="solar"){
      tmp= IM_center(reference1, image_noise, 7, 1, 1);
    }else{
      tmp= IM_center(reference1, image_noise, 1, 1, 1);
    }
    DATA(,,i) = PeakIM;
    fma;pli, PeakIM;
    plg, array(8.0, 100), span(0, 16, 100) , color="red";
    plg, span(0, 16, 100) , array(8.0, 100), color="red";
    limits;
    f=swrite(format="CPeak_shift_X_%.1f_"+object, shx);
    png, f, dpi=300;
  }
  
  return DATA;
}


// possible objects are "point", "LGS", "GC", "solar"

objects= ["point", "LGS", "GC", "solar"];

for(i=1;i<=4; i++){
  //tmp= Correlation_Peak( objects(i) );
}

tmp= Correlation_Peak(  "LGS" );
