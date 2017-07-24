#include "libCorrelationAlgorithmsGenerationV4.i"
#include "libImages.i"

//make fig plot: Fig3_SubPixel_fitting_alg_performance_20160921.i

func write_subpixel_data(object)
/* DOCUMENT func write_subpixel_data(object)

   Writes the bias data for the four object cases 
   object is one of "point", "LGS", "GC", "solar".

   Compare bias error for sub-pixel algorithms.

   Flux for point, LGS and GC = 5e3/image
   for solar 5e4/image  (line 30 and 31)

   flux_i= 5e3;//total flux in image | SNR = 5e3/(sqrt(5e3+16*16*RON=1)) = 15.0;
   flux_solar=5e4; //solar image flux 
   
 */
{
  //setup object specific variables
  datafile= "Fits/Data_performance_subpixel_"+object+"V2_5_5.fits";

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

  N= 40;//number of shifts
  NR= 500;//number of random realizations
  
  IN= Parabola= CoG= OCoG= Equiangular= Pyramid= QuadraticInterpolationD=
    GaussData= SubPixelSpline= SubPixelBilinear= FFT= array(0.0, 4, N);
  
  RParabola= RCoG= ROCoG= REquiangular= RPyramid= RQuadraticInterpolationD= RGaussData= 
    RSubPixelSpline= RSubPixelBilinear= RFFT= array(0.0, 2, NR);
  
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
    
    if(object=="solar") {
        image_sh*= flux_solar; //5e4
    }else{
        image_sh*= flux_i; // 5e3
    }
    
    
    for(j=1; j<=NR; j++){//for reach random realization
      image_noise= poisson(image_sh) + RON*random(dimsof(image_sh)); //noisify image
      
      //call IM_center to output external variable PeakIM an image with the
      //correlation peak
      reference1=reference;
      if(object=="solar"){
          tmp= IM_center(reference1, image_noise, 7, 1, 1);
      }else{
          tmp= IM_center(reference1, image_noise, 1, 1, 1);
      }
      
      RPyramid(,j)= PyramidApp(PeakIM) - Nx/2.;//subtract Nx/2 to have 0 instead of 8 pix
      RParabola(,j)= ParabolaApp(PeakIM) - Nx/2.;
      REquiangular(,j)= Equiangular_line(PeakIM) - Nx/2.;
      RGaussData(,j)= GaussApp(PeakIM) - Nx/2.;
      RQuadraticInterpolationD(,j)= QIApp(PeakIM) - Nx/2.;
      RCoG(,j)= CoGApp(PeakIM) - Nx/2.;
      ROCoG(,j)= TCoG2(image_noise, 1, 1, 5) - Nx/2.;
    }
    
    IN(,i)= [shx, shy, shx, shy];//save input shifts
    
    Parabola(,i)= [ RParabola(1,avg), RParabola(2,avg),
      RParabola(1,rms)/sqrt(NR), RParabola(2,rms)/sqrt(NR)]; //standard error on the mean
    
    Equiangular(,i)= [ REquiangular(1,avg), REquiangular(2,avg),
        REquiangular(1,rms)/sqrt(NR), REquiangular(2,rms)/sqrt(NR)];
    
    Pyramid(,i)= [ RPyramid(1,avg), RPyramid(2,avg),
                  RPyramid(1,rms)/sqrt(NR), RPyramid(2,rms)/sqrt(NR)];
    
    GaussData(,i)= [ RGaussData(1,avg), RGaussData(2,avg),
        RGaussData(1,rms)/sqrt(NR), RGaussData(2,rms)/sqrt(NR)];
    
    QuadraticInterpolationD(,i)= [ RQuadraticInterpolationD(1,avg), RQuadraticInterpolationD(2,avg),
        RQuadraticInterpolationD(1,rms)/sqrt(NR), RQuadraticInterpolationD(2,)(rms)/sqrt(NR) ];
    
    CoG(,i)= [ RCoG(1,avg),  RCoG(2,avg), RCoG(1,rms)/sqrt(NR), RCoG(2,rms)/sqrt(NR) ];

    OCoG(,i)= [ ROCoG(1,avg),  ROCoG(2,avg), ROCoG(1,rms)/sqrt(NR), ROCoG(2,rms)/sqrt(NR) ];

    write, i, sh, GaussApp(PeakIM)(1)-Nx/2., GaussApp(PeakIM)(1)-Nx/2.0-sh;
  }
  
  //initialize output data
  Data= array(0.0, N, 13);
  
  Data(,1)= IN(1,); //x value
  
  // x value
  Data(,2)= QuadraticInterpolationD(1,)-IN(1,);   //QPF
  Data(,3)= Parabola(1,)-IN(1,);                  //PF
  Data(,4)= OCoG(1,)-IN(1,);                      // CoGApp(image_noise)
  Data(,5)= Pyramid(1,)-IN(1,);                   //PYF
  Data(,6)= GaussData(1,)-IN(1,);                 //GF
  Data(,7)= CoG(1,)-IN(1,);                       //TCoG
  
  //error
  Data(,8)= QuadraticInterpolationD(3,);
  Data(,9)= Parabola(3,);
  Data(,10)= OCoG(3,);
  Data(,11)= Pyramid(3,);
  Data(,12)= GaussData(3,);
  Data(,13)= CoG(3,);
  
  //this file is used in Fig. 3
  fits_write, datafile, Data, overwrite=1;
  
  return Data;
}


// possible objects are "point", "LGS", "GC", "solar"

objects= ["point", "LGS", "GC", "solar"];

for(i=1;i<=2; i++){
  tmp= write_subpixel_data( objects(i) );
}

//tmp= write_subpixel_data( "LGS" );

