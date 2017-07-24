require, "libCorrelationAlgorithmsGenerationV4.i";
require, "libImages.i";

/* used by these figs
Fig4_SNR_Point1.i
Fig6_SNR_GC1.i
Fig5_SNR_LGS1.i
Fig7_SNR_Solar1.i
*/

func write_SNR_data(object)
/* DOCUMENT func write_SNR_data(object)

   Writes the bias data for the four object cases 
   object is one of "point", "LGS", "GC", "solar".

   bias vs SNR
   bias vs sigma
   sigma vs SNR

   True shifts applied for point, GC and solar 0.3 pixel
                       for LGS 0.45 pixel

 */
{
  //setup object specific variables
  datafile= "Fits/Two_step_vs_SNR_"+object+".fits";

  if (object=="point")reference= image= MakePoint(8.0, 8.0);
  if (object=="LGS") reference= image= MakeLGS(8.0, 8.0);
  if (object=="GC") reference= image= MakeGC(8.0, 8.0);
  if (object=="solar") reference= image= MakeSolar(8.0, 8.0);

  Nx= dimsof(reference)(2);

  
  flux_i= 5e3;//total flux in image
  flux_solar=5e4; //solar flux per image
  flux_r= 16*16*flux_i;//total flux in reference = number of subaperture times flux_i
  RON= 1.0;//read out noise
  Interp=5;
  
  if(object=="solar") {
      reference= poisson(16*16*flux_solar*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
  }else {
      reference= poisson(flux_i*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
  }

  N= 40;//number of shifts
  NR= 500;//number of random realizations
  
  IN= Parabola= CoG= Equiangular= Pyramid= QuadraticInterpolationD=
    GaussData= SubPixelSpline= SubPixelBilinear= FFT= array(0.0, 4, N);
  
  RParabola= RCoG= REquiangular= RPyramid= RQuadraticInterpolationD= RGaussData= 
    RSubPixelSpline= RSubPixelBilinear= RFFT= array(0.0, 2, NR);

  ROCoG=array(0.0, 2, NR);
  OCoG=array(0.0, 4, N);
  
  
  VaryingFlux=spanl(50.0, 1e4, N); /* Flux varying point, LGS, GC*/
  SolarFlux=spanl(1e4, 1e7, N); /* Flux varying solar */
  
  SNR=array(0.0, N);/* computed SNR will be stored in here*/

  
  
  
  for(i=1; i<=N; i++){//for each shift
  
      if(object=="LGS"){
          shx= 0.4;
      }else if (object == "solar"){
          shx= 0.25;
        }else{
          shx=0.15;
        }
  
    xc= Nx/2.0 + shx;
    yc= Nx/2.0 + shx;
    
    if (object=="point") image_sh= MakePoint(xc, yc);
    if (object=="LGS") image_sh= MakeLGS(xc, yc);
    if (object=="GC") image_sh= MakeGC(xc, yc);
    if (object=="solar") image_sh= MakeSolar(xc, yc);

    
    /* Flux changed here */
    if(object=="solar") {
        CorrAlg=7;
        image_sh*= SolarFlux(i); //select solar flux
    }else{
        CorrAlg=1;
        image_sh*= VaryingFlux(i)//selects flux for point, LGS and GC objects
    }
    
    SNR(i)=sum(image_sh)/sqrt( sum(image_sh) + Nx*Nx*RON^2);//signal to noise ratio
     
    for(j=1; j<=NR; j++){//for reach random realization
      image_noise= poisson(image_sh) + RON*random(dimsof(image_sh)); //noisify image
      
      //call IM_center to output external variable PeakIM an image with the
      //correlation peak
      
      
      RParabola(,j)= SPcenter(reference, image_noise, CorrAlg, 1, Interp, 1)(1:2)- Nx/2.;
      RGaussData (,j)= SPcenter(reference, image_noise, CorrAlg, 1, Interp, 2)(1:2)- Nx/2.;
      RQuadraticInterpolationD(,j)= SPcenter(reference, image_noise, CorrAlg, 1, Interp, 3)(1:2)- Nx/2.;
      RPyramid(,j)= SPcenter(reference, image_noise, CorrAlg, 1, Interp, 4)(1:2)- Nx/2.;
      RCoG(,j)= SPcenter(reference, image_noise, CorrAlg, 1, Interp, 5)(1:2)- Nx/2.;
      ROCoG(,j)=
      //subtract Nx/2 to have 0 instead of 8 pix
      
/*
      window, 0;
      plp, RCoG(,j)(1)-sh, SNR(i);
*/

    }

    IN(,i)= [shx, shx, shx, shx];//save input shifts
    
    Parabola(,i)= [ RParabola(1,avg), RParabola(2,avg),
      RParabola(1,rms), RParabola(2,rms)]; //standard error on the mean
    
    Equiangular(,i)= [ REquiangular(1,avg), REquiangular(2,avg),
        REquiangular(1,rms), REquiangular(2,rms)];
    
    Pyramid(,i)= [ RPyramid(1,avg), RPyramid(2,avg),
                  RPyramid(1,rms), RPyramid(2,rms)];
    
    GaussData(,i)= [ RGaussData(1,avg), RGaussData(2,avg),
        RGaussData(1,rms), RGaussData(2,rms)];
    
    QuadraticInterpolationD(,i)= [ RQuadraticInterpolationD(1,avg), RQuadraticInterpolationD(2,avg),
        RQuadraticInterpolationD(1,rms), RQuadraticInterpolationD(2,)(rms) ];
    
    CoG(,i)= [ RCoG(1,avg),  RCoG(2,avg), RCoG(1,rms), RCoG(2,rms) ];
    OCoG(,i)= [ ROCoG(1,avg),  ROCoG(2,avg), ROCoG(1,rms), ROCoG(2,rms) ];

    write, i, shx, CoG(,i)(1)-shx;
  }
  
  //initialize output data
  Data= array(0.0, N, 13);  
  Data(,1)= SNR; //x value
  
  // x value
  Data(,2)= QuadraticInterpolationD(1,)-IN(1,);   //QPF
  Data(,3)= Parabola(1,)-IN(1,);                  //PF
  Data(,4)= OCoG(1,)-IN(1,);               //-- not used
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
    tmp=  write_SNR_data( objects(i) );
}


