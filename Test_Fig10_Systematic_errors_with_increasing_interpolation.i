#include "libCorrelationAlgorithmsGenerationV4.i"
require, "libImages.i";


//used by Fig9_bias_vs_interpolation_performanceV4_plotV2.i

func write_bias_error_vs_interpolation(void)
/* DOCUMENT func write_subpixel_data(object)

   Writes the bias data vs interpolation
   for the four object cases object is
   one of "point", "LGS", "GC", "solar".

   Algorithm:
   1. Shifted point, GC and solar images with 0.3 px with their reference images. Whereas LGS shifted with 0.45 pixels.

   2. Measured back the input shifts with SPcenter using different interpolation.  
 */
{
  //setup object specific variables
  datafile= "Fits/Data_bias_vs_interpolation.fits";
  objects= ["point", "LGS", "GC", "solar"];
  
  flux_i= 5e3;//total flux in image
  flux_solar=5e4; //solar
  
  flux_r= 16*16*flux_i;//total flux in reference = number of subaperture times flux_i
  RON= 1.0;//read out noise
 
  N= 10; //10 interpolations
  NR= 500;//number of random realizations
  Nx=16;

  IN=array(0.0, N);
  Point=Point2=LGS=LGS2=GC=Solar= array(0.0, 4, N);
  RPoint=RLGS=RPoint2=RLGS2=RGC=RSolar= array(0.0, 2, NR);
  

  /* for each Object */
  for(OBJ=1; OBJ<=4; OBJ++){
      object=objects(OBJ);
      write, object;
      
      if(object=="LGS"){
      LGSxc=Nx/2.0 + 0.4;
      LGSyc=Nx/2.0 + 0.4;
      } else if(object == "solar"){
        Solarxc = Nx/2.0 + 0.26;
        Solaryc = Nx/2.0 + 0.26;
      } else {
        xc= Nx/2.0 + 0.15;
        yc= Nx/2.0 + 0.15;
      }
          
      
      
      
      /*Prepare reference and a shifted image */
      if (object=="point") {
          reference= MakePoint(8.0, 8.0); 
          image_sh= MakePoint(xc, yc);
      }
      
      if (object=="LGS") {
          reference= MakeLGS(8.0, 8.0);
          image_sh= MakeLGS(LGSxc, LGSyc);
          image_sh_lgs= MakeLGS(xc+0.15,yc+0.15);
          
      }
      if (object=="GC") {
          reference= MakeGC(8.0, 8.0);
          image_sh= MakeGC(2*xc, 2*yc);
      }

      //for any of the point, GC and LGS
      reference= poisson(16.0*16.0*flux_i*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
      image_sh*= flux_i; //flux for point, GC and LGS
      if(object == "LGS")image_sh_lgs*= flux_i*10; //flux for point, GC and LGS
      
      if (object=="solar") {
          reference= MakeSolar(8.0, 8.0);
          image_sh= MakeSolar(Solarxc, Solaryc);

          reference= poisson(16.0*16.0*flux_solar*reference)+RON*random(dimsof(reference)); //fluxify and noisify reference
          image_sh*= flux_solar; //flux for solar
      }
    

      /*For different interpolation */
  for(interpol=1; interpol<=N; interpol++){ //interpolation
      IN(interpol)=interpol;
      
      /* For different random noise read */
      for(Random=1; Random<=NR; Random++){//for reach random realization
          image_noise= poisson(image_sh) + RON*random(dimsof(image_sh)); //noisify image
          image_noise_lgs= poisson(image_sh_lgs) + RON*random(dimsof(image_sh)); //noisify image
          
          //call IM_center to output external variable PeakIM an image with the
          //correlation peak
          reference1=reference;
          if (object=="point") {
            RPoint(, Random)= SPcenter(reference1, image_noise, 1, 1, interpol, 5)(1:2);
            RPoint2(, Random)= TCoG2(image_noise, 1, interpol, 5)(1:2);
          } //CoG
          
          if (object=="LGS"){
            RLGS(, Random)= SPcenter(reference1, image_noise, 1, 1, interpol, 5)(1:2);
            RLGS2(, Random)= TCoG2(image_noise_lgs, 1, interpol, 5)(1:2);
          }//CoG
         if (object=="GC") RGC(, Random)= SPcenter(reference1, image_noise, 1, 1, interpol, 5)(1:2); //CoG
         if (object=="solar") RSolar(, Random)= SPcenter(reference1, image_noise, 7, 1, interpol, 5)(1:2); //CoG
      }//Random

      if (object=="point"){
        Point(,interpol)= [ RPoint(1,avg), RPoint(2,avg),
                            RPoint(1,rms), RPoint(2,rms)]; //standard error on the mean
        Point2(,interpol)= [ RPoint2(1,avg), RPoint2(2,avg),
                            RPoint2(1,rms), RPoint2(2,rms)]; //standard error on the mean
      }
      
      if (object=="LGS") {
        LGS(,interpol)= [ RLGS(1,avg), RLGS(2,avg),
                          RLGS(1,rms), RLGS(2,rms)]; //standard error on the mean
        LGS2(,interpol)= [ RLGS2(1,avg), RLGS2(2,avg),
                          RLGS2(1,rms), RLGS2(2,rms)]; //standard error on the mean
      }
      
      if (object=="GC") GC(,interpol)= [ RGC(1,avg), RGC(2,avg),
                                         RGC(1,rms), RGC(2,rms)]; //standard error on the mean
      
      if (object=="solar") Solar(,interpol)= [ RSolar(1,avg), RSolar(2,avg),
                                               RSolar(1,rms), RSolar(2,rms)]; //standard error on the mean
     
      write, interpol;
  }//interpolation

  
  } //OBJ
  
  //initialize output data
  Data= array(0.0, N, 11);
  Data(,1)= IN; //x value
  
  // x value
  Data(,2)= (abs(Point(1,)-xc) + abs(Point(2,)-xc))/2.0;
  Data(,3)= (abs(LGS(1,)-LGSxc) + abs(LGS(2,)-LGSxc) )/2.0;
  Data(,4)= ( abs(GC(1,)-xc)  + abs(GC(2,)-xc) )/2.0;
  Data(,5)= ( abs(Solar(1,)-Solarxc) + abs(Solar(2,)-Solarxc)  )/2.0;
  Data(,6)= (abs( Point2(1,)-xc) + abs( Point2(2,)-xc) )/2.0;
  Data(,7)= (abs(LGS2(1,)-xc-0.15) + abs(LGS2(2,)-xc-0.15) )/2.0;
  
  //error
  Data(,8)= Point(3,);  
  Data(,9)= LGS(3,);
  Data(,10)= GC(3,);
  Data(,11)= Solar(3,);
  
 
  //this file is used in Fig. 7
  fits_write, datafile, Data, overwrite=1;
  return Data;
}


Data=write_bias_error_vs_interpolation(DO);


write, "Improvement in percentage for point, LGS, GC and solar =", (Data(,2)(1)-Data(,2)(5))/Data(,2)(1)*100, (Data(,3)(1)-Data(,3)(5))/Data(,3)(1)*100, (Data(,4)(1)-Data(,4)(5))/Data(,4)(1)*100, (Data(,5)(1)-Data(,5)(5))/Data(,5)(1)*100, (Data(,6)(1)-Data(,6)(5))/Data(,6)(1)*100, (Data(,7)(1)-Data(,7)(5))/Data(,7)(1)*100;

write, "Improvement in factor for point, LGS, GC and solar =", (Data(,2)(1))/Data(,2)(5), (Data(,3)(1))/Data(,3)(5), (Data(,4)(1))/Data(,4)(5), (Data(,5)(1))/Data(,5)(5), Data(,6)(1)/Data(,6)(5), Data(,7)(1)/Data(,7)(5);



