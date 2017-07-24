func CoG2(image, key, samp, CurveFit)
/* DOCUMENT func SPcenterCoG(image)

   does the SPcenter version for standerd CoG
 */
{
  Nx= dimsof(image)(2);//template X size
  P=FindPeak(image, CurveFit);
  
  if( int(samp)==1 ) return P;
 
  //check if IPcenter converged
  if(P(1) <= 3 || P(2) <= 3 || P(1) >= Nx || P(2) >= Nx )
    nerror("Got problem with func CoG()");

  //at frac =0.0 or shift ==0 (from CoG);
  Cx=P(1);
  Cy=P(2);
  
  //padding image
  PadrightIm= PadImage2x(image); //pad zeros double to its size

  for (k=1; k<=int(samp-1); k++){
    frac= double(k)/double(samp);

    //subpixel interpolation extraction
    temp= ExtractData(PadrightIm, key, Nx/2+1+frac, Nx+Nx/2+frac, Nx/2+1+frac, Nx+Nx/2+frac, Nx);
    
    //apply CoG algorithm
    Center=FindPeak(double(temp), CurveFit);
    grow, Cx, Center(1)+frac;
    grow, Cy, Center(2)+frac;
  }// end k loop
  
  return  [ [avg(Cx), avg(Cy)]]; //pixels
}
