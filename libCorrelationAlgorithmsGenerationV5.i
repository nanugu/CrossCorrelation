require, "yao.i";
require, "libGauss3.i";
require, "libSubPixelAlgV1.i";

func ShiftImage(image, x, y)
/* DOCUMENT ShiftImage(image, x, y)

   shifts an image with yao fftshift
   the image is previously padded to twice its size.
 */
{
    temp=PadImage2x(image);
    shifted=fftshift(temp, double(x), double(y));
    return UnPadImage2x(shifted);
}

func PadImage2x(im)
/* DOCUMENT PadImage2x(im)

   takes an image (with even number of pixel) and doubles its
   size padding with zeroes.

 */
{
    Nx=dimsof(im)(2);
    Ny=dimsof(im)(3);
    if (Nx != Ny) nerror("PadImage2x: image must be square");
    left= Nx/2;
    win= Nx;

    Padded=array(0.0, 2*Nx, 2*Ny);
    Padded(left+1:left+win, left+1:left+win)= im;
    return Padded;
}

func UnPadImage2x(Padded)
/* DOCUMENT UnPadImage2x(im)

   takes an even images previously padded with zeros and returns center.
   
   N.B.: if Nx in the size of padded image; unpadded image starts at  Nx/2-Nx/4.
*/
{
    Nx=dimsof(Padded)(2);
    Ny=dimsof(Padded)(3);
    if (Nx != Ny) nerror("UnPadImage2x: image must be square");
    left= Nx/2;
    win= Nx/2;
    
    return Padded(left-win/2+1:left+win/2, left-win/2+1:left+win/2);
}

func ExtractData(image, key, Is, Ie, Js,Je, sampling)
/* DOCUMENT ExtractData(image, key, Is, Ie, Js, Je, sampling)
   
   key==1 -> bilinear
   key==2 -> interp2
   key==3 -> spline2
   else  -> akima2
   
   image(Is:Ie,Js:Je) -> how many pixel samples
   ex: image(8.1:24.2, 8.1:24.2) to 16 pixels
*/
{
    result= image;
    sampling= int(sampling);
    x0 = span(Is, Ie, sampling)(,-:1:sampling);
    y0 = span(Js, Je, sampling)(-:1:sampling,);

    if(key == 1){
      
        result=bilinear(image, x0, y0);

    }else if (key==2){

      sizeX = dimsof(image)(2);
      X= array(0., 2, sizeX, sizeX);
      X(1,)= span(1,sizeX,sizeX)(,-:1:sizeX);
      X(2,)= span(1,sizeX,sizeX)(-:1:sizeX,);
      result= interp2(y0,x0,image,X(2,),X(1,));

    }else if(key==3){

      result= double(spline2(float(image), float(x0), float(y0)));

    } else nerror("key is an int from 1 to 4.");
    
    return result;
}

func CrossCorrFFT(template, image)
/* DOCUMENT CrossCorrFFT(template, image)
   template -- template image, 2d, real, squared image
   image  -- live image, 2d, real, same size as template

   returns normalised cross correlated image
 */
{
    //ZP is zero padded
    ZPtemplate=  PadImage2x(template);
    ZPimage=  PadImage2x(image);
    
    dims = dimsof(ZPtemplate);
    ZPtemplate = fft(ZPtemplate, +1);
    ZPimage = fft(ZPimage, +1);
    
    CrossC = fft(ZPimage*conj(ZPtemplate), -1);
    CrossC = double(CrossC);//take real part
    offset= dims(2:3)/2 - 1;//put center at left to find peak at 0
    CrossC= roll(CrossC, offset);
    
    Peak= UnPadImage2x(CrossC);//unpadding the cross correlation peak
    denom = sqrt(sum(template^2) * sum(image^2));

    return  Peak/denom;
}

func CrossCorrFFT1(template, image)
/* DOCUMENT CrossCorrFFT(template, image)
   template -- template image, 2d, real, squared image
   image  -- live image, 2d, real, same size as template

   returns normalised cross correlated image
 */
{
    //ZP is zero padded
    ZPtemplate=  template;
    ZPimage=  image;
    
    dims = dimsof(ZPtemplate);
    ZPtemplate = fft(ZPtemplate, +1);
    ZPimage = fft(ZPimage, +1);
    
    CrossC = fft(ZPimage*conj(ZPtemplate), -1);
    CrossC = double(CrossC);//take real part
    offset= dims(2:3)/2 - 1;//put center at left to find peak at 0
    CrossC= roll(CrossC, offset);
    
    Peak= CrossC;//unpadding the cross correlation peak
    denom = sqrt(sum(template^2) * sum(image^2));

    return  Peak/denom;
}

func CorrelationAlgorithm(image, template, Flag)
/* DOCUMENT CorrelationAlgorithm(image, template, Flag)

   Uses image and template to compute a "correlation" accordig to Flag
   
   Flag== 1 -> cross correlation
   Flag== 2 -> SSD
   Flag== 3 -> SAD
   Flag== 4 -> ZSAD
   Flag== 5 -> ZSSD
   Flag== 6 -> ZASSD

   cf. Lofdhal 2010, AA, 524, A90 and notes
*/
{
  TINY= 1e-50;

  if (Flag==1) { result = sum( image*template ); //cross correlation

  } else if (Flag==2) {//SSD

    result= 1.0/(sum( (image-template)^2 ) + TINY);

  } else if (Flag==3) { //SAD

    result= 1.0/(sum( abs(image-template) ) + TINY);

  } else if (Flag==4){//ZSAD

    result= 1.0/( sum( abs( (image-avg(image)) - (template-avg(template)) ) )
                  + TINY);

  } else if (Flag==5){//ZSSD

    result= 1.0/( sum( ( (image-avg(image)) - (template-avg(template)) )^2 )
                  + TINY);

  } else if (Flag==6) {//ZASSD

    result= 1.0/( sum(abs(image-template))^2 + TINY);

  }  else if (Flag==7){
      result = sum( image*template );
  }
  else nerror("Flag is a int from 1 to 6.")

  return result;
}

func FindPeak(image, CurveFit)
/* DOCUMENT FindPeak(image, CurveFit)

   returns the position of the maximum of an image using a variety
   of algorithms defined in flag CurveFit

   CurveFit == 1 -> ParabolaApp
   CurveFit == 2 -> GaussFit
   CurveFit == 3 -> QI
   CurveFit == 4 -> Pyramid
   CurveFit == 5 -> CoGApp
*/
{
    if (CurveFit ==1 ){

      C= ParabolaApp(image);

    } else if(CurveFit == 2){

      C= GaussApp(image);

    } else if (CurveFit == 3){
      
      C= QIApp(image);

    } else if (CurveFit == 4){

      C= PyramidApp(image);

    } else if (CurveFit==5) {

      C=CoGApp(image);

    } else nerror("CurveFit is integer from 1 to 5");

    return C;
}

func IM_center(template, image,  Flag, key, samp)
/* DOCUMENT IM_center(template, image,  Flag, key, samp)
   template-> template image
   image -> to match
   both the images should be the same size and squared
   returns the correlated position array

   Flag==1 -> cross correlation
   Flag==2 -> SSD
   Flag==3- > SAD
   Flag==4 -> ZSAD
   Flag==5 -> ZSSD
   Flag==6 -> ZASSD
   
   Flag=7 -> FFT correlation
   Flag=8 -> GaussFit
 
   samp-> the number of smaple points b/w two pixels
          to make interpolation

   The method is not fast when samp > 1;
   See: help,  SPcenter
          
   returns the centroids in pixels
 */
{
  extern PeakIM;
  
  image=double(image);
  template=double(template);
  if(Flag==7) {
      image=image-avg(image);
      template=template-avg(template);
  }
 
  if( key != 1) nerror("key is 1 only!");
  if( samp != 1) nerror("samp is 1 only!");
  
  
  //image=sum(template)*image/sum(image);
  tempX = dimsof(template)(2); //X size of template image
  win =tempX/2;

  PadrightIm= PadImage2x(image);
  PadX= dimsof(PadrightIm)(2);
  PadY= dimsof(PadrightIm)(3);
  
  
  Corr= array(0.0, dimsof(PadrightIm));
  
  for(i= win+1; i<= PadX-(win+1); i++){//scan image from left to right
    
    for(j= win+1; j<= PadY-(win+1); j++){//idem from bottom to top
      
      temp=PadrightIm(i-win+1:i+win, j-win+1:j+win);
      
      Corr(i,j)= CorrelationAlgorithm(temp, template, Flag);
      
    }
  }
  
  PeakIM= UnPadImage2x( Corr );
  if (Flag==7) {
      if(abs(min(PeakIM)) > abs(max(PeakIM))) PeakIM=-PeakIM;
  }
  
  C= where2( PeakIM == max(PeakIM) );
  
  return  [ [C(1), C(2)], compfwhm(Corr)]; //pixels
}

func IPcenter(template, image,  Flag, key, samp, CurveFit)
/* DOCUMENT IP_center(template, image,  Flag, key, samp)
   template-> template image
   image -> to match
   both the images should be the same size and squared
   returns the correlated position array

   Flag==1 -> cross corrlation
   Flag==2 -> SSD
   Flag==3 -> SAD
   Flag==4 -> ZSAD
   Flag==5 -> ZSSD
   Flag==6 -> ZASSD
   Flag==7 -> FFT correlation
   Flag==8 -> GaussFti
 
   key==1 -> bilinear interpolation
   key==2 -> interp2
   key==3 -> spline2
   else  -> akima2

   CurveFit == 1 -> ParabolaApp
   CurveFit == 2 -> GaussFit
   CurveFit == 3 -> QI

   
   samp-> the number of smaple points b/w two pixels
          to make interpolation
          
   returns the centroids in pixels
 */
{
  extern PeakCI;
  
  image=double(image);
  template=double(template);
  
  extern PeakFFT_;

  //sandard approach no subsampling
  if (Flag == 7) {
    
    Peak = CrossCorrFFT(template, image);
    PeakFFT_=Peak;

    C= FindPeak(Peak, CurveFit);
    
    return [[C(1), C(2)], [compfwhm(Peak), compfwhm(Peak)]];
  }
    
  if (Flag == 8) { //idem
    
    G=GaussFit(image);
    
    return [ [G(2), G(3)], [G(4)*2.355, G(5)*2.355]]; //pixels
  }
  
  image=sum(template)*image/sum(image);
  tempX = dimsof(template)(2);//size of template in X
  win =tempX/2;
  
  PadrightIm= PadImage2x( image ); //pad zeros double to its size
  lefti= dimsof(PadrightIm)(2)/4+1;
  leftj= dimsof(PadrightIm)(3)/4+1;
  
  //setup cube fo each sampling
  Corr= array(0.0, dimsof(PadrightIm))(,,-:1:int(samp));
  
  for(i=lefti; i<=lefti+(2*win-1); i++){ //scanning from left to right
    
      for(j=leftj; j<=leftj+(2*win-1); j++){
        
        //interpolation  k
        for(k=1; k<=int(samp); k++){
          
          if(samp == 1 || k == 1) {
            
            temp= PadrightIm(i-win+1:i+win, j-win+1:j+win); //extract scan image
            
          } else {
            
            frac= (k-1.0)/double(samp);
            
            samp=double(samp);
            
            temp= ExtractData(PadrightIm, key, i-win+1+frac, i+win+frac,\
                               j-win+1+frac, j+win+frac, tempX);
            temp= double(temp); //Interpolated data output
            
          }
          
          Corr(i,j,k)= CorrelationAlgorithm(temp, template, Flag);
          
        }//for k
        
      }
    }
    
    Cx=0.0;
    Cy=0.0;

    for(k=1; k<=int(samp); k++){
      
      frac= (k-1.0)/double(samp);
      
      temp= UnPadImage2x(Corr(,,k));
      temp/= max(temp);

      if (Flag==7) {
          if(abs(min(temp)) > abs(max(temp))) temp=-temp;
      }
      
      C= FindPeak(temp, CurveFit) + frac;
      
      Cx+= C(1);
      Cy+= C(2);   
    }
    
    return  [ [Cx/double(samp), Cy/double(samp)], compfwhm(Corr(,,1))]; //pixels
}


func IP_center(template, image,  Flag, key, samp)
/* DOCUMENT IP_center(template, image,  Flag, key, samp)

   the same as IPcenter but with CurveFit=1 (ParabolaApp)
   cf. IPcenter
*/
{
  return IPcenter(template, image, Flag, key, samp, 1);
}


func SPcenter(template, image,  Flag, key, samp, CurveFit)
/* DOCUMENT SPcenter(template, image,  Flag, key, samp, CurveFit)
   template-> template image
   image -> to match
   both the images should be the same size and squared
   returns the correlated position array

   Flag==1 -> cross corrlation
   Flag==2-> SSD
   Flag==3-> SAD
   Flag==4-> ZSAD
   Flag==5 -> ZSSD
   Flag==6 -> ZASSD
   Flag=7-> CFI, see  Lofdahl 2010 formulas
   Flag=8-> FFT correlation

   key==1 -> bilinear interpolation
   key==2 -> interp2
   key==3 -> spline2
   else  -> akima2

   CurveFit ==1 -> ParabolaApp
   CurveFit ==2 -> GaussFit
   CurveFit ==3 -> QI
   CurveFit ==4 -> Pyramid
   CurveFit ==5 -> CoGApp
  
   
   
   samp-> the number of sample points between two pixels
          to make interpolation
          
   returns the centroids in pixels
 */
{
  extern PeakFFT, PeakC;
  image=double(image);
  template=double(template);
  
  if(Flag==7){
  image=image-avg(image);
  template=template-avg(template);
  }
  image=image*max(template)/max(image);
  tempX = dimsof(template)(2);//template X size
  
  if (Flag == 8) { //standard approach
    
    PeakFFT= CrossCorrFFT(template, image);
    
    
    C= FindPeak(PeakFFT, CurveFit);
    
    //save centre and FWHM
    return [[C(1), C(2)], [compfwhm(PeakFFT), compfwhm(PeakFFT)]];
  }
  
  
  if(Flag == 8) {//idem

    G=GaussFit(image);

    //save centre and FWHM
    return [ [G(2), G(3)], [G(4)*2.355, G(5)*2.355]]; //pixels
  }
  
  //first iteraction
  //obtain integer shift cross correlation centre 
  P=IM_center(template, image,  Flag, 1, 1);
  P=FindPeak(PeakIM, CurveFit);
  
  //check if IPcenter converged
  if(P(1) <= 3 || P(2) <= 3 || P(1) >= tempX || P(2) >=tempX )
    nerror("Got problem with func IPcenter()");
  
  if( int(samp)==1 ) return P;

  //at frac =0.0 or shift ==0 (from IPcenter);
  Cx=P(1);
  Cy=P(2);
  
  
  win= tempX/2;
  
  //position of peak in the padded image
  PX= int( P(1) + tempX/2.0 );
  PY= int( P(2) + tempX/2.0 );

  //padding image
  PadrightIm= PadImage2x(image); //pad zeros double to its size
  
  //set up cube for subsampling
  Corr= array(0.0, dimsof(PadrightIm))(,,-:1:int(samp-1));

  //scan only in the 5 pixel window around maximum of first iteraction
  for (i= PX-2; i<= PX+2; i++){
    //define window X limits in scan
    leftX= i-win+1;
    rightX= i+win;

    for (j= PY-2; j<= PY+2; j++){
        //define window Y limits in scan
      leftY= j-win+1;
      rightY= j+win;
      
      //interpolation  k
      for (k=1; k<=int(samp-1); k++){
          frac= double(k)/double(samp);
          //subpixel interpolation extraction
          temp= ExtractData(PadrightIm, key, leftX+frac,rightX+frac,
                            leftY+frac,rightY+frac, tempX);
          temp = double(temp);
          
          
          //apply correlation algorithm
          Corr(i,j,k)= CorrelationAlgorithm(temp, template, Flag);
          
      }// end k loop
    }//j
  }//i

 
  
  PeakC=Corr;
  
  leftX= leftY= win+1;
  rightX= rightY= 3*win;

  for (k=1; k<=int(samp-1); k++){
      frac= double(k)/double(samp); //put back fractional sampling
      
      temp= Corr(leftX:rightX, leftY:rightY, k);
      
      if (Flag==7) {
          if(abs(min(temp)) > abs(max(temp))) temp=-temp;
      }
      
      C= FindPeak(temp, CurveFit) + frac;
      
      grow, Cx, C(1);
      grow, Cy, C(2);   
  }

  return  [ [avg(Cx), avg(Cy)], compfwhm(Corr(,,1))]; //pixels
}


func PlotOverPeak(Peak, title){
    window, 1;
    fma; pli, Peak;
    plg, span(0, 16, 20), array(8, 20);
    plg, array(8, 20), span(0, 16, 20);
    pltitle, title;
    return 0;
}

func padzero(image){
    im=image;
    nx=dimsof(im)(2);
    n=2;

    im(1:n, :)=0;
    im(:, 1:n)=0;
    im(nx-n+1:nx, :)=0;
    im(:, nx-n+1:nx)=0;

    return im;
}


func Method1(Template, Target, Flag, interpolation, samp, PeakAlg)
/* DOCUMENT help, SPcenter; input parameters are same as SPcenter function.

   It uses IM_center function

   It measures image shift at 1 pixel grid
   and sub-pixel grid
 */
{

    if(Flag==7){
    Template=Template-avg(Template);
    Target=Target-avg(Target);
    }
    
    Cx=Cy=[];
 
    Nx=dimsof(Template)(2);

    //compute image shift at 1 pixel grid
    Peak=IM_center(Template, Target, Flag, interpolation, 1);
    C= FindPeak(PeakIM, PeakAlg);
    grow, Cx, C(1);
    grow, Cy, C(2);

    PadrightIm= PadImage2x(Target);

    //compute image shift at sub-pixel grid
    //shift image at 1/samp and compute cross-correlation
    for(i=2; i<=samp; i++){
        frac=1/double(samp)*(i-1);

        temp= ExtractData(PadrightIm, 1, Nx/2+1-frac, Nx/2+Nx-frac, Nx/2+1-frac, Nx/2+Nx-frac, Nx);

        /* Obtain cross correlation */
        Peak=IM_center(Template, temp, Flag, interpolation, 1);
        /* Apply peak finding algorithm*/
        C= FindPeak(PeakIM, PeakAlg)-frac;
        
        grow, Cx, C(1);
        grow, Cy, C(2);   
    }

    /*Return average of all image shifts (they are anti-symmetric) which reduces the systematic error*/
    return [ [avg(Cx), avg(Cy)] ]; //pixels
}

func Center(Template, Target, Flag, interpolation, samp, PeakAlg)
/* DOCUMENT help, SPcenter; input parameters are same as SPcenter function.

   It uses IM_center function

   It measures image shift at 1 pixel grid
   and sub-pixel grid
 */
{

    if(Flag==7){
    Template=Template-avg(Template);
    Target=Target-avg(Target);
    }
    
 
    Nx=dimsof(Template)(2);
    Cx=Cy=[];
    
    //compute image shift at 1 pixel grid
    Peak=IM_center(Template, Target, Flag, interpolation, 1);
    C= FindPeak(PeakIM, PeakAlg);
    if(samp==1) return C;
    
    grow, Cx, C(1);
    grow, Cy, C(2);   
   
    
    PadrightIm= PadImage2x(Target);
    Zero=0;
    Fracarray=span(-0.45, 0.45, int(samp));

    if( odd(int(samp)) ) {
        Zero=int(samp/2.0+1);
    }
    
    
    
    //compute image shift at sub-pixel grid
    //shift image at 1/samp and compute cross-correlation
    for(i=1; i<=int(samp); i++){

        
        if(Zero != i){//forget the zero shift frac
        frac=Fracarray(i);
        /* Shift image using interpolation */
        temp= ExtractData(PadrightIm, 1, Nx/2+1+frac, Nx/2+Nx+frac, Nx/2+1+frac, Nx/2+Nx+frac, Nx);

        /* Obtain cross correlation */
        Peak=IM_center(Template, temp, Flag, interpolation, 1);
        /* Apply peak finding algorithm*/
        C= FindPeak(PeakIM, PeakAlg)+frac;
        
        grow, Cx, C(1);
        grow, Cy, C(2);   
        }
    }
    /*Return average of all image shifts (they are anti-symmetric) which reduces the systematic error*/
    return [ [avg(Cx), avg(Cy)] ]; //pixels
}


func CoG_step2(image, key, samp, G_sigma)
/* DOCUMENT func SPcenterCoG(image)

   does the SPcenter version for standerd CoG
 */
{
  Nx= dimsof(image)(2);//template X size
  P=WCoGXY_Gauss(image, G_sigma);
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
    Center=WCoGXY_Gauss(double(temp));
    grow, Cx, Center(1)+frac;
    grow, Cy, Center(2)+frac;
  }// end k loop
  
  return  [ [avg(Cx), avg(Cy)], [compfwhm(image), compfwhm(image)] ]; //pixels
}
