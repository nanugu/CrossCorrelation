/*
  Provides functions that generates the apperture images.

  From oasis.i

  - GaussC2D
  - Gauss2D

 */

require, "convol.i";
require, "yao.i";

func GaussC2D(x,a,&pder,deriv=)
  /* DOCUMENT GaussC2D

     DESCRIPTION
        2D gaussian circle with amplitude a(1), center (a(2),a(3)),
	and radius a(4). x is of the form [[x1,y1],...].
   */

{
  xx=x(1,)-a(2);
  yy=x(2,)-a(3);

  e=exp(-.5*((xx^2+yy^2)/a(4)^2));
  y=a(1) * e;
  
  if (deriv) {
    pder=array(0.,dimsof(x)(numberof(dimsof(x))),4);
    pder(,1)=             e;
    pder(,2)= xx/a(4)^2*y;
    pder(,3)= yy/a(4)^2*y;
    pder(,4)= (xx^2+yy^2)/a(4)^3*y;
  }
  
  return y;
} 

func Gauss2D(x,a,&pder,deriv=)
  /* DOCUMENT Gauss2D

     DESCRIPTION
        2D gaussian elipse with amplitude a(1), center (a(2),a(3)),
	a axis a(4), b axis a(5) and theta=a(6).

	The elipse x^2/a^3+y^2/b^2 and in general the x,y axes are shifted
	by (a(2),a(3)) and rotated by theta=a(6).
   */
{
  xx=x(1,)-a(2);
  yy=x(2,)-a(3);

  xxr= cos(a(6))*xx+sin(a(6))*yy;
  yyr=-sin(a(6))*xx+cos(a(6))*yy;
  
  e=exp(-.5*(xxr^2/a(4)^2+yyr^2/a(5)^2));
  y=a(1) * e;
  
  if (deriv) {
    pder=array(0.,dimsof(x)(numberof(dimsof(x))),6);
    pder(,1)=            e;
    pder(,2)=-(-cos(a(6))*xxr/a(4)^2+sin(a(6))*yyr/a(5)^2)*y;
    pder(,3)= ( sin(a(6))*xxr/a(4)^2+cos(a(6))*yyr/a(5)^2)*y;
    pder(,4)= xxr^2/a(4)^3*y;
    pder(,5)= yyr^2/a(5)^3*y;
    pder(,6)=-(xxr*yyr/a(4)^2+xxr*yyr/a(5)^2)*y;
  }
  
  return y;
}

func MakePoint(xc,yc)
/* DOCUMENT func MakePoint(xc,yc)

   Makes a 16x16 image with a circular gaussian centred at position [xc,yc] and with
   fwhm = 2 pixel. Total flux in image is 1.   
*/
{
  x= double(indgen(16)(,-:1:16));
  y= transpose(x);

  a= [1.0, xc, yc, 2/2.35];

  image= GaussC2D(transpose([y,x]),a);
  
  return image/sum(image);
}

func MakeLGS(xc,yc)
/* DOCUMENT func MakeLGS(xc,yc)

   Makes a 16x16 image with an eliptical gaussian centred at position [xc,yc] and with
   fwhms = 2 pixel and 5 pixel. Total flux in image is 1.
*/
{
  x= double(indgen(16)(,-:1:16));
  y= transpose(x);

  a= [1.0, xc, yc, 3.0/2.35, 6./2.35, -pi/4];
  result=Gauss2D(transpose([y,x]),a);

  return result/sum(result);
}

func MakeGC(xc,yc)
/* DOCUMENT func MakeGC(xc,yc)

   Makes a 16x16 pix image, "centred" at position [xc,yc] and with 10 stars
   "randomly" located in the field.

   The stars are circular gaussians with fwhm = 2 pixel.
   Total flux in image is 1.
*/
{

  f= span(5.0,10.0,10);//fluxes of 10 stars
  d= span(0.0,6.0,10); //distances to centre
  q= spanl(0.1, 4*pi, 10); //position angle wrt centre
  image=array(0.,16,16);
  
  for (i=1;i<=10;i++){
    dx= d(i)*cos(q(i));
    dy= d(i)*sin(q(i));
    image+= f(i)*MakePoint(xc+dx,yc+dy);
  }

  return image/sum(image);
}

SolarImageOrig= double(fits_read("Fits/gband_22May2002_AR9957_1.fits"));


func MakeSolar(xc,yc)//old function
{
  extern SolarImageOrig; //for efficiency reasons is only read once from disk

  //scale= 100;//pixels of SolarImage used by one final pixel
  scale= 10;//pixels of SolarImage used by one final pixel

  /*
  nSolar=2000;//original image is 2010x2029, cropping to 2000x2000
  nSolarX= dimsof(SolarImageOrig)(2);
  nSolarY= dimsof(SolarImageOrig)(3);
  nSolar= min(nSolarX,nSolarY);
  dx= int(nSolarX-nSolar)/2;//xshift
  dy= int(nSolarY-nSolar)/2;//yshift
  SolarImage= SolarImageOrig(1+dx:dx+nSolar,1+dy:dy+nSolar);//crop to square
  */
  
  SolarImage=SolarImageOrig(250:451, 250:451);

  //make PSF
  //PSF should be odd (cf. convoln help)
  //we want to build a 2000x2000 image to crop to 1600x1600 and then to 16x16
  //nPSF= 1601;//PSF size
  nPSF= 161;//PSF size
  x= double(indgen(nPSF)(,-:1:nPSF));
  y= transpose(x);
  a= [1.0, double(xc)*scale, double(yc)*scale, 2.0*scale/2.35];
  PSF= GaussC2D(transpose([y,x]),a);

  //convolve and crop
  image= convoln(SolarImage, PSF);
/*
  //colapse and crop to 16x16 pix
  width= int(dimsof(image)(2)/20); //image is collapsed to 20x20
  result= array(0., 16, 16); //cropped result
  for (i=3; i<=18; i++){ //we crop 2 pixels in x and y
    for (j=3; j<=18; j++){
      result(i-2,j-2)= sum(image( (i-1)*width+1:i*width,  (j-1)*width+1:j*width ));
    }
  }
*/

  result=bin2d(image, 10);
  result=result(3:18, 3:18);
  return result/sum(result);
}

func RandomFlux(im, percentageMultiply){
    image=im;
    n=dimsof(image)(2);
    p=random(100000);
    p=p*2;
    list=where(  (p > (100.0-percentageMultiply)/100.0) & (p < (100.0+percentageMultiply)/100.0)  );
    P=p(list);
    Random=array(1.0, n, n);

    k=0;
    for(i=1; i<=n*n; i++){
        k++;
        Random(k)=P(k);
    }

    image=image * Random;
    return image;
}

#include "libCorrelationAlgorithmsGenerationV4.i"
func test(sh){
   
    Ref=MakeSolar(8,8);
    Tar=MakeSolar(sh, sh);
    
    write, "Test 1:  correlation type 1";
    X=IM_center(Ref, Tar, 1, 1,1);
    write, "Error=",ParabolaApp(PeakIM)(2)-sh;

    write, "Test 2:  correlation type 2";
    X=IM_center(Ref, Tar, 2, 1,1);
    write, "Error=",ParabolaApp(PeakIM)(2)-sh;

    write, "Test 3:  correlation type 3";
    X=IM_center(Ref, Tar, 3, 1,1);
    write, "Error=",ParabolaApp(PeakIM)(2)-sh;

    write, "Test 4:  correlation type 4";
    X=IM_center(Ref, Tar, 4, 1,1);
    write, "Error=",ParabolaApp(PeakIM)(2)-7.2;

    write, "Test 5:  correlation type 5";
    X=IM_center(Ref, Tar, 5, 1,1);
    write, "Error=",ParabolaApp(PeakIM)(2)-sh;

   
    write, "Test 6:  correlation type 7";
    X=IM_center(Ref, Tar, 7, 1,1);
    write, "Error=", ParabolaApp(PeakIM)(1:2)-sh;

    write, "Test 7:  Step 2 method";
    write, "Error=", SPcenter(Ref, Tar, 7, 1, 3, 5)(1:2)-sh;
     
    return 0;
}
p=test(7.1);

