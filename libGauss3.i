require, "yao.i";
func GaussIM(Np,a)
/* DOCUMENT GaussIM(Np,a)

   returns a normalized (total flux == 1) gaussian image.
   
   x is an 2D array(0., Npix, Npix) with x and y (in pixel units)
   a is an array with:
   a(1) -- I0
   a(2) -- x center in pixels
   a(3) -- y center in pixels
   a(4) -- gaussian sigma
*/
{
    if(dimsof(a)(2) <4 || dimsof(a)(2) > 4) {
        print, " a  dimension should equal to 4:  help, GaussIM"; 
        1/0;
    }
    Np=int(Np);
    X= array(0., 2, Np, Np);
    X(1,)= span(1,Np,Np)(,-:1:Np);
    X(2,)= span(1,Np,Np)(-:1:Np,);

    result= a(1)*exp(-0.5* ((X(1,)-a(2))^2+(X(2,)-a(3))^2)/(a(4)^2));
    result/=2*pi*a(4)^2;
    return result;
}

func GaussIM2Sigma(Np,a)
/* DOCUMENT GaussIM(Np,a)

   returns a normalized (total flux == 1) gaussian image.
   
   x is an 2D array(0., Npix, Npix) with x and y (in pixel units)
   a is an array with:
   a(1) -- I0
   a(2) -- x center in pixels
   a(3) -- y center in pixels
   a(4) -- gaussian sigma
*/
{
    
    X= array(0., 2, Np, Np);
    X(1,)= span(1,Np,Np)(,-:1:Np);
    X(2,)= span(1,Np,Np)(-:1:Np,);

    result= a(1)*exp( -0.5* (   (X(1,)-a(2))^2/( 2*a(4)^2 ) +  (X(2,)-a(3))^2/ (2*a(5)^2)  ) ) ;
    result/=2*pi*a(4)*a(5);
    return result;
}


func GaussFunc(x,a)
/* DOCUMENT gaussfunc(x,a)
   The model function  to  fit (which creates the initial guess gaussian)
   
   Narsi 28/06/2012
*/
{
    X = x(1,)-a(2);
    Y=  x(2,)-a(3);
    ans= a(1)* exp(-0.5*(X^2+Y^2)/a(4)^2);
    ans = ans/(2*pi*a(4)^2);
    return ans;
}


func GaussFunc2d(x,a)
/* DOCUMENT gaussfunc(x,a)
   The model function  to  fit (which creates the initial guess gaussian)
   
   Narsi 28/06/2012
*/
{
    X = x(1,)-a(2);
    Y=  x(2,)-a(3);
    ans= a(1)* exp(-0.5*(X^2/ (2*a(4)^2)+ Y^2/(2*a(5)^2) ));
    ans = ans/(2*pi*a(4)*a(5));
    return ans;
}

/** -----------------------------oOo------------------------------------**/
func GaussFit(y, weight)
/* DOCUMENT GaussFit(y)
   
   fits gaussian for an input image data
   fma;pli, fitted_gauss
   
   Narsi 28/06/2012
  */
{
    extern fitted_gauss;

    if(y==[]) {
        print, "image is NULL, help, GaussFit"; 
        1/0;
    }
    
    Nx= dimsof(y)(2); //width
    Ny= dimsof(y)(3); //height
    
    x= array(0., 2, Nx, Ny); 
    x(1,)= (indgen(Nx)-Nx/2)(,-:1:Ny);
    x(2,)= (indgen(Ny)-Ny/2)(-:1:Nx,);
    
    sigma =compfwhm(y)/2.35;
    if(weight==[]) {
        w= 1.0;
    }else w=weight;
    
    B= median(y);
    pos=where2(y==max(y))-dimsof(y)(2)/2.0;
    
    a=[sum(y), pos(1), pos(2), sigma, sigma];
    // if (!silent) { write,"Fitting the Guassian ...";};
    result=lmfit(GaussFunc2d, x,a,y,w,stdev=1);
    fitted_gauss=GaussFunc2d(x,a);
    
    a(2) +=dimsof(y)(2)/2.;
    a(3) +=dimsof(y)(2)/2.;
    return [a, *result.stdev]; 
}

func Gauss1d(Np,a)
/* DOCUMENT GaussIM(Np,a)

   returns a normalized (total flux == 1) gaussian image.

   X= span(1,Np,Np);
   result= a(1)*exp(-0.5*((X-a(2))/a(3))^2);
   result/=2*pi*a(3)^2;
   
   X is an 1D array(0., Npix) with X (in pixel units)
   a is an array with:
   a(1) -- I0
   a(2) -- center in pixels
   a(4) -- gaussian sigma
*/
{
    X= span(1,Np,Np);
    
    result= a(1)*exp(-0.5*((X-a(2))/a(3))^2);
    result/=sqrt(2*pi)*a(3);
    return result;
}


func Eclipse1d(Np,a)
/* DOCUMENT Eclipse1d(Np,a)

   returns a normalized (total flux == 1) gaussian image.

   X= span(1,Np,Np);
   result= a(1)*x^2 + a(2)*xy + a(3)*y^2 + a(4)*x + a(5)*y + a(6);
   
   X is an 1D array(0., Npix) with X (in pixel units)
   
   X= span(1,Np,Np);
   x=X-x;
   y=Y-y;
   
   result= a(1)*x^2 + a(2)*xy + a(3)*y^2 + a(4)*x + a(5)*y + a(6);
   return result;
    
*/
{
    
          
    X= span(1,Np,Np);
    Y= transpose(span(1,Np,Np));
    x=X;
    y=Y;
    
    result= a(1)*x^2 + a(2)*x*y + a(3)*y^2 + a(4)*x + a(5)*y + a(6);
    return result;
}
