#include "libGauss3.i"
func ParabolaApp(image)
/* DOCUMENT func ParabolaApp(image)

   Measures the center of the an image with sub-pixel parabola fitting.
   
   Sub pixel accuracy measurement with Parabola approximation
   Sub-pixel estimation of local estrema,

   DG Bailey, Image and Vision Computing NZ
   
   Eq No. 4
*/
{
    image=image/max(image);
  IndicesMax = where2(image == max(image));
  im =  IndicesMax(1);
  jm =  IndicesMax(2); 
  
  size=dimsof(image)(2);
    
  if(im==1){
    im1= size; //im-1
    Xm = im + 0.5*(image(im1,jm)- image(im+1,jm))/\
      ( image(im1,jm) + image(im+1,jm) -2*image(im,jm));
  }

  if (jm==1){
    jm1=size; //jm-1
    Ym = jm + 0.5*( image(im,jm1)- image(im,jm+1) )/\
      ( image(im,jm1) + image(im,jm+1) -2*image(im,jm));
  }

  if (jm==size){
    jm16=1;  //jm+1
    Ym = jm + 0.5*( image(im,jm-1)- image(im,jm16) )/\
      ( image(im,jm-1) + image(im,jm16) -2*image(im,jm));
  }

  if (im==size){
    im16=1; //jm+1
    Xm = im + 0.5*(image(im-1,jm)- image(im16,jm))/\
      ( image(im-1,jm) + image(im16,jm) -2*image(im,jm));
  }

  if (im !=1 && im !=size){
    Xm = im + 0.5*(image(im-1,jm)- image(im+1,jm))/\
      ( image(im-1,jm) + image(im+1,jm) -2*image(im,jm));
    
  }
  
  if (jm!=1 && jm !=size){
    Ym = jm + 0.5*( image(im,jm-1)- image(im,jm+1) )/           \
      ( image(im,jm-1) + image(im,jm+1) -2*image(im,jm));
  }
  
  return [Xm, Ym];
}


func QIApp(image)
/* DOCUMENT func Polynomial2Order(image)

   M. G. Löfdahl: Evaluation of image-shift measurement algorithms for solar SH WFS, A&A, (2010)
   Adapted equation from Table 2.
*/
{

   image=image/max(image);
    IndicesMax = where2(image == max(image));
    im =  IndicesMax(1);
    jm =  IndicesMax(2);

    a2 = ( image(im+1,jm)-image(im-1,jm) )/2.;
    a3 =  ( image(im+1,jm)-2*image(im,jm)+ image(im-1,jm) )/2.;
    a4 = ( image(im,jm+1)-image(im,jm-1) )/2.;
    a5 = ( image(im,jm+1) -2*image(im,jm)+ image(im,jm-1) )/2.;
    a6 = ( image(im+1,jm+1) - image(im-1,jm+1) - image(im+1,jm-1) + image(im-1,jm-1))/4.;
  
    Xm = im - a2/(2*a3);
    Ym = jm - a4/(2*a5);

    X2m = im + (2*a2*a5 -a4*a6)/(a6^2-4*a3*a5);
    Y2m = jm + (2*a3*a4 -a2*a6)/(a6^2-4*a3*a5);
    
    return [X2m, Y2m];
}


func CoGApp(image)
/* DOCUMENT func CoGApp(image)
   Sub pixel accuracy measurement with Pyramid approximation
   Sub-pixel estimation of local estrema, DG Bailey, Image and Vision Computing NZ
   
   Eq No. 8
*/
{
    image=image/max(image);
    IndicesMax = where2(image == max(image));
    im =  IndicesMax(1);
    jm =  IndicesMax(2);

    if (image(im-1,jm) < image(im+1,jm))
        {
            minIndicX = image(im-1,jm);
        }else minIndicX = image(im+1,jm);

    if (image(im,jm-1) < image(im,jm+1))
        {
            minIndicY = image(im,jm-1);
        }else minIndicY = image(im,jm+1);

    
    P=minn( image, [im-1, jm], [im+1, jm]);
    Xm = im +(image(im+1,jm)- image(im-1,jm))/(image(im,jm) +image(im+1,jm) + image(im-1,jm) - 3*image(P(1), P(2) ));   

    P=minn( image, [im, jm-1], [im, jm+1]);
    Ym = jm + ( image(im,jm+1)- image(im,jm-1) )/(image(im,jm)+image(im,jm+1) + image(im,jm-1) -3*image(P(1), P(2)) );
    
    return [Xm, Ym];
}

func CoGApp0(image)
/* DOCUMENT func CoGApp(image)
   Sub pixel accuracy measurement with Pyramid approximation
   Sub-pixel estimation of local estrema, DG Bailey, Image and Vision Computing NZ
   
   Eq No. 8
*/
{
    image=image/max(image);
    IndicesMax = where2(image == max(image));
    im =  IndicesMax(1);
    jm =  IndicesMax(2);

    size=dimsof(image)(2);

    if(im==1)
        {
            im1= size; //im-1
            Xm = im + 0.5*(image(im1,jm)- image(im+1,jm))/( image(im1,jm) + image(im+1,jm) -2*image(im,jm));
        }

    if (jm==1)
        {
            jm1=size; //jm-1
            Ym = jm + 0.5*( image(im,jm1)- image(im,jm+1) )/( image(im,jm1) + image(im,jm+1) -2*image(im,jm));
        }

    if (jm==size)
        {
            jm16=1;  //jm+1
            Ym = jm + 0.5*( image(im,jm-1)- image(im,jm16) )/( image(im,jm-1) + image(im,jm16) -2*image(im,jm));
        }

    if (im==size)
        {
            im16=1; //jm+1
            Xm = im + 0.5*(image(im-1,jm)- image(im16,jm))/( image(im-1,jm) + image(im16,jm) -2*image(im,jm));
        }


    if (im !=1 && im !=size)
        {
            Xm = (     (image(im-1,jm)*(im-1) + image(im,jm)*(im) + image(im+1,jm)*(im+1) )/( image(im-1,jm) +image(im,jm)+image(im+1,jm) )                 );
           
        }

    if (jm!=1 && jm !=size)
        {
            Ym = (     (image(im,jm-1)*(jm-1) + image(im,jm)*(jm) + image(im,jm+1)*(jm+1) )/( image(im,jm-1) +image(im,jm)+image(im,jm-1) ));
                             
           
        }

    return [Xm, Ym];
}




func TCoG_multi_app(im, T, window)
/* DOCUMENT TCoGXY_multi_app(im, T, window)

   returns centriods [x,y]

   T ranges from 0.2-1.0;
   window can be adjusted; Image selected from peak-window:peak+window.
 */
{
    extern TCoGImage;
    if(window ==[]){
      window = 2;
    }

    if(T == []){
      T=0.2;
    }
    window=int(window);

    
    Ixy=where2(im==max(im));
    image=im(Ixy(1)-window/2:Ixy(1)+window/2, Ixy(2)-window/2:Ixy(2)+window/2);
    
    size = dimsof(image)(2);
    sumxI=sumyI=sumI1=sumI2=0.0;

    I_T= T*max(image);
    
    for(i=1; i<=size*size; i++)
        {
            if(image(i) < I_T) image(i)=0.0;
        }

    i=[];
    for(i=1; i<=size; i++)
        {
            sumxI += sum(i*image(i,));
            sumyI += sum(i*image(,i));
            sumI1 += sum(image(i,));
            sumI2 += sum(image(,i));
        }
    TCoGImage=image;
    Xc= Ixy(1)-(window/2+1)+sumxI/sumI1;
    Yc= Ixy(2)-(window/2+1)+sumyI/sumI2;
  
    return [Xc,Yc] ;
}


func GaussApp(image)
/* DOCUMENT GaussApp(image)

 */
{
    image=image/max(image);
    z=where(image <= 0);
    image(z)=1e-2;
    
    IndicesMax = where2(image == max(image));
    im =  IndicesMax(1);
    jm =  IndicesMax(2); 

    size=dimsof(image)(2);
    
    if(im==1)
        {
            im1= size; //im-1
            Xm = im + 0.5*(image(im1,jm)- image(im+1,jm))/( image(im1,jm) + image(im+1,jm) -2*image(im,jm));
        }

    if (jm==1)
        {
            jm1=size; //jm-1
            Ym = jm + 0.5*( image(im,jm1)- image(im,jm+1) )/( image(im,jm1) + image(im,jm+1) -2*image(im,jm));
        }

    if (jm==size)
        {
            jm16=1;  //jm+1
            Ym = jm + 0.5*( image(im,jm-1)- image(im,jm16) )/( image(im,jm-1) + image(im,jm16) -2*image(im,jm));
        }

    if (im==size)
        {
            im16=1; //jm+1
            Xm = im + 0.5*(image(im-1,jm)- image(im16,jm))/( image(im-1,jm) + image(im16,jm) -2*image(im,jm));
        }



     if (im !=1 && im !=size)
        {
            Xm = im + 0.5*( log(image(im-1,jm))- log(image(im+1,jm)) )/( log(image(im-1,jm)) + log(image(im+1,jm)) -2*log(image(im,jm)) );
           
        }

     if (jm!=1 && jm !=size)
        {
            Ym = jm + 0.5*( log(image(im,jm-1))- log(image(im,jm+1)) )/( log(image(im,jm-1)) + log(image(im,jm+1)) -2*log(image(im,jm)) );
        }
    
    return [Xm, Ym];

};

func minn(image, P1, P2)
{
   
    return P= (image(P1(1), P1(2)) <= image (P2(1), P2(2))) ? P1:P2;     
}



func maxx(image, P1, P2)
{
   
    return P= (image(P1(1), P1(2)) >= image (P2(1), P2(2))) ? P1:P2;     
}


func PyramidApp(image)
/* DOCUMENT func Pyramid(image)
   Sub pixel accuracy measurement with Pyramid approximation
   Sub-pixel estimation of local estrema, DG Bailey, Image and Vision Computing NZ
   
   Eq No. 4
*/
{
    image=image/max(image);
    IndicesMax = where2(image == max(image(3:13, 3:13)));
    im =  IndicesMax(1);
    jm =  IndicesMax(2); 

    size=dimsof(image)(2);

    
    if(im==1)
        {
            im1= size; //im-1
            Xm = im + 0.5*(image(im1,jm)- image(im+1,jm))/( minn( image, [im1, jm], [im+1, jm]) - image(im, jm) )  ;
        }

    if (jm==1)
        {
            jm1=size; //jm-1
            Ym = jm + 0.5*( image(im,jm1)- image(im,jm+1) )/( minn( image, [im, jm1], [im, jm+1]) - image(im, jm) );
        }

    if (jm==size)
        {
            jm16=1;  //jm+1
            Ym = jm + 0.5*( image(im,jm-1)- image(im,jm16) )/( image(im,jm-1) + image(im,jm16) -2*image(im,jm));
        }

    if (im==size)
        {
            im16=1; //jm+1
            Xm = im + 0.5*(image(im-1,jm)- image(im16,jm))/( image(im-1,jm) + image(im16,jm) -2*image(im,jm));
        }

    if (im !=1 && im !=size)
        {
            P=minn( image, [im-1, jm], [im+1, jm]);
            Xm = im + 0.5*(image(im-1,jm)- image(im+1,jm))/( image(P(1), P(2)) - image(im, jm) );
           
        }

    if (jm!=1 && jm !=size)
        {
            P=minn( image, [im, jm-1], [im, jm+1]);
            Ym = jm + 0.5*( image(im, jm-1)- image(im, jm+1) )/( image(P(1), P(2)) - image(im, jm) );
        }
    
    
    return [Xm, Ym];
}


func Equiangular_line(image)
/* DOCUMENT Equiangular_line(image)
   Shimizu_Masao_2005
 */
{
    
    IndicesMax = where2(image == max(image));
    im =  IndicesMax(1);
    jm =  IndicesMax(2); 

    size=dimsof(image)(2);
    
    if(im==1)
        {
            im1= size; //im-1
            Xm = im + 0.5*(image(im1,jm)- image(im+1,jm))/( image(im1,jm) + image(im+1,jm) -2*image(im,jm));
        }

    if (jm==1)
        {
            jm1=size; //jm-1
            Ym = jm + 0.5*( image(im,jm1)- image(im,jm+1) )/( image(im,jm1) + image(im,jm+1) -2*image(im,jm));
        }

    if (jm==size)
        {
            jm16=1;  //jm+1
            Ym = jm + 0.5*( image(im,jm-1)- image(im,jm16) )/( image(im,jm-1) + image(im,jm16) -2*image(im,jm));
        }

    if (im==size)
        {
            im16=1; //jm+1
            Xm = im + 0.5*(image(im-1,jm)- image(im16,jm))/( image(im-1,jm) + image(im16,jm) -2*image(im,jm));
        }



    if( image(im+1, jm) <= image(im-1, jm) ){

        Xm = im + 0.5*(  (image(im+1, jm)-image(im-1, jm))  / (image(im, jm)-image(im-1, jm)) );
    }else {
        Xm = im + 0.5*(  (image(im+1, jm)-image(im-1, jm))  / (image(im, jm)-image(im+1, jm)) );
    }



    if( image(im, jm+1) <= image(im, jm-1) ){

        Ym = jm + 0.5*(  (image(im, jm+1)-image(im, jm-1))  / (image(im, jm)-image(im, jm-1)) );
    }else {
        Ym = jm + 0.5*(  (image(im, jm+1)-image(im, jm-1))  / (image(im, jm)-image(im, jm+1)) );
    }
    
    return [Xm, Ym];
}






func WCoGXY_Gauss(image, Gauss_sigma)
/* DOCUMENT WCoGXY_Gauss(image, Gauss_sigma)

   Centroding with Gauss weighing function.

   Gauss_sigma [sigma_x, sigma_y] can be varied.
   
   returns centriods [x,y]

 */
{
    size = dimsof(image)(2);
    sumxI=sumyI=sumI1=sumI2=0.0;
    
    if (Gauss_sigma ==[]) Gauss_sigma = [2.0/2.35, 2.0/2.35];
    C=where2(image==max(image));
    
    //C=CoGApp(image);
    ax = [1.0, C(1), Gauss_sigma(1)];
    ay = [1.0, C(2), Gauss_sigma(2)];

     W = span(1,size,size);
     
    //Gaussin weights
    Gx = gauss(W,ax);
    Gy = gauss(W,ay);
    
    for(i=1; i<=size; i++)
        {
          //sum(x*I), sum(y*I)
            sumxI += sum(image(i,)*i*Gx(i));
            sumyI += sum(image(,i)*i*Gy(i));

            //sum(I)s
            sumI1 +=sum(image(i,)*Gx(i));
            sumI2 +=sum(image(,i)*Gy(i));
        }
   
    return [double(sumxI/sumI1), double(sumyI/sumI2)];
}





func SCoG(image)
/* DOCUMENT CoG(image)

   x=sum(x(i)*I(i))/sum(I(i));
   y=sum(y(i)*I(i))/sum(I(i));
   
   returns centriods [x,y]

 */
{
  //image=image(8-3:8+3, 8-3:8+3);
    size = dimsof(image)(2);
    sum_xI=sum_yI=sumI1=sumI2=0.0;
    
    for(i=1; i<=size; i++)
        {
          //sum(x*I), sum(y*I)
          sum_xI += sum(image(i,)*i);
          sum_yI += sum(image(,i)*i);
          
            //sum(I)s
          sumI1 +=sum(image(i,));
          sumI2 +=sum(image(,i));
        }
   
    return [double(sum_xI/sum(image)), double(sum_yI/sum(image))];
}


