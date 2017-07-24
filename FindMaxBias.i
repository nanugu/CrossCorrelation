Data1=fits_read("Fits/Data_performance_subpixel_pointV2_5_5.fits");
Data2=fits_read("Fits/Data_performance_subpixel_LGSV2_5_5.fits");
Data3=fits_read("Fits/Data_performance_subpixel_GCV2_5_5.fits");
Data4=fits_read("Fits/Data_performance_subpixel_solarV2_5_5.fits");


func FindMaxBias(Data){
  D=array(0.0, 6);
  for(i=2; i<=7; i++){
    j=where(Data(,i)(10:30) == max(Data(,i)(10:30)));
    D(i-1) = Data(,1)(10:30)(j);
  }
  write, "QI =", D(1), "PA=", D(2), "PY=", D(4), "GA=", D(5), "TCoG=", D(6);
  return D;
}

FindMaxBias(Data1)
FindMaxBias(Data2)
FindMaxBias(Data3)
FindMaxBias(Data4)
