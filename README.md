Convetional method: func IM_center
Method1: func Method1
Method2: func SPcenter



Yorick files:

1) Libraries

libCorrelationAlgorithmsGenerationV4.i
libCorrelation_algorithmV9.i
libGauss3.i
libSubPixelAlgV1.i
libImages.i


2) Scripts that create fits files

//scripts to generate fits
Test_Fig1_SystematicErros.i
Test_Fig3_SubPixel_Fitting_all_2017Jan29.i
Test_Fig4_7_SNR_vs_bias_vs_sigma_for_all_objects.i
Test_Fig8_TwoStep_vs_Traditional_all_2017Feb21.i
Test_Fig9_Two_step_SNR.i     
Test_Fig10_Systematic_errors_with_increasing_interpolation.i


functions generates fits files to plot later:
write_subpixel_data(object)
write_SNR_data(object)
write_bias_error_vs_interpolation(void)
write_two_step_vs_traditional(object)
write_bias_error_vs_interpolation(void)




3) Scripts that plot figures from fits file

PlotAll.i
Fig1_SystematicErrorPlot.i                         
Fig2_SubapertureImagesV3.i                         
Fig3_SubPixel_fitting_alg_performance_20160921.i   
Fig4_SNR_Point1.i                                  
Fig5_SNR_LGS1.i                                   
Fig6_SNR_GC1.i                                     
Fig7_SNR_Solar1.i                                  
Fig8_TwoMethods_vs_Traditional.i  
Fig9_Two_step_SNRNew.i
Fig10_bias_vs_interpolation_performanceV4_plotV2.i  

4) Computational cost test: 

ComputationalTime.i   

5) 5e3 e-/image to star near-infrared star magnitude

FluxSH.i                             

5) fits files

For solare image: gband_22May2002_AR9957_1.fits

Data_bias_vs_interpolation.fits             Data_Two_step_traditional_GC.fits     GC_SNR.fits
Data_performance_subpixel_GCV2_5_5.fits     Data_Two_step_traditional_LGS.fits    LGS_SNR.fits
Data_performance_subpixel_LGSV2_5_5.fits    Data_Two_step_traditional_point.fits  point_SNR.fits
Data_performance_subpixel_pointV2_5_5.fits  Data_Two_step_traditional_solar.fits  solar_SNR.fits
Data_performance_subpixel_solarV2_5_5.fits  





