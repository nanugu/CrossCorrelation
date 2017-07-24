
mH=9.5;
NsubApertures=9;
IntTime=10e-3;


//lib functions
func MeasureFlux(mH, NsubApertures, IntTime)
/* DOCUMENT MeasureFlux(mH, NsubApertures, IntTime)
   
   mH-> magnitude of a star 
   Zero point H band (1.654um) flux density N0-> 9.56e9/s/m^2/um
   (from Tokunaga table 7.5, page 156/7, allen astrophysical quantities)
   
   dlambda = width of Hband chosen =0.4um
   A = telescope collecting area = 49.29^2 for VLT
   (due to central obstruction  {(D-Dobs)*(D-Dobs)}^2)

   T=80%
   
   Nsub = effective number of subapertures = 
   = Area of subaperture / telescope area 
   = A/ (8m/NsubApertures)^2
 */
    
{
  
    if(is_void(IntTime)) IntTime=1e-3;
    
    T = 0.8; /* Mirror reflections and detector efficiency */
    A = 49.22; // meters^2
    N_EffectiveSubs = A/(8.0/double(NsubApertures))^2; //no units
    N0 = 9.56e9; //  1/s/m^2/um
    dlambda = 0.5; //um Hband window width
    
    return N0*10^(-0.4*mH)*A*dlambda*T/N_EffectiveSubs*IntTime; // e-/time
}


write, "flux (e-/image) for 8-m tel., 9x9 lenslet = ", MeasureFlux(mH, NsubApertures, IntTime)
