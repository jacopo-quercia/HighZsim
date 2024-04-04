function trackLenghts = computeTracks(obj,massAttCoeff,count)
%DESCRIPTION:
%
    if massAttCoeff~=0 && count~=0
        density = obj.materialProperties.density;
        tau = 1/(density.*massAttCoeff./0.01); 
        trackLenghts = exprnd(tau,1, count); 
    else if massAttCoeff==0 && count~=0
        trackLenghts = zeros(1, count);
    else if count==0
        trackLenghts = [];  
    end 

    end 
end 

