function [trackLengths, perElementPerEnergyTracks] = sampleTrackLength(obj)
    %DESCRIPTION: 
    
    % log-log Interpolate XCOM data mass attenuation coefficient for each
    % spectrum energy and each element of the detector material 

    
    for i = 1 : length(obj.materialProperties.elementList)
        
        for j = 1 : length(obj.materialProperties.elementList)
            xdata = obj.materialProperties.XCOMdata(j).energy;
            xdata = modifyRepeatingElements(obj, xdata, 5); 
            ydata = obj.materialProperties.XCOMdata(j).photoelectric;
            queryEnergy = unique(obj.energyEmissionMap{end}.energy);
            queryEnergy = queryEnergy(queryEnergy~=0); 
            interpValuePar(j,:) = exp(interp1(log(xdata),log(ydata),log(queryEnergy.*(10^-6)),"linear"));
           
        end 
         %weighted sum of all elements cross-sections
         weigth = transpose(str2double(obj.materialProperties.abundanceList)./sum(str2double(obj.materialProperties.abundanceList), 'all')); 
         interpValue(i,:) = sum(times(interpValuePar, weigth),1); 
    end 
    
    %depending of relavtive abundace of each element and the photon energy
    %of the spectrum, sample the track length from a poisson distribution
    %with the interpolated mass attenuation coefficient

    % (str2double(obj.materialProperties.abundanceList)); 

    %For each photon energy upload the per element probability based o mass
    %attenuation coefficient data at its energy 
    
    perElementProbability = zeros(length(queryEnergy), length((str2double(obj.materialProperties.abundanceList))));
    
    fractionByWeigth = (str2double(obj.materialProperties.abundanceList)); 

    perEnergyPhoton = histc(obj.energyEmissionMap{end}.energy, queryEnergy);

    % Normalize 
    perEnergyPhoton = perEnergyPhoton./sum(perEnergyPhoton, 'all');

    for i = 1 : length(queryEnergy)
        for j = 1 : length((str2double(obj.materialProperties.abundanceList)))
            perElementProbability(i,j) = fractionByWeigth(1,j)*((interpValuePar(j,i))./interpValue(j,i)).*perEnergyPhoton(1,i); 
        end 
    end 

    totalPhotons = length(obj.energyEmissionMap{end}.energy);
    
    perElementperEnergyPhoton = round(perElementProbability.*totalPhotons); 


    if sum(perElementperEnergyPhoton, 'all') ~= totalPhotons 
        diff = totalPhotons - sum(perElementperEnergyPhoton, 'all');
        rand_idx = randi(length(queryEnergy)); 
        rand_idy = randi(length((str2double(obj.materialProperties.abundanceList)))); 
        perElementperEnergyPhoton(rand_idx,rand_idy) = abs(perElementperEnergyPhoton(rand_idx,rand_idy) + diff); 
    end 
 
    % Initialize matrix to store combinations count
    perElementperEnergyPhoton = perElementperEnergyPhoton.';

    trackLengths = arrayfun(@(massAttCoeff,count) computeTracks(obj,massAttCoeff,count), interpValuePar, perElementperEnergyPhoton, 'UniformOutput', false); 
    

    
    perElementPerEnergyTracks = trackLengths; 
    
    trackLengths = ConcatenateCellArrayColumnRow(obj,trackLengths); 

    
 

end
