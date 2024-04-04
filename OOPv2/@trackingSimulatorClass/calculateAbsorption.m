function absorbedInShell = calculateAbsorption(obj, perElementPerEnergyAbsorbed)
    % Get the size of the input array
    [n, m] = size(perElementPerEnergyAbsorbed); 
    
    absorbedInShell = cell(n, m);
    for i = 1 : n
        for j = 1 : m 
            if ~isempty(perElementPerEnergyAbsorbed{i,j}) && perElementPerEnergyAbsorbed{i,j}(1,1)~=0
                value = perElementPerEnergyAbsorbed{i,j}(1,1); 
                idx = findCellIndexWithChar(obj,obj.materialProperties.IXASdata{1,i}(3:end,2), ' Energy (eV)');
                Kedges = str2double(obj.materialProperties.IXASdata{1,i}(3:idx,2));
                absorbedInShell{i,j} = zeros(size(perElementPerEnergyAbsorbed{i,j}));
            
                % Find closer K-edge where energy value is still higher 
                [~, higherIndex] = findIndicesAroundValue(obj, Kedges(1:end), value);
                jumpFactor = str2double(obj.materialProperties.IXASdata{1,i}(higherIndex+2,5));
                absorptionProbability = (jumpFactor-1)/jumpFactor; 
                nAbsorbed = round(absorptionProbability.*length(perElementPerEnergyAbsorbed{i,j}(1 : end))); 
                absorbedInShell{i,j}(1:nAbsorbed) = higherIndex; 
                nAbsorbedPrev = 0; 
                for k = 1 : length(Kedges(higherIndex:end))
                    jumpFactor = str2double(obj.materialProperties.IXASdata{1,i}(k+higherIndex+2,5));
                    absorptionProbability = (jumpFactor-1)/jumpFactor; 
                    nAbsorbedPrev = nAbsorbed+nAbsorbedPrev;
                    nAbsorbed = round(absorptionProbability.*length(perElementPerEnergyAbsorbed{i,j}(nAbsorbedPrev+1 : end))); 
                    if nAbsorbedPrev + nAbsorbed <= length(absorbedInShell{i,j})
                        absorbedInShell{i,j}(nAbsorbedPrev+ 1:nAbsorbedPrev + nAbsorbed) = higherIndex+k;   
                    end  
                end 

            end 
            for kk = 1 : length(absorbedInShell{i,j})
                if absorbedInShell{i,j}(kk) == 0
                    absorbedInShell{i,j}(kk) = length(Kedges);
                end 
            end           
        end 
   end
end 
 