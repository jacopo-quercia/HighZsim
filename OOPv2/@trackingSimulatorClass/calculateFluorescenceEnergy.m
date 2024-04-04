function emittedEnergy = calculateFluorescenceEnergy(obj, absorbedInShell, perElementPerEnergyAbsorbed)
    %DESCRIPTION: 
    %Compute energy of K-shell fluorescence photons 
    
    [n, m] = size(absorbedInShell); 
    emittedEnergy = cell(n,m);
    for i = 1 : n
        for j = 1 : m
            emittedEnergy{i,j} = zeros(size(absorbedInShell{i,j}));
            shellValues = unique(absorbedInShell{i,j}); 
            shellValues = shellValues(shellValues~=0); 
            if ~isempty(shellValues)
                for k = 1 : length(shellValues)
                    idx = find(absorbedInShell{i,j} == shellValues(k));
                    shellString = "("+obj.materialProperties.IXASdata{1,i}(2+shellValues(k),1);
                    startidx = find(contains(obj.materialProperties.IXASdata{1,i}(:,1),'Line'));
                    strindx = find(contains(obj.materialProperties.IXASdata{1,i}(startidx+1:end,1),shellString)); 
                    strindx = strindx + startidx; 
                    intensity = str2double(obj.materialProperties.IXASdata{1,i}(strindx,3));
                    intensity = intensity/sum(intensity,'all'); 
                    fluorescenceYield = str2double(obj.materialProperties.IXASdata{1,i}(2+shellValues(k),4)); 
                    %Find indices in absorbedInShell related to the shell 
                    values = str2double(obj.materialProperties.IXASdata{1,i}(strindx,2)); 
                    counts = round(intensity.*fluorescenceYield.*length(idx));
                    % Initialize an empty array to store the results
                    resultArray = zeros(1,length(idx));
                    concArray = [];
                    % Loop through each value and its corresponding count
                    for l = 1:length(values)
                    concArray = [concArray,repmat(values(l), 1, counts(l))];
                    end
                    resultArray(1:sum(counts,'all')) = concArray; 
                    emittedEnergy{i,j}(idx) = resultArray;
                end 
            end
 
        end 
    end 
    

end 