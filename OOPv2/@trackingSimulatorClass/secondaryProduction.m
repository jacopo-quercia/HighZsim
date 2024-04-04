function obj = secondaryProduction(obj, perElementPerEnergyTracks)
%DESCRIPTION:
%

%Divide all absorbed energy in subsets based on material and energy
[n, m] = size(perElementPerEnergyTracks);

% Initialize the lengths array with the same dimensions
lengths = zeros(n, m);

% Iterate over each cell in the cell array
for i = 1:n
    for j = 1:m
        % Get the length of the array in the current cell
        % and store it in the corresponding position in the lengths array
        lengths(i, j) = length(perElementPerEnergyTracks{i, j});
    end
end

lengths = adjustArrayToTargetSum(obj, lengths, length(obj.energyDepositionMap{end}.x));

perElementPerEnergyAbsorbed = reverseConcat(obj, lengths, n, m);

obj = computeFluorescence(obj, perElementPerEnergyAbsorbed);



end

