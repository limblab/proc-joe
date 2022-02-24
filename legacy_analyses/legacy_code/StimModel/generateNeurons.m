function [neurons] = generateNeurons(numNeurons, minDistance, maxDistance)
% generates a list of neurons by randomly assigning preferred directions and
% distances. neurons are represented as (distance, x_dir, y_dir) where
% x_dir and y_dir have a magnitude of 1 (Norm2(x_dir,y_dir)^2 = 1)

directions = rand(numNeurons,1)*2*pi; % random numbers [0, 2pi)
neurons = [rand(numNeurons,1)*(maxDistance-minDistance)+minDistance, ... % random distance
    cos(directions),sin(directions)]; % directions as x and y components

end