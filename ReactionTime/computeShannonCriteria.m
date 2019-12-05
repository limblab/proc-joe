function [ k ] = computeShannonCriteria( I,T,area )
    % I in A, T in us, area in cm^2

    k = log10(I*T/area) + log10(I*T);
end

