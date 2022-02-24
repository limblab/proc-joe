function [ p ] = poissTestHuffman( lambda1, lambda2, s1, s2 )
% implements poisson two sample rate test in Huffman 1984
% lambda1; occurances observed for pop 1 (rate1*s1)
% lambda2, occurences observed for pop 2 (rate2*s2)
% s1; observation time length
% s2; observation time length 

    

    d = s2/s1;
    z = 2*(sqrt(lambda2 + 3/8) - sqrt(d*(lambda1 + 3/8)))/sqrt(1+d);

    p = normcdf(z,'upper');
end

