function [ p ] = pearsonChiSquareBinomial( p1,p2,n1,n2,tail )
% implements unpooled pearson chi-square test (expresed as a z-stat)
% p1 and n1 are the binomial prob and number of trials for sample 1
% p2 and n2 are the binomial prob and number of trials for sample 2

% tail determines tail for normal distribution. upper is p1 > p2;
% lower if p1 < p2
% both if p1 ~= p2

    sigmaD = sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2);

    z = (p1-p2)/sigmaD;

    p = normcdf(z);
    
    switch tail
        case 'upper'
            p=1-p;
        case 'lower'
            p=p;
        case 'both'
            if(p1 > p2)
                p = 2*(1-p);
            else
                p=2*p;
            end
        otherwise
            error('tail is inputted improperly');
    end

end

