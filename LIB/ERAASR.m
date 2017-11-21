function [ Xcleaned ] = ERAASR( Xraw, opts )
% implements O'Shea 2017 ERAASR method according to the paper
% INPUT: Xraw is the raw data tensor with size (C,T,P,R) with C = channels, T =
% timepoints per pulse, P = pulses, R = trials

% OUTPUT: Xcleaned with same size

% PARAMETERS (opts): K = number of principal components to
% describe artifact structure over channels, pulses trials (3,1) array
%  lambda = number of adjacent channels to remove (3,1) array
% beta = perform cleaning over pulses separately on each channel? true or
% false, (2,1) array

    %% configure options
    opts = configureOpts(opts);
    
    %% clean across channels
    stepIdx = 1;
    % reshape array
    Mc = reshape_ERAASR(Xraw,'channels');
    McProj = Mc;
    % get PCA weights
    Wc = getPCAProjection(Mc,opts.K(stepIdx));
    % remove projection for each channel
    X = zeros(size(Xraw));
    for c = 1:size(Mc,2)
        WcZeroed = Wc;
        WcZeroed(max(1,c-opts.LAMBDA(stepIdx)):min(size(Mc,2),c+opts.LAMBDA(stepIdx)),:) = 0;
        Ac = Mc*WcZeroed;
        McProj(:,c) = Mc(:,c)-inv(Ac'*Ac)*Ac'*Mc(:,1);
    end
    X = unshape_ERAASR(McProj);
    % X is the output from this step
    
    %% clean across pulses
    stepIdx = 2;
    Mp = reshape_ERAASR(X,'pulses');
    
    % X is the output from this step
    
    %% clean across trials
    stepIdx = 3;
    Mr = reshape_ERAASR(X,'trials');
end

function [opts] = configureOpts(optsInput)

    opts.K = [4,2,1];
    opts.LAMBDA = [0,0,0];
    opts.BEAT = [0,1]; % true or false

    %% check if in opts and optsInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end

% this function will do all of our reshaping. 
function [out] = reshape_ERAASR(X,step)
    % time points are in the second column of X, so when we reshape, we
    % need to keep those aligned
    % step is the reshaping step -- 'channels','pulses','trials'
    switch step
        case 'channels'
            temp = permute(X,[2,1,3,4]);
            out = reshape(temp,[size(temp,1)*size(temp,3)*size(temp,4),size(temp,2)]);
        case 'pulses'
            temp = permute(X,[2,3,1,4]);
            out = reshape(temp,[size(temp,1)*size(temp,3)*size(temp,4),size(temp,2)]);
        case 'trials'
            temp = permute(X,[2,4,3,1]);
            out = reshape(temp,[size(temp,1)*size(temp,3)*size(temp,4),size(temp,2)]);
    end
end

function [out] = unshape_ERAASR(X,chan,sizes,step)
    out = [];    
    switch step
        case 'channels'
            temp = reshape(X,[sizes(2),sizes(1),sizes(3),sizes(4)]);
            out = permute(temp,[2,1,3,4]);
        case 'pulses'
            temp = reshape(X,[sizes(2),sizes(3),sizes(1),sizes(4)]);
            out = permute(temp,[3,1,2,4]);
    end
    
end

function [weights] = getPCAProjection(X,K)
    % this function gets the K pca weights for the matrix X
    
    [coeff,score] = pca(X);
    
    if(K > size(score,2))
        weights = coeff(:,:);
    else
        weights = coeff(:,1:K);
    end
end
