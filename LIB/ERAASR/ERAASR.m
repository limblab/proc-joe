function [ Xcleaned ] = ERAASR( Xinput, opts )
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
    STEP_NAMES = {'Channels','Pulses','Trials'};
    Xraw = Xinput;
    Xcleaned.input = Xinput;
    %% temporally align stimulation trials
    if(opts.ALIGN_TRIALS)
        Xcleaned.inputAligned = alignTrials_ERAASR(Xcleaned.input,opts);
        Xraw = Xcleaned.inputAligned;
    end
    %% clean across channels, then pulses, then trials
    for stepIdx = [1,2,3]  
        if(opts.K(stepIdx) > 0)
            % reshape array
            Mc = reshape_ERAASR(Xraw,STEP_NAMES{stepIdx});
            McProj = Mc;
            % get PCA weights
            Wc = getPCAProjection(Mc,opts.K(stepIdx));
            % remove projection for each channel
            X = zeros(size(Xraw));
            for c = 1:size(Mc,2)
                WcZeroed = Wc;
                WcZeroed(max(1,c-opts.LAMBDA(stepIdx)):min(size(Mc,2),c+opts.LAMBDA(stepIdx)),:) = 0;
                Ac = Mc*WcZeroed;
                WcProj = Ac\Mc(:,c);
                McProj(:,c) = Mc(:,c) - Ac*WcProj;
            end
            Xraw = unshape_ERAASR(McProj,size(Xinput),STEP_NAMES{stepIdx});
            Xcleaned.(strcat('post',STEP_NAMES{stepIdx})) = Xraw;
        end
    end
end

function [opts] = configureOpts(optsInput)

    opts.K = [4,2,1];
    opts.LAMBDA = [0,0,0];
    opts.BETA = [0,1]; % true or false
    opts.ALIGN_TRIALS = 0; % not implemented currently
    opts.UPSAMPLE_FACTOR = 10;
    opts.MAX_LAG = 10;
    
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
        case 'Channels'
            temp = permute(X,[2,1,3,4]);
            out = reshape(temp,[size(temp,1)*size(temp,3)*size(temp,4),size(temp,2)]);
        case 'Pulses'
            temp = permute(X,[2,3,1,4]);
            out = reshape(temp,[size(temp,1)*size(temp,3)*size(temp,4),size(temp,2)]);
        case 'Trials'
            temp = permute(X,[2,4,3,1]);
            out = reshape(temp,[size(temp,1)*size(temp,3)*size(temp,4),size(temp,2)]);
    end
end

function [out] = unshape_ERAASR(X,sizes,step)
    out = [];    
    switch step
        case 'Channels'
            temp = reshape(X,[sizes(2),sizes(1),sizes(3),sizes(4)]);
            out = permute(temp,[2,1,3,4]);
        case 'Pulses'
            temp = reshape(X,[sizes(2),sizes(3),sizes(1),sizes(4)]);
            out = permute(temp,[3,1,2,4]);
        case 'Trials'
            temp = reshape(X,[sizes(2),sizes(4),sizes(3),sizes(1)]);
            out = permute(temp,[4,1,3,2]);
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


function [Xalign] = alignTrials_ERAASR(Xinput,opts)
    % this function aligns trials by upsampling, computing an
    % crosscorrelation, shifting, then downsampling.
    
    % upsample
    Xalign = Xinput;
    Xupsample = zeros(size(Xinput,1),size(Xinput,2)*opts.UPSAMPLE_FACTOR,size(Xinput,3),size(Xinput,4));
    for channel = 1:size(Xinput,1)
        for pulse = 1:size(Xinput,3)
            for trial = 1:size(Xinput,4)
               Xupsample(channel,:,pulse,trial) = interp(squeeze(Xinput(channel,:,pulse,trial)),opts.UPSAMPLE_FACTOR); 
               % compute an crosscorrelation comparing each trial to the first one
                if(trial > 1)
                    % shift
                    [~,shift] = max(xcorr(squeeze(Xupsample(channel,:,pulse,1)),squeeze(Xupsample(channel,:,pulse,trial)),opts.MAX_LAG));
                    shift = (shift-opts.MAX_LAG); 
                    Xupsample(channel,:,pulse,trial) = circshift(squeeze(Xupsample(channel,:,pulse,trial)),shift);
                    % downsample
                    Xalign(channel,:,pulse,trial) = Xupsample(channel,1:opts.UPSAMPLE_FACTOR:end,pulse,trial);
                end
            end
        end
    end

        
end