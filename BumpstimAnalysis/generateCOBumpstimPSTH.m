function [  ] = generateCOBumpstimPSTH(cds,neuronNumber,varargin)
% generates a raster for the COBumpstim experiments -- wrapper for the plot
% function

% separate trials into no bump no stim, only bump, both bump and stim
% also remove abort and incomplete trials
noBumpNoStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 0 & ...
    cds.trials.delayBump == 0 & cds.trials.moveBump == 0 & isnan(cds.trials.stimCode);
ctrHoldBumpNoStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 1 & ...
    cds.trials.delayBump == 0 & cds.trials.moveBump == 0 & isnan(cds.trials.stimCode); 
ctrHoldBumpStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 1 & ...
    cds.trials.delayBump == 0 & cds.trials.moveBump == 0 & ~isnan(cds.trials.stimCode); 
delayBumpNoStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 0 & ...
    cds.trials.delayBump == 1 & cds.trials.moveBump == 0 & isnan(cds.trials.stimCode); 
delayBumpStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 0 & ...
    cds.trials.delayBump == 1 & cds.trials.moveBump == 0 & ~isnan(cds.trials.stimCode); 
moveBumpNoStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 0 & ...
    cds.trials.delayBump == 0 & cds.trials.moveBump == 1 & isnan(cds.trials.stimCode); 
moveBumpStimMask = cds.trials.result == 'R' & ~isnan(cds.trials.ctrHold) & cds.trials.ctrHoldBump == 0 & ...
    cds.trials.delayBump == 0 & cds.trials.moveBump == 1 & ~isnan(cds.trials.stimCode); 

% plot rasters for ctr hold bumps if there are those
if(sum(ctrHoldBumpNoStimMask) ~= 0 && sum(ctrHoldBumpStimMask) ~= 0)
    plotPSTHCtrHoldBump(cds,neuronNumber,ctrHoldBumpNoStimMask,ctrHoldBumpStimMask,varargin{:});
end

% plot rasters for ctr hold bumps if there are those
if(sum(delayBumpNoStimMask) ~= 0 && sum(delayBumpStimMask) ~= 0)
    % currently empty because I don't do these
end

% plot rasters for ctr hold bumps if there are those
if(sum(moveBumpNoStimMask) ~= 0 && sum(moveBumpStimMask) ~= 0)
    % currently empty because I don't do these
end

end

