%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getMoveOnset(trial_data,min_ds)
%
%   This will find a time bin representing movement onset and peak speed
% Currently for kinematics but could be easily adapted to force, etc.
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : parameter struct
%     .which_method : (string) how to compute
%                           'peak' : uses acceleration and peak speed
%                           'thresh' : uses a basic velocity threshold
%                       Note: defaults to peak. In thresh, will not return
%                       a peak speed as a field.
%     .min_ds       : minimum diff(speed) to find movement onset
%     .s_thresh     : % speed threshold in cm/s (secondary method if first fails)
%     .peak_idx_offset  :  indices after start_idx to find max speed. 
%     .start_idx_offset :  indices after start_idx to find movement onset
%     .which_field     : which field to find movement onset from 
%     .field_idx       : idx of the above field to find movement onset from
%     .threshold_mult  : threshold multiplier
% OUTPUTS:
%   trial_data : same struct, with fields for:
%       idx_peak_speed
%       idx_movement_onset
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getMoveOnset(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
which_method  =  'peak';
min_s        =  0.3;
s_thresh      =  7;
peak_idx_offset = [0,1000];
start_idx_offset = 0;
max_rt_offset = 400;
which_field = 'speed';
field_idx = 1;
threshold_mult = 0.5;
pre_move_thresh = 0.6;
max_rt = 0.35;
be_aggressive = 0;
use_emg = 0;
emg_idx = [];
threshold_acc = -1; % less than 0, not used

% these parameters aren't documented because I expect them to not need to
% change but you can overwrite them if you need to.
start_idx     =  'idx_goCueTime';
end_idx       =  'idx_endTime';
onset_name    =  'movement_on';
if nargin > 1, assignParams(who,params); end % overwrite defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some pre-processing
td = getSpeed(trial_data);


if(be_aggressive) % find a threshold based on all trials
    s_all = [];
    if(threshold_acc <= 0)
        for trial = 1:length(trial_data)
            if(isfield(td(trial),'tgtDir') && strcmpi(which_field,'speed')~=1)
                % project (which_field) onto the target axis
                s_all = [s_all;sum(td(trial).(which_field)(td(trial).idx_tgtOnTime+35:td(trial).(start_idx),:)*[cos(td(trial).tgtDir*pi/180);sin(td(trial).tgtDir*pi/180)],2)];
            else
                s_all = [s_all;td(trial).(which_field)(td(trial).idx_tgtOnTime+35:td(trial).(start_idx),1)];
            end
        end
        s_all = s_all(~isnan(s_all));
        thresh_all = std(s_all);
    end
    
    % if emg exists, do emg as well
    if(use_emg && isfield(td(1),'emg'))
        emg_all = [];
        for trial = 1:length(trial_data)
                emg_all = [emg_all;td(trial).emg(td(trial).idx_tgtOnTime+35:td(trial).(start_idx),emg_idx)];
        end
        emg_all = emg_all(~isnan(emg_all));
        emg_thresh_all = std(emg_all);

    end
end

for trial = 1:length(trial_data)
    if(isfield(td(trial),'tgtDir') && strcmpi(which_field,'speed')~=1)
%         project (which_field) onto the target axis
        s = sum(td(trial).(which_field)*[cos(td(trial).tgtDir/180*pi);sin(td(trial).tgtDir/180*pi)],2);
    else
        s = td(trial).(which_field)(:,1);
%         s = sqrt(sum(td(trial).(which_field).^2,2));
    end
    
    if(use_emg && isfield(td(1),'emg'))
        emg = td(trial).emg(:,emg_idx);
    end
    
    % find the time bins where the monkey may be moving
    move_inds = false(size(s));

    move_inds(td(trial).(start_idx)+start_idx_offset-3:td(trial).(end_idx)) = true;
    
    [on_idx] = deal(NaN);
    
    % detect moving before go cue
    flag_move = 0;
    if(sum(s(max(1,td(trial).(start_idx)-30):td(trial).(start_idx)) > pre_move_thresh) > 2)
        flag_move = 1;
        on_idx = [];
    end
    if ~flag_move && strcmpi(which_method,'peak')
        ds = [0; diff(s)];
        peaks = [ds(1:end-1)>0 & ds(2:end)<0; 0];
        mvt_peaks = find(peaks & (1:length(peaks))' > td(trial).(start_idx)+peak_idx_offset(1) & ...
            (1:length(peaks))' < td(trial).(start_idx)+peak_idx_offset(2) & s > min_s & move_inds);

        if ~isempty(mvt_peaks)
            [~,mvt_peak] = max(s(mvt_peaks));
%             mvt_peak = 1;
            mvt_peak = mvt_peaks(mvt_peak);
            if(threshold_acc > 0)
                thresh = threshold_acc;
            elseif(~be_aggressive)
                thresh = s(mvt_peak)*threshold_mult; % default threshold is half max acceleration peak
            else
                thresh = thresh_all*threshold_mult;
            end
            on_idx = find(s<thresh & (1:length(s))'<mvt_peak & move_inds,1,'last')+1;

%             on_idx = find(s > thresh & (1:length(s))'<mvt_peak & move_inds,1,'first');
            % check to make sure the numbers make sense
            if isempty(on_idx) || on_idx <= td(trial).(start_idx)+start_idx_offset || on_idx > td(trial).(start_idx)+max_rt_offset
                % something is fishy. Fall back on threshold method
                on_idx = NaN;
            end
        end
        if(isempty(mvt_peaks))
           disp('ugh'); 
        end
        % peak is max velocity during movement
        temp = s; temp(~move_inds) = 0;
        [~,peak_idx] = max(temp);
    end
    
    if isempty(on_idx) || isnan(on_idx)
%         on_idx = find(s > s_thresh & move_inds,1,'first');
        if isempty(on_idx) % usually means it never crosses threshold
            on_idx = NaN;
        end
        if(~isnan(td(trial).tgtDir))
            warning('Could not identify movement onset');
        end
    end
    trial_data(trial).(['idx_' onset_name]) = on_idx;

    
    % get on_idx using emg if requested
    if(use_emg && isfield(td(trial),'emg'))
        % detect moving before go cue
        flag_move = 0;
        if(sum(emg(max(1,td(trial).(start_idx)-30):td(trial).(start_idx)) > pre_move_thresh) > 2)
            flag_move = 1;
            emg_on_idx = [];
        end
        if ~flag_move && strcmpi(which_method,'peak')
            demg = [0; diff(emg)];
            peaks = [demg(1:end-1)>0 & demg(2:end)<0; 0];
            mvt_peaks = find(peaks & (1:length(peaks))' > td(trial).(start_idx)+peak_idx_offset(1) & (1:length(peaks))' < td(trial).(start_idx)+peak_idx_offset(2) & move_inds);
            [~,mvt_peak] = max(emg(mvt_peaks));
            mvt_peak = mvt_peaks(mvt_peak);
            if ~isempty(mvt_peak)
                if(~be_aggressive)
                    thresh = emg(mvt_peak)*threshold_mult; % default threshold is half max acceleration peak
                else
                    thresh = emg_thresh_all*threshold_mult;
                end
                emg_on_idx = find(emg>=thresh & (1:length(emg))'<mvt_peak & move_inds,1,'first');

                % check to make sure the numbers make sense
                if isempty(emg_on_idx) || emg_on_idx <= td(trial).(start_idx)+start_idx_offset || emg_on_idx > td(trial).(start_idx)+max_rt_offset
                    % something is fishy. Fall back on threshold method
                    emg_on_idx = NaN;
                end
            end
            % peak is max velocity during movement
            temp = emg; temp(~move_inds) = 0;
            [~,peak_idx] = max(temp);
        end

        if isempty(emg_on_idx) || isnan(emg_on_idx)
    %         on_idx = find(s > s_thresh & move_inds,1,'first');
            if isempty(emg_on_idx) % usually means it never crosses threshold
                emg_on_idx = NaN;
            end
            warning('Could not identify movement onset');

        end
        trial_data(trial).(['idx_' onset_name,'_emg']) = emg_on_idx;
    end
end

% restore logical order
trial_data = reorderTDfields(trial_data);

end



