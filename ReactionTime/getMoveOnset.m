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
peak_idx_offset = [0,100];
start_idx_offset = 0;
which_field = 'speed';
field_idx = 1;
threshold_mult = 0.5;
pre_move_thresh = 0.6;
max_rt = 0.5;
% these parameters aren't documented because I expect them to not need to
% change but you can overwrite them if you need to.
start_idx     =  'idx_goCueTime';
end_idx       =  'idx_endTime';
onset_name    =  'movement_on';
if nargin > 1, assignParams(who,params); end % overwrite defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some pre-processing
td = getSpeed(trial_data);

for trial = 1:length(trial_data)
    % use which_field to find bin corresponding to movement onset, movement offset, and peak speed
    s = td(trial).(which_field)(:,field_idx);
    
    % find the time bins where the monkey may be moving
    move_inds = false(size(s));
    move_inds(td(trial).(start_idx)+start_idx_offset:td(trial).(end_idx)) = true;
    
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
        mvt_peaks = find(peaks & (1:length(peaks))' > td(trial).(start_idx)+peak_idx_offset(1) & (1:length(peaks))' < td(trial).(start_idx)+peak_idx_offset(2) & s > min_s & move_inds);
        [~,mvt_peak] = max(s(mvt_peaks));
        mvt_peak = mvt_peaks(mvt_peak);
        if ~isempty(mvt_peak)
            thresh = s(mvt_peak)*threshold_mult; % default threshold is half max acceleration peak
            on_idx = find(s<thresh & (1:length(s))'<mvt_peak & move_inds,1,'last');
            
            % check to make sure the numbers make sense
            if on_idx <= td(trial).(start_idx)
                % something is fishy. Fall back on threshold method
                on_idx = NaN;
            elseif(on_idx > td(trial).(start_idx) + max_rt/mode([td.bin_size]))
                on_idx = NaN;
            end
        end
        % peak is max velocity during movement
        temp = s; temp(~move_inds) = 0;
        [~,peak_idx] = max(temp);
    end
    
    if isempty(on_idx) || isnan(on_idx)
%         on_idx = find(s > s_thresh & move_inds,1,'first');
        if isempty(on_idx) % usually means it never crosses threshold
            warning('Could not identify movement onset');
            on_idx = NaN;
        end
    end
    trial_data(trial).(['idx_' onset_name]) = on_idx;

end

% restore logical order
trial_data = reorderTDfields(trial_data);

end



