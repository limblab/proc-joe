% GETTUNINGCURVES Gets tuning curves for firing rate against
% move_var. Outputs curve in 8 bins, along with high and low CI.
% INPUTS - 
%   trial_data - trial_data struct on which to operate
%   params - parameters struct
%       .out_signals - signal to get tuning curves for (usually firing rates)
%       .out_signal_names : names of signals to be used as signalID pdTable
%                           default - empty
%       .use_trials    : trials to use.
%                         DEFAULT: 1:length(trial_data
%       .move_corr - movement correlate to find tuning to
%                    options:
%                           'vel' : velocity of handle (default)
%                           'acc' : acceleration of handle
%                           'force'  : force on handle
%       .num_bins - number of directional bins (default: 8)
%       .prefix : prefix to add onto columns of weight table
%       .meta   : meta parameters for makeNeuronTableStarter
% OUTPUTS -
%   curves - table of tuning curves for each column in signal, with 95% CI
%
% Written by Raeed Chowdhury. Updated Jul 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [curves] = getTuningCurves(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
out_signals      =  [];
use_trials        =  1:length(trial_data);
move_corr      =  'vel';
num_bins        = 8;
prefix = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters?
calc_CIs = true;
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
possible_corrs = {'vel','acc','force','dlc_vel_handxy','dlc_acc_handxy'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
trial_data = trial_data(use_trials);
if isempty(out_signals), error('Need to provide output signal'); end
if isempty(move_corr), error('Must provide movement correlate.'); end
if ~any(ismember(move_corr,possible_corrs)), error('Correlate not recognized.'); end
out_signals = check_signals(trial_data(1),out_signals);
response_var = get_vars(trial_data,out_signals);
move_corr = check_signals(trial_data(1),move_corr);
move_var = get_vars(trial_data,move_corr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get bins
bins = linspace(-pi,pi,num_bins+1);
bins = bins(2:end);
bin_spacing = mean(unique(diff(bins)));
assert(numel(bin_spacing)==1,'Something went wrong...')

% bin directions
dir = atan2(move_var(:,2),move_var(:,1));
dir_bins = round(dir/bin_spacing)*bin_spacing;
dir_bins(dir_bins==-pi) = pi;

% find response_var in each bin, along with CI
curve = zeros(size(response_var,2),num_bins);
curve_CIhigh = zeros(size(response_var,2),num_bins);
curve_CIlow = zeros(size(response_var,2),num_bins);
for i = 1:num_bins
    % get response_var when move_var is in the direction of bin
    % Also transpose response_var so that rows are neurons and columns are observations
    response_var_in_bin = response_var(dir_bins==bins(i),:)';

    % Mean binned response_var has normal-looking distribution (checked with
    % bootstrapping on a couple S1 neurons)
    curve(:,i) = mean(response_var_in_bin,2); % mean firing rate

    if calc_CIs
        curve_stderr = std(response_var_in_bin,0,2)/sqrt(size(response_var_in_bin,2)); % standard error
        tscore = tinv(0.975,size(response_var_in_bin,2)-1); % t-score for 95% CI
        curve_CIhigh(:,i) = curve(:,i)+tscore*curve_stderr; %high CI
        curve_CIlow(:,i) = curve(:,i)-tscore*curve_stderr; %low CI
    end
end

% set up output struct
if ~isempty(prefix)
    prefix = strcat(prefix,'_');
end
if calc_CIs
    var_names = [{'bins'} strcat(prefix,move_corr{:,1},{'Curve','CurveCIlow','CurveCIhigh'})];
    tab_append = table(repmat(bins,size(response_var,2),1),curve,curve_CIlow,curve_CIhigh,...
                'VariableNames',var_names);
else
    var_names = [{'bins'} strcat(prefix,move_corr{:,1},{'Curve'})];
    tab_append = table(repmat(bins,size(response_var,2),1),curve,...
                'VariableNames',var_names);
    tab_append.Properties.VariableDescriptions = {'circular','linear'};
end
curves = horzcat(makeNeuronTableStarter(trial_data,params),tab_append);

end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
