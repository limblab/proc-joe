%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pdTable = getTDPDs(trial_data,params)
%
%   Gets PD table for given out_signal. You need to define the out_signal
% and move_corr parameters at input.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .out_signals  : which signals to calculate PDs for
%       .out_signal_names : names of signals to be used as signalID pdTable
%                           default - empty
%       .trial_idx    : trials to use.
%                         DEFAULT: 1:length(trial_data
%       .in_signals   : which signals to calculate PDs on
%                           note: each signal must have only two columns for a PD to be calculated
%                           default - 'vel'
%       .distribution : distribution to use. See fitglm for options
%       .bootForTuning : whether to bootstrap for tuning significance and CI
%           (default: true)
%       .num_boots    : # bootstrap iterations to use (default: 1000)
%       .do_plot      : plot of directions for diagnostics, not for general
%                       use.
%       .prefix       : prefix to add before column names (will automatically include '_' afterwards)
%       .meta   : meta parameters for makeNeuronTableStarter
%
% OUTPUTS:
%   pdTable : calculated velocity PD table with CIs
%
% Written by Raeed Chowdhury. Updated Nov 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pdTable = getTDPDs3D(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
out_signals      =  [];
trial_idx        =  1:length(trial_data);
in_signals      = 'vel';
num_boots        =  1000;
distribution = 'Poisson';
bootForTuning = true;
do_plot = false;
prefix = '';
verbose = true;
doing_3D_pd=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
assert(~isempty(out_signals),'Need to provide output signal')

out_signals = check_signals(trial_data(1),out_signals);
response_var = get_vars(trial_data(trial_idx),out_signals);

in_signals = check_signals(trial_data(1),in_signals);
num_in_signals = size(in_signals,1);
for i = 1:num_in_signals
    assert(length(in_signals{i,2})==2 | length(in_signals{i,2})==3,'Each element of in_signals needs to refer to only two or 3 column covariates')
end
input_var = get_vars(trial_data(trial_idx),in_signals);

if ~isempty(prefix)
    if ~endsWith(prefix,'_')
        prefix = [prefix '_'];
    end
end

if num_boots<2
    bootForTuning = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate final table
dirArr = zeros(size(response_var,2),1);
dirCIArr = zeros(size(response_var,2),2);
moddepthArr = zeros(size(response_var,2),1);
moddepthCIArr = zeros(size(response_var,2),2);
isTuned = false(size(response_var,2),1);
tab_append = cell(1,size(in_signals,1));
for in_signal_idx = 1:size(in_signals,1)
    if bootForTuning
        tab_append{in_signal_idx} = table(dirArr,dirCIArr,moddepthArr,moddepthCIArr,isTuned,...
                        'VariableNames',strcat(prefix,in_signals{in_signal_idx,1},{'PD','PDCI','Moddepth','ModdepthCI','Tuned'}));
        tab_append{in_signal_idx}.Properties.VariableDescriptions = {'circular','circular','linear','linear','logical'};
    else
        tab_append{in_signal_idx} = table(dirArr,moddepthArr,...
                        'VariableNames',strcat(prefix,in_signals{in_signal_idx,1},{'PD','Moddepth'}));
        tab_append{in_signal_idx}.Properties.VariableDescriptions = {'circular','linear'};
    end

    if do_plot
        h{in_signal_idx} = figure;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate PD
bootfunc = @(data) transpose(glmfit(data(:,2:end),data(:,1),distribution));
scramblefunc = @(data) transpose(glmfit(input_var,data,distribution));
unit_tic = tic;
for uid = 1:size(response_var,2)
    %bootstrap for firing rates to get output parameters
    data_arr = [response_var(:,uid) input_var];
    % check if actually bootstrapping
    if bootForTuning
        boot_coef = bootstrp(num_boots,@(data) bootfunc(data),data_arr);

        % get scrambled moddepth
        scramble_coef = bootstrp(num_boots,@(response) scramblefunc(response), response_var(:,uid));
    else
        % don't bootstrap
        boot_coef = bootfunc(data_arr);
    end

    for in_signal_idx = 1:size(in_signals,1)
        move_corr = in_signals{in_signal_idx,1};

        dirs = atan2(boot_coef(:,1+in_signal_idx*2),boot_coef(:,in_signal_idx*2));
        tab_append{in_signal_idx}.([prefix move_corr 'PD'])(uid,:)=circ_mean(dirs);
        if(doing_3D_pd) % compute z-angle as well
            z_ang = acos(boot_coef(:,2+in_signal_idx*2)/sqrt(boot_coef(:,in_signal_idx*2).^2 + boot_coef(:,1+in_signal_idx*2).^2 + boot_coef(:,2+in_signal_idx*2).^2));
            tab_append{in_signal_idx}.([prefix move_corr 'PD_zAng'])(uid)=circ_mean(z_ang);
        end
        %handle wrap around problems:
        centeredDirs=minusPi2Pi(dirs-circ_mean(dirs));
        
        if bootForTuning
            tab_append{in_signal_idx}.([prefix move_corr 'PDCI'])(uid)=prctile(centeredDirs,[2.5 97.5])+circ_mean(dirs);
        end
        
        % moddepth is poorly defined for GLM context, but for this case, let's use sqrt(sum(squares))
        moddepths_2D = sqrt(sum(boot_coef(:,(2*in_signal_idx):(2*in_signal_idx+1)).^2,2));
        tab_append{in_signal_idx}.([prefix move_corr 'Moddepth'])(uid)= mean(moddepths_2D);
        if(doing_3D_pd)
            moddepths_3D = sqrt(sum(boot_coef(:,(2*in_signal_idx):(2*in_signal_idx+2)).^2,2));
            tab_append{in_signal_idx}.([prefix move_corr 'Moddepth3D'])(uid)= mean(moddepths_2D);
        end
        
        if bootForTuning
            tab_append{in_signal_idx}.([prefix move_corr 'ModdepthCI'])(uid,:)= prctile(moddepths_2D,[2.5 97.5]);

            % get scrambled moddepths and check tuning
            scramble_moddepths = sqrt(sum(scramble_coef(:,(2*in_signal_idx):(2*in_signal_idx+1)).^2,2));
            scramble_high = prctile(scramble_moddepths,95);
            tab_append{in_signal_idx}.([prefix move_corr 'Tuned'])(uid,:) = (mean(moddepths_2D) > scramble_high);
            % diagnostic info
            if do_plot
                figure(h{in_signal_idx})
                scatter(scramble_moddepths,uid*ones(size(scramble_moddepths,1),1),[],'k','filled')
                hold on
                scatter(mean(moddepths_2D),uid,[],'r','filled')
            end
        end
        
    end
    if verbose
        fprintf('  Bootstrapping GLM PD computation %d of %d (ET=%f s)\n',uid,size(response_var,2),toc(unit_tic))
    end
end
starter = makeNeuronTableStarter(trial_data,params);
pdTable = tab_append{:};
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
