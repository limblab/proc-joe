function [ psych_data ] = getPsychometricCurveData( td,input_data )


    %% seperate td based on target axis
    axis_angles = unique([td.tgtDir]);
    axis_angles = axis_angles(~isnan(axis_angles));
    
    % get psych data for actual data and all bootstraps
    for boot_idx = 1:input_data.num_bootstrap+1 % 1 extra to do work on actual data
        if(boot_idx == 1) % no random sample, use actual data
            td_boot = td;
        else % randomly sample data from td
            td_boot = datasample(td,numel(td),'Replace',true);
        end

        %% get data for each axis angle
        for axis_num = 1:numel(axis_angles)
            td_mask = [td_boot.tgtDir] == axis_angles(axis_num);
            td_axis = td_boot(td_mask==1);
            % handle the bumps
            input_data.isBump = 1; 
            input_data.stim_code = 0;
            input_data.axis_angle = axis_angles(axis_num);
            psych_data{boot_idx,axis_num} = getPsychDataWrapper(td_axis,input_data);

            psych_data{boot_idx,axis_num}.input_data = input_data;

            % handle the stims
            stim_codes = unique([td_axis.stimCode]);
            stim_codes = stim_codes(~isnan(stim_codes));
            for s = 1:numel(stim_codes)
                input_data.isBump = 0;
                input_data.stim_code = stim_codes(s);
                temp_psych_data = getPsychDataWrapper(td_axis,input_data);
                temp_psych_data.input_data = input_data;

                psych_data{boot_idx,axis_num}(end+1) = temp_psych_data;

            end
        end
        
        
    end

end

function [output_data] = getPsychDataWrapper(td,input_data)

    if(input_data.isBump)
        trial_mask = [td.isStimTrial] == 0 & [td.bumpMagnitude] > 0.1;
    else
        trial_mask = [td.isStimTrial] == 1 & [td.stimCode] == input_data.stim_code & [td.bumpMagnitude] > 0.1;
    end
    
    td_cond = td(trial_mask);
    td_all_trials = [td(trial_mask).trial_id];

    [temp_data] = getPsychData(td_cond, input_data);
    

    %% format output_data
    output_data.trial_ids = td_all_trials;
    output_data.bump_dirs = temp_data.bump_dirs;
    output_data.bump_correct = temp_data.bump_correct;
    output_data.bump_total = temp_data.bump_total;
    output_data.psych_curve_data = temp_data.psych_curve_data;
    output_data.psych_fit = temp_data.psych_fit;
    output_data.is_rewarded = temp_data.is_rewarded;
    output_data.bump_dir_list = temp_data.bump_dir_list;

end



function [output_data] = getPsychData(td_cond, input_data)

    bump_dirs = unique([td_cond.bumpDir]);
    bump_dirs = bump_dirs(~isnan(bump_dirs));
    bump_dirs = bump_dirs(bump_dirs <= 180);
    bump_correct = zeros(size(bump_dirs));
    bump_total = zeros(size(bump_dirs));
    
    bump_dir_list = [];
    is_rewarded_list = [];
    
    for t = 1:numel(td_cond)
        if(td_cond(t).bumpDir > 180)
            bump_idx = find(bump_dirs == -1*(td_cond(t).bumpDir-360));
            bump_dir_list(t) = -1*(td_cond(t).bumpDir-360);
        else
            bump_idx = find(bump_dirs == td_cond(t).bumpDir);
            bump_dir_list(t) = td_cond(t).bumpDir;
        end
        
        if(td_cond(t).idx_endTime - td_cond(t).idx_goCueTime < input_data.max_trial_time && ...
                td_cond(t).idx_endTime - td_cond(t).idx_goCueTime > input_data.min_trial_time)
            if(td_cond(t).result == 'R')
                bump_correct(bump_idx) = bump_correct(bump_idx) + 1;
                is_rewarded_list(t) = 1;
            else
                is_rewarded_list(t) = 0;
            end
            bump_total(bump_idx) = bump_total(bump_idx) + 1;
        end
        
        
    end

    psych_curve_data = bump_correct./bump_total;
    
    % 1-bump_percent_correct for dirs > 90 to make these figures
    psych_curve_data(bump_dirs > 90) = 1 - psych_curve_data(bump_dirs > 90);
    
    
    % fit psychometric curve
    fit_opts = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 -0.5 .1 100],...
        'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8,...
        'Weights',bump_total/sum(bump_total));
    ft = fittype('a+b*(erf(c*(x-d)))','options',fit_opts);
    
    psych_fit = [];
    try
        [psych_fit.fitObj,psych_fit.gof] = fit(bump_dirs', psych_curve_data', ft);
    catch
    end
    
    % format output_data
    output_data.is_rewarded = is_rewarded_list;
    output_data.bump_dir_list = bump_dir_list;
    output_data.bump_dirs = bump_dirs;
    output_data.bump_correct = bump_correct;
    output_data.bump_total = bump_total;
    output_data.psych_curve_data = psych_curve_data;
    output_data.psych_fit = psych_fit;
end
