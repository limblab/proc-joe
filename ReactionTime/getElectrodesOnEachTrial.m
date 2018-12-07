function [output_data] = getElectrodesOnEachTrial(input_data,opts)
% we have a list of electrodes and stim codes from MATLAB, need to match
% with the trial data (from behavior code)

%% configure opts and other parameters
    opts = configureOpts(opts);
    output_data = [];
    plots = [];
    map_data = loadMapFile(input_data.map_file_name);
    
%% translate stim_code_all back to c++ code indexing
    stim_code_matlab = input_data.stim_code_all' - 1; % translate back to c++ code
    stim_code_td = [[input_data.td_reward.stimCode]',[input_data.td_reward.trial_id]'];
    
    stim_code_td = stim_code_td(~isnan(stim_code_td(:,1)),:);
    EL_all = input_data.EL_all;

% %% align stim_code_all and stim_code_td    
% 
%     [align_data,lag] = xcorr(stim_code_matlab,stim_code_td(:,1));
%         
%     [~,lag_idx] = max(align_data);
%     stim_code_matlab_aligned = stim_code_matlab(max(1,lag(lag_idx)):end);
%     EL_all = input_data.EL_all(max(1,lag(lag_idx)):end);
    
%% go through each stim_code_td and find corresponding list of electrodes
% compute distance, store distance, find td_reward_idx based on trial id
% and store the corresponding rt as well
    output_data.EL_list = {};
    output_data.pos = {};
    output_data.rt = [];
    output_data.stim_code = [];
    
    matlab_idx = 1;
    for s = 1:size(stim_code_td,1)
        while(matlab_idx < numel(stim_code_matlab) && ...
                stim_code_matlab(matlab_idx) ~= stim_code_td(s,1))
            
            matlab_idx = matlab_idx + 1;
        end
        
        % reward_idx
        td_reward_idx = find([input_data.td_reward.trial_id] == stim_code_td(s,2));
        
        if(~isempty(td_reward_idx) && matlab_idx < numel(EL_all))
            % get EL data
            elec_list = EL_all{matlab_idx};
            output_data.EL_list{end+1} = elec_list;

            % get (x,y) pos of elec_list
            x_pos = []; y_pos = [];
            for e = 1:numel(elec_list)
                map_idx = find(map_data.chan == elec_list(e));
                x_pos(end+1,1) = map_data.row(map_idx);
                y_pos(end+1,1) = map_data.col(map_idx);
            end
            output_data.pos{end+1,1} = [x_pos,y_pos];
            
            % get rt data from input_data.data struct
            cueInfo_idx = find(input_data.data.cueInfo(stim_code_td(s,1)+1).td_idx_reward == td_reward_idx);
            if(~isempty(cueInfo_idx))
                output_data.rt(end+1,1) = input_data.data.cueInfo(stim_code_td(s,1)+1).rt(cueInfo_idx);
            else
                output_data.rt(end+1,1) = -1000;
                warning(['wtf',num2str(s)]);
            end
            
            output_data.stim_code(end+1,1) = stim_code_td(s,1);
        end
        
        % update our iterator
        matlab_idx = matlab_idx+1;

    end



end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldNames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldNames)
           if(isfield(opts,inputFieldNames{fn}))
               opts.(inputFieldNames{fn}) = optsInput.(inputFieldNames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    

end



            