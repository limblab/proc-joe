function [ temp_td ] = removeTrials( trial_data )
% this function removes trials from a trial data, resulting in a 1x1 struct
% with all of the data. This is a quick fix while the 'noTrials' option in
% parseFileByTrial does not work....

    temp_td = trial_data(:,1);
    
    % go through each field and deal with it. Meta data stays the same
    % (monkey, date, task). all idx_Times need to be updated. Append pos,
    % vel, acc, force, dlc_data, _spikes, _unit_guide if they exist
    
    td_names = fieldnames(trial_data);
    
    
    for i_trial = 2:numel(trial_data)
        for i_name = 1:numel(td_names)
            % if idx_ is present, increment and append to larger list 
            % if it's a cell, leave alone (likely a name array)
            % if it's _unit_guide, leave alone
            % else append to corresponding temp_td array

            if(~isempty(strfind(td_names{i_name},'idx_')))
                % increment then append
                temp_val = trial_data(i_trial).(td_names{i_name}) + temp_td.idx_endTime(i_trial-1);
                temp_td.(td_names{i_name}) = [temp_td.(td_names{i_name}); temp_val];
            elseif(~isempty(strfind(td_names{i_name},'_ts')))
                % for the spikes, we need to go through each cell and
                % append
                for i_unit = 1:numel(temp_td.(td_names{i_name}))
                    temp_td.(td_names{i_name}){i_unit} = [temp_td.(td_names{i_name}){i_unit}; trial_data(i_trial).(td_names{i_name}){i_unit}];
                end
            elseif(iscell(trial_data(1).(td_names{i_name})) || ~isempty(strfind(td_names{i_name},'_unit_guide')) || ...
                    ~isempty(strfind(td_names{i_name},'monkey')) || ~isempty(strfind(td_names{i_name},'date')) ||  ~isempty(strfind(td_names{i_name},'task')))
                % do nothing!
            
            else
                % append
                temp_td.(td_names{i_name}) = [temp_td.(td_names{i_name}); trial_data(i_trial).(td_names{i_name})];
            end
        end
    end


end

