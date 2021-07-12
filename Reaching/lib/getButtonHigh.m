function [trial_start, trial_end, mask] = getButtonHigh(td, varname)

    trial_start=find(diff(td.(varname)-mean(td.(varname))>3)>.5);
    trial_end=find(diff(td.(varname)-mean(td.(varname))<-3)>.5);

    % remove short trials
    trial_len = trial_end-trial_start;
    keep_trial = trial_len > mean(trial_len) - 2*std(trial_len);
    
    trial_start = trial_start(keep_trial==1);
    trial_end = trial_end(keep_trial==1);
    
    mask = zeros(numel(td.(varname)),1);
    
    keep_flag = 0;
    for i = 1:numel(mask)
        if(any(i == trial_start))
            keep_flag = 1;
        elseif(any(i == trial_end+1))
            keep_flag = 0;
        end
        mask(i) = keep_flag;
    end
end