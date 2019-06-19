function [] = plotReaches(td,trial_idx,start_field,end_field,window,colors)
    % plot reaches (trial_idx) between start_field+window(1) and end_field + window(2)
    % colors is either a 1x3 matrix, or a numel(trial_idx)x3 matrix
    
    if(size(colors,1) == 1)
        colors = repmat(colors,numel(trial_idx),1);
    end
    
    color_counter = 1;
    for t = trial_idx
        plot(td(t).pos((td(t).(start_field)+window(1)):(td(t).(end_field)+window(2)),1),...
            td(t).pos((td(t).(start_field)+window(1)):(td(t).(end_field)+window(2)),2),'color',colors(color_counter,:),'linewidth',2);
        hold on
        color_counter = color_counter + 1;
    end

end

