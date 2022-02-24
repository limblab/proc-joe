function [output_data] = getArtifactMetrics(artifact)
    % takes in a matrix of artifact data. # stimulations by # channels by
    % data points (3D array)
    
    % outputs metrics for each channel
    % 1. Time off rails (last data point above some value)
    % 2. Settling time??
    
    max_amp = 1000;
    output_data.idx_off_rails = nan(size(artifact,1),size(artifact,2));
    
    for stim = 1:size(artifact,1)
        for chan = 1:size(artifact,2)
            temp = find(abs(artifact(stim,chan,:)) > max_amp,1,'last');
            if(~isempty(temp))
                if(temp > 100) % presample
                    output_data.idx_off_rails(stim,chan) = temp;
                end
            end
        end
    end
    
   

end