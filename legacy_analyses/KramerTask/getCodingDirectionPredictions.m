function [output_data] = getCodingDirectionPredictions(neuro_data,cd_data,input_data)

    % project data from data onto coding direction
    output_data = [];
    figure();
    f.Name = '';
    
    x_data = 0:1:(size(neuro_data.spike_data,3)-1);
    cd_proj = zeros(size(neuro_data.spike_data,1),size(neuro_data.spike_data,3));    
    
    for trial = 1:numel(neuro_data.bump_dirs) % for each trial, plot projection as trial progresses
        cd_proj(trial,:) = (squeeze(neuro_data.spike_data(trial,:,:))'*cd_data.mean_cd)';
%         cd_proj(trial,:) = sum((squeeze(neuro_data.spike_data(trial,:,:)).*cd_data.mean_cd)',2);
    end

    % only use bump directions within a specified range
    if(isfield(input_data,'bump_range'))
        bump_mask = neuro_data.bump_dirs > input_data.bump_range(1) & neuro_data.bump_dirs < input_data.bump_range(2);
    else
        bump_mask = ones(size(data.bump_dirs));
    end
    
    % plot 0 deg projections
    plot(x_data,cd_proj(bump_mask == 1 & neuro_data.reached_0_deg==1,:),'b')
    hold on
    % plot 180 deg projections
    plot(x_data,cd_proj(bump_mask == 1 & neuro_data.reached_0_deg==0,:),'r')
    
    output_data.cd_proj = cd_proj;
    output_data.used_mask = bump_mask;
    
%     axis = 1;
%     bump_range = [80,100];
%     pred = mean(neurometric_data{axis}.spike_data(:,:,1:100),3)*cd_data{axis}.mean_cd;
%     pred(pred < 0) = 0;
%     pred(pred > 0) = 1;
%     
%     sum(pred(cd_data{axis}.training_trial_mask==1) == neurometric_data{axis}.reached_0_deg(cd_data{axis}.training_trial_mask==1))/...
%         numel(pred(cd_data{axis}.training_trial_mask==1))
% 
%     sum(pred(cd_data{axis}.training_trial_mask==1 & neurometric_data{axis}.bump_dirs > bump_range(1) & neurometric_data{axis}.bump_dirs < bump_range(2)) ...
%         == neurometric_data{axis}.reached_0_deg(cd_data{axis}.training_trial_mask==1 & neurometric_data{axis}.bump_dirs > bump_range(1) & neurometric_data{axis}.bump_dirs < bump_range(2)))/...
%         numel(pred(cd_data{axis}.training_trial_mask==1 & neurometric_data{axis}.bump_dirs > bump_range(1) & neurometric_data{axis}.bump_dirs < bump_range(2)))
%     



end