%% this script is meant to check 2D tracking across all cameras for each marker
% first, provide folderpath to videos (and DLC outputs) plus the video
% numbers to use


folderpath = 'D:\Lab\Data\DLC_videos\Han_20201204_rwFreeReach\'; % include \ at end
use_filtered = 1; % use filtered data or not. 
vid_nums = [5,6,7,8]; % idx at end of video, leave empty if you want all of the videos in the folder


if(use_filtered)
    file_str = [folderpath,'*_filtered.csv'];
else
    file_str = [folderpath,'*0.csv'];
end

files = dir(file_str);

% keep only files based on vid_nums
if(~isempty(vid_nums))
    keep_mask = zeros(size(files));
    % vid num is before DLC
    for i_file = 1:numel(files)
        underscore_idx = strfind(files(i_file).name,'_');
        DLC_idx = strfind(files(i_file).name,'DLC');
        underscore_idx = underscore_idx(find(underscore_idx < DLC_idx,1,'last'));
        num = str2num(files(i_file).name(underscore_idx+1:DLC_idx-1));
        
        if(any(num==vid_nums))
            keep_mask(i_file)=1;
        end
    end
    files = files(keep_mask==1);
end


% load in each file and store marker positions, likelihoods, etc. for each
% video

for i_file = 1:numel(files)
    % Reads in table. 
    tbl = readtable([files(i_file).folder,filesep,files(i_file).name],'headerlines',1);
    % sanitize tbl. Remove first col (called 'bodyparts'). Also remove first
    % row, which got some extra headers. for each body part, it's x, y, then likelihood.
    tbl.bodyparts = [];
    tbl(1,:) = [];

    % convert each entry to a double from a char. Also change column names
    % (goes x, then y, then likelihood for each marker)
    for i_col = 1:size(tbl,2)
    % convert 
        var_name = tbl.Properties.VariableNames{i_col};
        tbl.(var_name) = str2double(tbl.(var_name));

        % change var_name
        iter = mod(i_col-1,3);
        switch iter
            case 0
                % add _x
                tbl.Properties.VariableNames{i_col} = [var_name,'_cam_',num2str(i_file),'_x'];
            case 1
                % remove _1 and add _y
                tbl.Properties.VariableNames{i_col} = [var_name(1:end-2),'_cam_',num2str(i_file),'_y'];
            case 2
                % remove _2 and add _likelihood
                tbl.Properties.VariableNames{i_col} = [var_name(1:end-2),'_cam_',num2str(i_file),'_likeli'];
        end
    end

    % put table in larger table
    if(i_file==1)
        tracking_data = tbl;
    else
        if(size(tracking_data,1) > size(tbl,1))
            tracking_data = tracking_data(1:size(tbl,1),:);
        elseif(size(tbl,1) > size(tracking_data,1))
            tbl = tbl(1:size(tracking_data,1),1);
        end
        tracking_data = [tracking_data,tbl];
    end
end
% get marker names
var_names = tracking_data.Properties.VariableNames;

for i_name = 1:numel(var_names)
    underscore_idx = strfind(var_names{i_name},'_');
    var_names{i_name} = var_names{i_name}(1:underscore_idx-1);
end

marker_names = unique(var_names,'stable');

% remove points from marker names
keep_mask = ones(size(marker_names));
for i_marker = 1:numel(marker_names)
    if(strfind(marker_names{i_marker},'point'))
        keep_mask(i_marker) = 0;
    end
end
marker_names = marker_names(keep_mask==1);

%% now that we have loaded in the 2D tracking data for each camera
% plot 2 best likelihood over time for each marker
best_like = nan(size(tracking_data,1),numel(marker_names),2);
var_names = tracking_data.Properties.VariableNames;

figure();
for i_marker = 1:numel(marker_names)
    subplot(ceil(numel(marker_names)/2),2,i_marker)
    % find marker entries in tracking data, and likelihood entries
    keep_mask = ~cellfun(@isempty,strfind(var_names,marker_names{i_marker})) & ~cellfun(@isempty,strfind(var_names,'likeli'));

    likeli_data = table2array(tracking_data(:,keep_mask==1));
    
    [best_like(:,i_marker,1),max_idx] = max(likeli_data,[],2);
    
    for i_like = 1:size(likeli_data,1)
        likeli_data(i_like,max_idx(i_like)) = -1;
    end
    [best_like(:,i_marker,2)] = max(likeli_data,[],2);
    
    plot(squeeze(best_like(:,i_marker,:)))
    title(marker_names{i_marker})
end

%%
sum(best_like(:,:,2) > 0.7)/size(best_like,1)
% plot speed of points over time (2D speed)

% plot distance between sets of points over time, and how that distance
% changes over time















