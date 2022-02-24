td = splitTD(...
            trial_data,...
            struct(...
                'split_idx_name','idx_startTime',...
                'linked_fields',{{...
                    'trialID',...
                    'result',...
                    'bumpDir',...
                    'tgtDir',...
                    }},...
                'start_name','idx_startTime',...
                'end_name','idx_endTime',...
                'extra_time',[-2,2]));
[~,td] = getTDidx(td,'result','R');
td = reorderTDfields(td);

td((~isnan([td.bumpDir]))) = [];


td = trimTD(td,{'idx_goCueTime',40},{'idx_goCueTime',260});
td = removeBadTrials(td);
td(isnan([td.tgtDir])) = [];
td_avg = trialAverage(td,{'tgtDir'});

%% set a bin size, save emg and tgt dir info for each reach direction
bin_sizes = [5,10,20,25,50]/1000; %ms
filename = 'D:\Lab\Data\StimPDs\Han_20170203_emg_reach_';
for i_bin = 1:numel(bin_sizes)
    td_avg_bin = binTD(td_avg,floor(bin_sizes(i_bin)/td_avg(1).bin_size));
    for i_tgt = 1:numel(td_avg)
        data = td_avg_bin(i_tgt).emg;
        data = data - data(1,:);
        csvwrite([filename,num2str(td_avg_bin(i_tgt).tgtDir),'deg_',num2str(bin_sizes(i_bin)*1000),'ms.txt'],data);
    end
end




%%
for i_musc = 1:size(td_avg(1).emg,2)
    figure(); hold on
    for i = 1:numel(td_avg)
        plot(td_avg(i).emg(:,i_musc))
    end
end

