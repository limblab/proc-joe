latencies = [];
for i = 1:size(arrayData{1}.stimData,1)
    for j = 1:size(arrayData{1}.stimData,2)
        for k = 1:numel(report.excitatoryUnits{i,j})
            if(arrayData{report.excitatoryUnits{i,j}(k)}.excitatoryLatency{i,j} > 0)
                figure();
                plot(arrayData{report.excitatoryUnits{i,j}(k)}.bC{i,j})
                xlim([100,300])
                hold on
%             latencies = [latencies;arrayData{report.excitatoryUnits{i,j}(k)}.excitatoryLatency{i,j}]; 
            end
        end
    end
end