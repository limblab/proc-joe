data = artifactData.artifact(:,:,100:end);

railTime = [];
for i = 33:size(data,1)
    for j = 1:size(data,2)
        dataCurr = data(i,j,:);
        dataCurr = squeeze(dataCurr - mean(dataCurr));
        
        % find smallest
        [small,smallIdx] = min(dataCurr);
        while(min(dataCurr(smallIdx+1:end)) <= small*0.999)
            [small,smallIdxStep] = min(dataCurr(smallIdx+1:end));
            smallIdx=smallIdx+smallIdxStep;
        end
        % find largest
        [big,bigIdx] = max(dataCurr);
        while(max(dataCurr(bigIdx+1:end)) >= big*0.999)
            [big,bigIdxStep] = max(dataCurr(bigIdx+1:end));
            bigIdx = bigIdx+bigIdxStep;
        end
        % pick between the 2
        time = max(smallIdx,bigIdx);
        
        railTime(end+1) = time;

    end
    
end

figure;
bE = (0.5:0.05:1.3)*30000/1000;
[bC,bE] = histcounts(railTime,bE);
bC = bC/sum(bC);
bE = bE/30000*1000; % in ms now
bar(bE(1:end-1)+(bE(2)-bE(1))/2,bC)
xlim([0 50/30000*1000])
ylim([0 1])
xlabel('Time after stimulation (ms)');
ylabel('Percent artifacts')
formatForLee(gcf);