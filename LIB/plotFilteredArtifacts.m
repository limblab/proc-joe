function [ figHandle ] = plotFilteredArtifacts( artifact, b, a )
    [bOld,aOld]=butter(6,500/(30000/2),'high');
    
    xData = (0:1:size(artifact,2)-1)/30;
%     figHandle = figure()
    
    ax1=subplot(3,1,1)
    plot(xData,artifact','linewidth',1.5)
    ylabel('Voltage (\muV)')
    ylim([-1000,1000])
    xlim([0,3])
    formatForLee(gcf)
    
    
    ax2=subplot(3,1,2)
    plot(xData,fliplr(filter(bOld,aOld,fliplr(artifact)')')','linewidth',1.5)
    ylabel('Voltage (\muV)')
    ylim([-1000,1000])
    xlim([0,3])
    formatForLee(gcf)
    
    ax3=subplot(3,1,3)
    plot(xData,fliplr(filter(b,a,fliplr(artifact)')')','linewidth',1.5)
    ylim([-1000,1000])
    xlim([0,3])
    ylabel('Voltage (\muV)')
    xlabel('Time after stimulation onset (ms)')

    linkaxes([ax1,ax2]);
    
    formatForLee(gcf)
    
%     %% plot what would happen with forward filter or filtfilt
%     figure();
% %     artifact = [repmat(artifact(:,1),1,100),artifact];
% %     xData = [(-100:1:-1)/30,xData];
%     ax1 = subplot(3,1,1);
%     plot(xData,artifact)
%     ylim([-1000,1000])
%     xlim([0,3])
%     ax2 = subplot(3,1,2);
%     plot(xData,filter(b,a,artifact')')
%     ylim([-1000,1000])
%     xlim([0,3])
%     ax3 = subplot(3,1,3);
%     plot(xData,filtfilt(b,a,artifact')')
%     ylim([-1000,1000])
%     xlim([0,3])
%     linkaxes([ax1,ax2,ax3]);
%     
    
end

