function [heatmapDataPD,alphaArray] = getHeatmapDataPD(td_all,pd_all,mapData)

    %create array that you're going to pass to imagesc
    heatmapDataPD = nan(10,10);
    alphaArray = zeros(10,10);

    %ensure all angles are between 0-360deg
    PDtemp = radtodeg(pd_all.velPD);
    
    %put the angles in their respective locations in the array
    for pd_idx =1:numel(PDtemp)
        for map_data_idx = 1:numel(mapData.chan)
            if td_all(1).LeftS1_unit_guide(pd_idx) == mapData.chan(map_data_idx)
                heatmapDataPD((11-mapData.row(map_data_idx)),mapData.col(map_data_idx)) = PDtemp(pd_idx);
                alphaArray(11-mapData.row(map_data_idx),mapData.col(map_data_idx))=1;
            end
        end
    end
    

end

