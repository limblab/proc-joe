function [output_data] = reachingFiringRates(td_list,task_list)


    output_data = [];
    
    if(strcmpi(task_list{1},'RT')==1)
        td_plane = td_list{1};
        td_freereach = td_list{2};
    else
        td_plane = td_list{2};
        td_freereach = td_list{1};
    end



    % compare FR between plane and freereach for each neuron

    max_fr = [prctile(td_plane.LeftS1_FR,95,1);prctile(td_freereach.LeftS1_FR,95,1)];
    
    
    figure();
    histogram(max_fr(2,:)-max_fr(1,:),10);
    



end