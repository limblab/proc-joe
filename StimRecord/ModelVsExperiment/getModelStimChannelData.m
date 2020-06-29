function [mdl_data_all, diam_all, amp_all, cell_id_all] = getModelStimChannelData(input_data)
% input data contains f
    %folderpath: data directory
    %diam_list : 1, 2, 3 for diameter
    %amp_list : amplitude of stimulation
    %get_axon_dendrite_locs : bool, do you want locations of all axon or
            %dendrites?
    %cell_id_list : 6,11,16,21 corresponds to type of cell
    curr_fpath = cd;
    
    diam_all = [];
    amp_all = [];
    cell_id_all = [];
    
    mdl_data_all = [];
    cond_counter = 1;
    for i_diam = 1:numel(input_data.diam_list)
        for i_amp = 1:numel(input_data.amp_list)
            for i_cell_id = 1:numel(input_data.cell_id_list)
                diam_all(end+1,1) = input_data.diam_list(i_diam);
                amp_all(end+1,1) = input_data.amp_list(i_amp);
                cell_id_all(end+1,1) = input_data.cell_id_list(i_cell_id);
                
                % make coord data
                coord_data = [];
                coord_data.diam = input_data.diam_list(i_diam);
                coord_data.amp = input_data.amp_list(i_amp);
                coord_data.cell_id = input_data.cell_id_list(i_cell_id);
                coord_data.get_axon_dendrite_locs = input_data.get_axon_dendrite_locs;
                coord_data.folderpath = input_data.folderpath;
                
                % get location and spike data
                data_path = [coord_data.folderpath, 'Diam',num2str(coord_data.diam),'_',num2str(coord_data.amp),'uA\'];
                cd(data_path);
                location_data = getTrueCoordinates(coord_data);
                load(['data_soma',num2str(coord_data.cell_id)]);
                
                mdl_data_all(cond_counter).soma = location_data.soma;
                mdl_data_all(cond_counter).axon = location_data.axon;
                mdl_data_all(cond_counter).loc_all = location_data.loc_all;
                mdl_data_all(cond_counter).spike_data = data_soma;
                mdl_data_all(cond_counter).diam = coord_data.diam;
                mdl_data_all(cond_counter).amp = coord_data.amp;
                mdl_data_all(cond_counter).cell_id = coord_data.cell_id;
                mdl_data_all(cond_counter).folderpath = data_path;
                mdl_data_all(cond_counter).get_axon_dendrite_locs = coord_data.get_axon_dendrite_locs;
                
                
                % update counter
                cond_counter = cond_counter + 1;
            end
        end
    end    
    
    % go back to original folder
    cd(curr_fpath);

end




function [ output_data ] = getTrueCoordinates( input_data )
% input data contains f
    %diam : 1, 2, 3 for diameter
    %amp : amplitude of stimulation
    %get_axon_dendrite_locs : bool, do you want locations of all axon or
            %dendrites?
    %cell_id : 6,11,16,21 corresponds to type of cell
    
% outputs location data for each cell
    
    load('realx.dat') %x-coordinate (in um) of randomly seeded soma within spherical volume
    load('realy.dat') %y-coordinate (in um) of randomly seeded soma within spherical volume
    load('realz.dat') %z-coordinate (in um) of randomly seeded soma within spherical volume

    load('realang.dat') %Angle (in radian) of random rotation of neuron in azimuthal direction 

    soma_coord=load(['soma_coord_' num2str(input_data.cell_id) '.dat']); %x,y,z coordinates (in um) of somatic section (of a single neuron) 

    %To get true location of somatic section of each neuron within the spherical
    %volume
    soma = []; axon = []; loc_all = [];
    
    for k=1:length(realx)

        net_ad_x=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)))-realx(k);
        net_ad_y=realy(k)+soma_coord(2)-realy(k);
        net_ad_z=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)))-realz(k);

        x1=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)));
        y1=realy(k)+soma_coord(2);
        z1=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)));

        soma(k).coord(1)=x1-net_ad_x;
        soma(k).coord(2)=y1-net_ad_y;
        soma(k).coord(3)=z1-net_ad_z;
        soma(k).meta = 'xyz';

    end

    if(input_data.get_axon_dendrite_locs==1)
        intx=load(['intx_' num2str(cell_id) '.dat']); %x-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron 
        inty=load(['inty_' num2str(cell_id) '.dat']); %y-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron
        intz=load(['intz_' num2str(cell_id) '.dat']); %z-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron


        x_axon=load(['x_axon_' num2str(cell_id) '.dat']); %x-coordinates (in um) of all axonal sections (of a single neuron)
        y_axon=load(['y_axon_' num2str(cell_id) '.dat']); %y-coordinates (in um) of all axonal sections (of a single neuron)
        z_axon=load(['z_axon_' num2str(cell_id) '.dat']); %z-coordinates (in um) of all axonal sections (of a single neuron)


    %     To get true location of axon section of each neuron within the spherical
    %     volume
        for k=1:length(realx)

            net_ad_x=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)))-realx(k);
            net_ad_y=realy(k)+soma_coord(2)-realy(k);
            net_ad_z=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)))-realz(k);

            for i=1:length(x_axon)

                x1=((realx(k)+x_axon(i))*cos(realang(k)))-((realz(k)+z_axon(i))*sin(realang(k)));
                y1=realy(k)+y_axon(i);
                z1=((realx(k)+x_axon(i))*sin(realang(k)))+((realz(k)+z_axon(i))*cos(realang(k)));

                axon(k).coord(i,1)=x1-net_ad_x;
                axon(k).coord(i,2)=y1-net_ad_y;
                axon(k).coord(i,3)=z1-net_ad_z;
            end
        end



    %     To get true location of all neural sections (i.e. dendrite, axon and soma) of each neuron within the spherical
    %     volume
        for k=1:length(realx)

            net_ad_x=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)))-realx(k);
            net_ad_y=realy(k)+soma_coord(2)-realy(k);
            net_ad_z=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)))-realz(k);

            for i=1:length(intx)

                x1=((realx(k)+intx(i))*cos(realang(k)))-((realz(k)+intz(i))*sin(realang(k)));
                y1=realy(k)+inty(i);
                z1=((realx(k)+intx(i))*sin(realang(k)))+((realz(k)+intz(i))*cos(realang(k)));

                loc_all(k).coord(i,1)=x1-net_ad_x;
                loc_all(k).coord(i,2)=y1-net_ad_y;
                loc_all(k).coord(i,3)=z1-net_ad_z;
            end
        end
    end
    
    % package outputs
    output_data = [];
    
    output_data.soma = soma;
    output_data.axon = axon;
    output_data.loc_all = loc_all;
    output_data.diam = input_data.diam;
    output_data.folderpath = input_data.folderpath;
    output_data.amp = input_data.amp;
    output_data.cell_id = input_data.cell_id;
    output_data.get_axon_dendrite_locs = input_data.get_axon_dendrite_locs;
    
    
end

