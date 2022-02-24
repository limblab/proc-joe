%Plot a few consecutive frames on the whole arm for either RT2D and RT3D task
%kinematics data: 
wa_f_st = 5000;                 %whole arm frame start
wa_f_len = 30;                  %whole arm frame length (frames)
wa_f_end = wa_f_st + wa_f_len; %whole arm frame end
marker_size = 75;
marker_width = 3;
line_width = 0.5;

shoulder_col = 1:3;
elbow2_col = 7:9;
wrist2_col = 13:15;
hand3_col = 22:24;

one_in_n_frames = 1;

shoulder_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , shoulder_col);
elbow2_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , elbow2_col);
wrist2_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , wrist2_col);
hand3_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , hand3_col);

shoulder_spd = calculateSpeed(shoulder_pos,td_list{i_task}.bin_size);
elbow2_spd = calculateSpeed(elbow2_pos,td_list{i_task}.bin_size);
wrist2_spd = calculateSpeed(wrist2_pos,td_list{i_task}.bin_size);
hand3_spd = calculateSpeed(hand3_pos,td_list{i_task}.bin_size);

%since the hand speed is roughly less than 2m/s, I'm just gonna set
%the color map to 200 colors, and plot the colors accordingly
shoulder_clr = round(shoulder_spd.*1000);
elbow2_clr = round(elbow2_spd*1000);
wrist2_clr = round(wrist2_spd*1000);
hand3_clr = round(hand3_spd*1000);

shoulder_clr(shoulder_clr > 200) = 200;
elbow2_clr(elbow2_clr > 200) = 200;
wrist2_clr(wrist2_clr > 200) = 200;
hand3_clr(hand3_clr > 200) = 200;

myColorMap = parula(250);
figure();
%pbaspect([2 1 1.5]);
view(3)
view(50,30);
hold on
grid on
%plot the dots of the markers/landmarks

%plot the link between the dots
for i = 1:length(shoulder_pos) %iterate through each frame
    %plot each single landmark on each frame
    %scatter3(shoulder_pos(i,1),shoulder_pos(i,2),shoulder_pos(i,3),marker_size,'filled','CData',myColorMap(shoulder_clr(i),:));
    %scatter3(elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3),marker_size,'filled','CData',myColorMap(elbow2_clr(i),:));
    %scatter3(wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3),marker_size,'filled','CData',myColorMap(wrist2_clr(i),:));
    %scatter3(hand3_pos(i,1),hand3_pos(i,2),hand3_pos(i,3),marker_size,'filled','CData',myColorMap(hand3_clr(i),:));

    %plot the lines between adjacent body landmarks for each frame
    shoulder_elbow = [shoulder_pos(i,1),shoulder_pos(i,2),shoulder_pos(i,3);
                      elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3)];
    elbow_wrist = [elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3);
                   wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3)];
    wrist_hand = [wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3);
                  hand3_pos(i,1),hand3_pos(i,2),hand3_pos(i,3)];
    plot3(shoulder_elbow(:,1),shoulder_elbow(:,2),shoulder_elbow(:,3),'LineWidth',line_width,'Color','k');
    plot3(elbow_wrist(:,1),elbow_wrist(:,2),elbow_wrist(:,3),'LineWidth',line_width,'Color','k');
    plot3(wrist_hand(:,1),wrist_hand(:,2),wrist_hand(:,3),'LineWidth',line_width,'Color','k');

    %plot the lines (link) for each landmark between two adjacent
    %frames
    if i < length(shoulder_pos)
        adjacent_shoulder = [shoulder_pos(i,1),shoulder_pos(i,2),shoulder_pos(i,3);
                             shoulder_pos(i+1,1),shoulder_pos(i+1,2),shoulder_pos(i+1,3)];
        adjacent_elbow = [elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3);
                          elbow2_pos(i+1,1),elbow2_pos(i+1,2),elbow2_pos(i+1,3)];
        adjacent_wirst = [wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3);
                          wrist2_pos(i+1,1),wrist2_pos(i+1,2),wrist2_pos(i+1,3)];
        adjacent_hand = [hand3_pos(i,1),hand3_pos(i,2),hand3_pos(i,3);
                         hand3_pos(i+1,1),hand3_pos(i+1,2),hand3_pos(i+1,3)];
        %DASPECT
        %daspect
        %pbaspect
        %pbaspect('auto')
%                 if i_task == 1 %RT2D
%                     pbaspect([5 1 4]);
         pbaspect([1 1 1]);
%                 else
%                     pbaspect([2 1 1.5]);
         pbaspect([1 1 1]);
%                 end
        %pbaspect([2 1 10]);
        %arrow3(adjacent_shoulder(1,:),adjacent_shoulder(2,:),'r2',4);
        %arrow3(adjacent_elbow(1,:),adjacent_elbow(2,:),'r2',3);
        %arrow3(adjacent_wirst(1,:),adjacent_wirst(2,:),'r2',4);           
        %arrow3(adjacent_hand(1,:),adjacent_hand(2,:),'r2',1.5);   
        %plot3(adjacent_shoulder(:,1),adjacent_shoulder(:,2),adjacent_shoulder(:,3),'Color','r');
        %plot3(adjacent_elbow(:,1),adjacent_elbow(:,2),adjacent_elbow(:,3),'Color','r');
        %plot3(adjacent_wirst(:,1),adjacent_wirst(:,2),adjacent_wirst(:,3),'Color','r');           
        %plot3(adjacent_hand(:,1),adjacent_hand(:,2),adjacent_hand(:,3),'Color','r');          
    end

end
if i_task == 1
    xlim([-30,30]);
    ylim([-30,30]);
    zlim([-30,5]);
else
    xlim([-15,20]);
    ylim([-15,20]);
    zlim([-30,5]);
end

xlabel('X Axis');
ylabel('Y Axis');
zlabel('Z Axis');