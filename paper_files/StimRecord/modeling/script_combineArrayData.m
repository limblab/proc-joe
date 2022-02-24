%% combine data from multiple sources array datas
% Need: binned spike counts per trial, averaged binned counts, num stims, 
% amplitude, channel stim, channel rec, distance, bank, kin data as well

% Really need to get to an X and Y matrix
% Y -- average binned response between 1 and 5ms
% X -- distance, amplitude

%% NEED TO IMPLEMENT
% percent responding cells as a function of distance (use kstat to define
%   responsive) DONE
% inhibitory duration for each example
% peak latency
% onset latency

input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\';
input_data.mapFileName = {'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp',...
    'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'};
input_data.monkey_names = {'Chips','Han'};
input_data.poststim_window = [1,5];
input_data.prestim_window = [-10,-3];
input_data.wave_duration = 0.45;

input_data.max_p_val = 1; % throw away data above this value
input_data.p_val = 0.01; % determines responsiveness
% load all array datas, parse into an X and Y matrix
[data] = combineArrayData(input_data);


%% setup GLM
X = [data.distance/1000,data.amplitude];
Y = [data.response];

% remove below 0 Y values and X values
X(Y <= 0,:) = [];
Y(Y <= 0) = [];
Y(X(:,1) <= 0) = [];
X(X(:,1) <= 0,:) = [];


%% make GLM
num_evoked_glm = fitglm(X,Y,'Distribution','Poisson')

x_data_distance = 0:0.01:5;
y_data_distance = exp(num_evoked_glm.Coefficients.Estimate(1)+num_evoked_glm.Coefficients.Estimate(2)*x_data_distance + num_evoked_glm.Coefficients.Estimate(3)*mean(X(:,2)));
x_data_amplitude = 10:10:60;
y_data_amplitude = exp(num_evoked_glm.Coefficients.Estimate(1)+num_evoked_glm.Coefficients.Estimate(2)*mean(X(:,1)) + num_evoked_glm.Coefficients.Estimate(3)*x_data_amplitude);

figure()
% subplot(2,1,1)
plot(X(:,1),Y,'.','markersize',12);
hold on
plot(x_data_distance,y_data_distance,'r--','linewidth',2)
formatForLee(gcf)
xlabel('Distance (mm)');
ylabel('Response prob (# spikes/stim/cell)');
set(gca,'fontsize',14)

%% bin % responsive based on distance
    
d_dist = 250; % um bins
min_dist = 0;
max_dist = 5000;
bin_idx = [];
dist_bins = min_dist:d_dist:max_dist;

for d = 1:numel(data.dists)
    bin_idx_temp = find(data.dists(d) >= dist_bins,1,'last');
    if(~isempty(bin_idx_temp))
        bin_idx(d,1) = bin_idx_temp;
    else
        bin_idx(d,1) = [];
    end
end

num_responsive = zeros(max(bin_idx),1);
total_cells = zeros(max(bin_idx),1);
bE = d_dist*unique(bin_idx);

for d = 1:numel(data.dists)
    if(bin_idx(d) > 0)
        num_responsive(bin_idx(d)) = num_responsive(bin_idx(d)) + data.num_responsive(d);
        total_cells(bin_idx(d)) = total_cells(bin_idx(d)) + data.total_cells(d);
    end
end

num_responsive(total_cells==0) = [];
total_cells(total_cells==0) = [];

figure();
plot(bE/1000,num_responsive./total_cells,'.','markersize',12)
formatForLee(gcf)
xlabel('Distance (mm)');
ylabel('% responding cells');
set(gca,'fontsize',14)
ax = gca;
ax.YLim(1) = 0;
% fit data with exponential
percent_respond_glm = fitglm(bE/1000,num_responsive./total_cells,'Distribution','Poisson')
x_data = 0:0.1:5;
y_data = exp(percent_respond_glm.Coefficients.Estimate(1) + percent_respond_glm.Coefficients.Estimate(2)*x_data);
hold on
plot(x_data,y_data,'--r','linewidth',1.5)

%% estimate evoked # of spikes
% units = mm
dr = 0.001;
max_r = 6;
r_data = 0:dr:max_r;
a_data = 1:1:100;

C = @(a)exp(num_evoked_glm.Coefficients.Estimate(1) + num_evoked_glm.Coefficients.Estimate(3)*a);
rho = 75000; % Tehovnik 1995, cells/mm^3

f_num_spikes = @(r) exp(num_evoked_glm.Coefficients.Estimate(2)*r).*(r.^2)*dr;

f_percent_respond = @(r)1;% exp(percent_respond_glm.Coefficients.Estimate(1) + percent_respond_glm.Coefficients.Estimate(2)*r);

num_spikes_cum = [];
for a = a_data
    num_spikes = C(a)*rho*4*pi*f_percent_respond(r_data).*f_num_spikes(r_data);
    num_spikes_cum(end+1,:) = cumsum(num_spikes);
end

figure
plot(r_data,num_spikes_cum([10,50],:),'linewidth',2)
formatForLee(gcf)
xlabel('Distance (mm)')
ylabel('Cumulative # indirect spikes')
set(gca,'fontsize',14)

figure
plot(a_data,num_spikes_cum(:,end)/4800,'linewidth',1.5)
formatForLee(gcf)
xlabel('Amplitude (\muA)')
ylabel('Number evoked spikes or activated cells');
set(gca,'fontsize',14)
%% directly activated # cells based on Stoney's work
rho = 75000; % Tehovnik 1995, cells/mm^3
k = 272; % uA/mm^2; for low threshold neurons
i = a_data;
num_directly_activated_cells = 4/3*pi*(i/k).^(3/2)*rho;
hold on
plot(i,num_directly_activated_cells,'linewidth',1.5,'color','r')