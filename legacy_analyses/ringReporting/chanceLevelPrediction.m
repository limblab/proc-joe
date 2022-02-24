%% simulate biomimetic data to get an idea for what MCE is
num_reaches = 25;
num_dir = 20;
std_max = 5;
directions = 2*pi*rand(num_dir,1);
std = std_max*rand(num_dir,1);

reaches = [];
reaches_shuffled = [];
MCE = [];
MCE_suffled = [];
for i = 1:numel(directions)
    reaches(i,:) = rand(1,num_reaches)*std(i) - std(i)/2 + directions(i);
    MCE(i) = mean_cosine_error(reaches(i,:),directions(i));
end

reaches = wrapTo2Pi(reaches);

figure();
for i = 1:numel(directions)
    plot(directions(i),reaches(i,:),'k.','markersize',12)
    hold on
end

xlim([-1,1+2*pi])
ylim([-1,1+2*pi])



reaches_shuffled = reshape(reaches(randperm(numel(reaches))),size(reaches));
MCE_shuffled = mean_cosine_error(reaches_shuffled,directions);

% for i = 1:numel(directions)
%     plot(directions(i),reaches_shuffled(i,:),'r.','markersize',12)
%     hold on
% end

figure();
plot(directions,MCE,'k.','markersize',12)
hold on
plot(directions,MCE_shuffled,'r.','markersize',12)
ylim([0,2])