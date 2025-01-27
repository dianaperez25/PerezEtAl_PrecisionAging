data_dir = '/Users/dianaperez/Desktop/FC_Parcels_333'; % the timecourses for the different parcels
output_dir = '/Users/dianaperez/Desktop/'; % where to store the results
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/Nifti/FC_Parcels_333/'; % the timecourses for the different parcels
output_dir = '/Users/dianaperez/Desktop/'; % where to store the results
%% OPTIONS
datasets = {'iNet-NU', 'Lifespan-FSU', 'iNet-FSU', 'Lifespan-NU'};%% 

% How many points to sample for "true" half
truehalf_datapts = 3808; %3808 -> ~70 min, 5454 -> ~100 min, 8181 -> ~150 min

% How much data to add at each step in minutes (how many data points that
% corresponds to is calculated below)
step_time = 2.5;
step_datapts = round((60*step_time)/1.1);% 136 -? ~2.5 min, 272 -> ~5 min, will add this number of frames each time it subsamples data

% How many iterations to run
iterations = 1000;

for d = 1:numel(datasets)
    dataset = datasets{d};
    t0=tic;
    %% VARIABLES
if strcmpi(dataset, 'Lifespan-NU')
    subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'}; %
    sessions = 5;
    %rgb colors for plotting results
    rgb_colors = [1 0 0;%LS02
                0, 1, 0;%LS03
                0, 0, 1;%LS05
                0, 1, 1;%LS08
                1, 0, 1;%LS11
                0.4660 0.6740 0.188;%LS14
                0.9290 0.6940 0.1250;%LS16
                0.4940 0.1840 0.5560];%LS17
elseif strcmpi(dataset, 'iNet-NU')
    subject = {'INET002', 'INET003', 'INET005', 'INET006','INET010',...
    'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
    'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
    'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
    'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
    'INET057', 'INET058', 'INET059', 'INET060', 'INET062', 'INET063',...
    'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %'INET061', 'INET001', 
    sessions = 4;
elseif strcmpi(dataset, 'iNet-FSU')
    subject = {'INET074', 'INET075', 'INET077', 'INET078', 'INET083', 'INET084',...
    'INET085', 'INET086', 'INET087', 'INET088',...
    'INET091', 'INET093', 'INET094', 'INET096', 'INET098', 'INET099', 'INET101',...
    'INET103', 'INET104', 'INET105', 'INET106', 'INET107', 'INET108', 'INET109',...
    'INET110', 'INET111', 'INET112', 'INET114', 'INET115', 'INET116', 'INET120',...
    'INET123', 'INET133', 'INET136', 'INET137', 'INET140', 'INET141', 'INET143',...
    'INET156', 'INET158', 'INET159', 'INET160', 'INET165', 'INET168'};
    sessions = 4;
elseif strcmpi(dataset, 'Lifespan-FSU')
    subject = {'LS31', 'LS32', 'LS33', 'LS39', 'LS43', 'LS44', 'LS45', 'LS46',...
    'LS47', 'LS54', 'LS61', 'LS62', 'LS63', 'LS68', 'LS70', 'LS71', 'LS72',...
    'LS76', 'LS77', 'LS79', 'LS85', 'LS89', 'LS94', 'LS108'};
    sessions = 5;
end

load(['/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/Reliability/' dataset '_Reliability_truehalf_3808_corrdata.mat'])
if strcmpi(dataset, 'Lifespan-FSU')
    means{1,8}=[];
end
figure;

if strcmpi(dataset, 'Lifespan-NU')

    for s = 1:numel(subject)
       plot(times_all(s,1:size(means{1,s},2)),means{1,s},'Color',rgb_colors(s,:),'LineWidth', 3)
       hold on
    end
else
    for s = 1:numel(subject)
       plot(times_all(s,1:size(means{1,s},2)),means{1,s}, 'LineWidth', 2)
       hold on
    end
end

mean_of_means = [];



for t = 1:size(times_all,2)
    tmp = [];
    for s = 1:numel(subject)
        if size(means{1,s},2)>=t
            tmp = [tmp;means{1,s}(t)];
        else
            continue;
        end
    end
    mean_of_means(t) = mean(tmp);
end
axis([0 100 .45 1])
yticks([.45:.10:.95])

if strcmpi(dataset, 'iNet-FSU')
    times_all = step_time:step_time:31*step_time;
else
    times_all = step_time:step_time:size(times_all,2)*step_time;
end
plot(times_all,smooth(smooth(smooth(smooth(mean_of_means(1:size(times_all,2)))))), ':', 'Color', [0,0,0], 'LineWidth',3) %average

plot([0 100], [.85 .85], 'color', [.5 .5 .5 .5], 'LineWidth',3, 'LineStyle', '--')
                       

ind_mark = 0;
for m = 1:length(mean_of_means)
    %disp(round(mean_of_means(m), 2))
    if mean_of_means(m) >= .85
        ind_mark = m;
        break;
    end
end
if strcmpi(dataset, 'iNet-NU')
    ind_mark = ind_mark + 1;
end
plot(times_all(ind_mark), .85, 'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','w','MarkerSize',10)
ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');

set(gca,'FontSize',16)

print(gcf,[output_dir dataset '_Reliability_truehalf_' num2str(truehalf_datapts) '.jpg'],'-dpng','-r300');
close all

end