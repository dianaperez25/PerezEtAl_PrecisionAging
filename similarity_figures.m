%% Script to make similarity matrix figures
% This script plots a similarity matrix 
% Input: .mat files containing timecourses for either the volume
% ROI's (Seitzman300) or the surface parcels (Gordon333)
% Output: a similarity matrix figure that will be saved, and the average
% between- and within-subject correlations (will not be saved)
% Written to work with both Lifespan and iNetworks subjects assumming 5 and
% 4 sessions for each dataset, respectively
% This script will force the same amount of data across the two datasets

clear all
%% ------------------------------------------------------------------------------------------------
%% PATHS
data_dir = '/Users/dianaperez/Desktop/';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

%% OPTIONS
datasets = {'Lifespan-NU', 'Lifespan-FSU', 'iNet-NU', 'iNet-FSU'};% 
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 1361; % if this is commented out or set to 0, then the script will calculate it
exclude_subs = {'LS46', 'INET108'};

%% LOAD SIMILARITY DATA
for d = 1:numel(datasets)
    disp(['dataset: ' datasets{d}])
    [subject, sessions, N] = get_subjects(datasets{d}, exclude_subs);
    % create a matrix mask 
    file_name = sprintf('%s/%s_similarityMat_matchedData.mat', data_dir, datasets{d});
    load(file_name);

    sub_labels = {};
    for s = 1:length(subject)
        sub_labels{s}=sprintf('sub-%d', s);
    end


%% MAKE FIGURE
figure('Position',[1 1 1000 800]);
load better_jet_colormap.mat
imagesc(simmat,[0 1]); colormap(better_jet_colormap_diff);%colormap("jet");%
tick_marks = [0:sessions:(5*numel(subject))]+0.5;
hline_new(tick_marks,'k',2);
vline_new(tick_marks,'k',2);
set(gca,'XTick',tick_marks(1:numel(subject))+(sessions/2), 'YTick', tick_marks(1:numel(subject))+(sessions/2), 'XTickLabel',...
    sub_labels, 'YTickLabel', sub_labels);
axis square;
colorbar;

set(gca,'FontSize',16)

title('Correlation Matrix Similarity');
if match_data
    saveas(gcf,[output_dir datasets{d} '_similarityMat_matchedData_' num2str(amt_data) '.jpg'],'jpg');
else
    saveas(gcf,[output_dir datasets{d} '_similarityMat_unMatchedData.jpg'],'jpg');
end
close('all');

%% CALCULATE average within- and between-subject correlations
count = 1;
within = [];
between = [];
sessions = ses_lims;
for s = 1:numel(subject)
    if any(contains(exclude_subs, subject{s}))
        continue;
    else
        lines = [count:(count+sessions(s)-1)];
        sub_vals = simmat(lines,:);
        maskmat = ones(sessions(s),sessions(s));
        maskmat = logical(triu(maskmat, 1));
        within_sub = sub_vals(:,lines);
        within = [within; within_sub(maskmat)];
        sub_means.within(s) = mean(within_sub(maskmat));
        maskmat = ones(size(sub_vals));
        maskmat(:,lines) = 0;
        between = [between; sub_vals(maskmat==1)];
        count = count+sessions(s);
        sub_means.between(s) = mean(sub_vals(maskmat==1));
        
    end
end
mean_between = mean(between);
mean_within = mean(within);
disp(['The average similarity between subjects is ' num2str(mean(between))])
disp(['The average similarity within subjects is ' num2str(mean(within))])
save([output_dir datasets{d} '_similarityMat_matchedData.mat'], 'simmat', 'amt_data', 'mean_between', 'mean_within', 'sub_means', 'subject', 'ses_lims')

% if ~strcmpi(datasets{d}, 'Lifespan-NU')
%     subsample = datasample(1:length(subject), 10, 'Replace', false);
%     sub_labels_subsample = sub_labels{subsample};
%     subject_subsample = subject{subsample};
%     simmat_subsample = simmat(subsample, subsample);
%     
%     %plot
%     figure('Position',[1 1 1000 800]);
%     imagesc(simmat_subsample,[0 1]); colormap(better_jet_colormap_diff);%colormap("jet");%
%     tick_marks = [0:sessions:(5*numel(subject_subsample))]+0.5;
%     hline_new(tick_marks,'k',2);
%     vline_new(tick_marks,'k',2);
%     set(gca,'XTick',tick_marks(1:numel(subject_subsample))+(sessions/2), 'YTick', tick_marks(1:numel(subject_subsample))+(sessions/2), 'XTickLabel',...
%         sub_labels_subsample, 'YTickLabel', sub_labels_subsample);
%     axis square;
%     colorbar;
% 
%     set(gca,'FontSize',16)
% 
%     title('Correlation Matrix Similarity');
%     if match_data
%         saveas(gcf,[output_dir datasets{d} '_similarityMat_matchedData_' num2str(amt_data) '_subsample.jpg'],'jpg');
%     else
%         saveas(gcf,[output_dir datasets{d} '_similarityMat_unMatchedData_subsample.jpg'],'jpg');
%     end
%     close('all');
% end
end
%% THE END

function [subject, sessions, N] = get_subjects(dataset, exclude_subs)

if strcmpi(dataset, 'Lifespan-NU')
    subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'}; %
    sessions = 5;
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
subject(:,find(contains(subject, exclude_subs))) = [];
N = numel(subject);

end