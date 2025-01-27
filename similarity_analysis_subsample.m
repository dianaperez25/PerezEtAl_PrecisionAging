%% Script for Similarity analysis, but with option to force the same amount of data per subject 
% This script calculates the correlation (i.e. the similarity) of rs-FC matrices
% across sessions and subjects. 
% Input: .mat files containing timecourses for either the volume
% ROI's (Seitzman300) or the surface parcels (Gordon333)
% Output: a similarity matrix figure that will be saved, and the average
% between- and within-subject correlations (will not be saved)
% Written to work with both Lifespan and iNetworks subjects assumming 5 and
% 4 sessions for each dataset, respectively
% This script will force the same amount of data across the two datasets
%% THIS SCRIPT ONLY LOOKS AT SUBSAMPLE OF SUBJECTS FOR ILLUSTRATION PURPOSES

clear all
%% ------------------------------------------------------------------------------------------------
%% PATHS
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/Nifti/FC_Parcels_333/';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

%% OPTIONS
datasets = {'Lifespan-FSU', 'iNet-NU', 'iNet-FSU'};% 'Lifespan-NU', 
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 1361; % if this is commented out or set to 0, then the script will calculate it
exclude_subs = {'LS46', 'INET108', 'LS72'};
sample_size = 10;

%% DATA MATCHING
if match_data && amt_data == 0
    allSubs_amtData = [];
    
    for d = 1:numel(datasets)
        count = 1;
        [subject, sessions, N] = get_subjects(datasets{d}, exclude_subs);            
        % get minimum amt of data per session
        for s = 1:N           
                for i = 1:sessions
                    allSubs_amtData(count, 1) = d;
                    allSubs_amtData(count, 2) = s;
                    allSubs_amtData(count, 3) = i;
                    data_file = [data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'];
                    if exist(data_file, 'file')
                        load([data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
                        masked_data = parcel_time(logical(tmask_concat),:)';
                        allSubs_amtData(count, 4) = size(masked_data,2);
                        count = count + 1;
                    else
                        disp(sprintf('sub-%s ses-%d file not found!', subject{s}, i))
                    end
                end
            end
    end
    amt_data = min(allSubs_amtData(:,4));
end

%% SIMILARITY CALCULATIONS
for d = 1:numel(datasets)
    disp(['dataset: ' datasets{d}])
    % create a matrix mask 
    maskmat = ones(333);
    maskmat = logical(triu(maskmat, 1));
    % sets number of subjects and sessions depending on the dataset being
    [subject, sessions, N] = get_subjects(datasets{d}, exclude_subs);
    subsample = datasample(subject, 10, 'Replace', false);
    subject = subsample;
    % main loop; starts analysis
    count = 1; ses_lims = []; matcheddata_corrlin = [];
    for s = 1:numel(subject)
        
            ses_count = 0;
            for i = 1:sessions

                % for each session and each subject, load the timeseries data...
                data_file = [data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'];
                if exist(data_file, 'file')
                    load(data_file)
                    masked_data = parcel_time(logical(tmask_concat), :)';

                    % ... then sample the pre-defined amount of data from the timeseries data...
                    if size(masked_data,2)<amt_data || match_data == 0
                       matched_data = masked_data;
                    else
                        matched_data = masked_data(:,1:amt_data);
                    end

                    disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(matched_data,1), size(matched_data,2)))

                    % ... calculate the correlation matrix...
                    corrmat_matched_data = paircorr_mod(matched_data');

                    % ... make it linear and store it in a variable...
                    matcheddata_corrlin(count,:) = single(FisherTransform(corrmat_matched_data(maskmat)));

                    % ... then onto the next session.
                    count = count + 1; ses_count = ses_count + 1;
                else
                    disp(sprintf('sub-%s ses-%d data file is missing!', subject{s}, i))
                end
            end
            ses_lims = [ses_lims; ses_count];
    end

% then, calculate the correlation/similarity across all of those linear matrices
simmat = corr(matcheddata_corrlin');

%% MAKE FIGURE
figure('Position',[1 1 1000 800]);
load better_jet_colormap.mat
imagesc(simmat,[0 1]); colormap(better_jet_colormap_diff);
tick_marks = [0:sessions:(5*numel(subject))]+0.5;
hline_new(tick_marks,'k',2);
vline_new(tick_marks,'k',2);
set(gca,'XTick',tick_marks(1:numel(subject))+(sessions/2), 'YTick', tick_marks(1:numel(subject))+(sessions/2), 'XTickLabel',...
    subject, 'YTickLabel', subject);
axis square;
colorbar;
title('Correlation Matrix Similarity');
if match_data
    saveas(gcf,[output_dir datasets{d} '_similarityMat_matchedData_' num2str(amt_data) '_subsample.jpg'],'jpg');
else
    saveas(gcf,[output_dir datasets{d} '_similarityMat_unMatchedData_subsample.jpg'],'jpg');
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
save([output_dir datasets{d} '_similarityMat_matchedData_subsample.mat'], 'simmat', 'amt_data', 'mean_between', 'mean_within', 'sub_means', 'subject', 'ses_lims')

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
    
    
    
    
   



