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

%clear all
%% ------------------------------------------------------------------------------------------------
%% PATHS
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/FC_Parcels_333/';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

%% OPTIONS
datasets = {'Lifespan-NU'};% 
exclude_subs = {'LS46', 'INET108'};

% How much data to add at each step in minutes (how many data points that
% corresponds to is calculated below)
step_time = 2.5;
% TR of the fMRI images
TR = 1.1;
% How many points to sample for steps
step_datapts = round((60*step_time)/TR);% 136 -? ~2.5 min, 272 -> ~5 min, will add this number of frames each time it subsamples data


%% SIMILARITY CALCULATIONS
for d = 1:numel(datasets)
    disp(['dataset: ' datasets{d}])
    % create a matrix mask 
    maskmat = ones(333);
    maskmat = logical(triu(maskmat, 1));
    % sets number of subjects and sessions depending on the dataset being
    [subject, sessions, N] = get_subjects(datasets{d}, exclude_subs);

    % main loop; starts analysis
    count = 1; matcheddata_corrlin = []; %sub_data = [];
%     for s = 1:numel(subject)
%         
%             ses_count = 0;
%             for i = 1:sessions
% 
%                 % for each session and each subject, load the timeseries data...
%                 data_file = [data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'];
%                 if exist(data_file, 'file')
%                     load(data_file)
%                     sub_data{s,i} = parcel_time(logical(tmask_concat), :)';
%                 else
%                     disp(sprintf('sub-%s ses-%d data file is missing!', subject{s}, i))
%                 end
%                 end
%             end

    num_data = 0; step = 1; 
    for p = 1:12
        ses_lims = 0; count = 1; matcheddata_corrlin =[];
        this_step_num_data = (step_datapts*step);

        for s = 1:numel(subject)
            ses_count = 0;
            for i = 1:5
                this_sub_data = sub_data{s,i};
                if this_step_num_data > size(this_sub_data,2)
                    disp(sprintf('subject %s ses %d is NOT included in step %d...', subject{s}, i, step))
                    continue;
                else
                    
                    matched_data = this_sub_data(:,1:this_step_num_data);
                    % ... calculate the correlation matrix...
                    corrmat_matched_data = paircorr_mod(matched_data');

                    % ... make it linear and store it in a variable...
                    matcheddata_corrlin(count,:) = single(FisherTransform(corrmat_matched_data(maskmat)));

                    % ... then onto the next session.
                    count = count + 1; ses_count = ses_count + 1;
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
%tick_marks = [0:sessions:40]+0.5;
%hline_new(tick_marks,'k',2);
%vline_new(tick_marks,'k',2);
axis square;
colorbar;
title('Correlation Matrix Similarity');

    saveas(gcf,[output_dir datasets{d} '_similarityMat_matchedData_' num2str(this_step_num_data) '.jpg'],'jpg');

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
        lines = [count:(count+sessions(s+1)-1)];
        sub_vals = simmat(lines,:);
        maskmat2 = ones(sessions(s+1),sessions(s+1));
        maskmat2 = logical(triu(maskmat2, 1));
        within_sub = sub_vals(:,lines);
        within = [within; within_sub(maskmat2)];
        sub_mean_within(s,p) = mean(within_sub(maskmat2));
        maskmat2 = ones(size(sub_vals));
        maskmat2(:,lines) = 0;
        between = [between; sub_vals(maskmat2==1)];
        sub_mean_between(s,p) = mean(sub_vals(maskmat2==1));
        count = count+sessions(s+1);
    end
end
mean_between = mean(between);
mean_within = mean(within);
disp(['The average similarity between subjects is ' num2str(mean(between))])
disp(['The average similarity within subjects is ' num2str(mean(within))])
save([output_dir datasets{d} '_similarityMat_matchedData_' num2str(this_step_num_data) '.mat'], 'simmat', 'this_step_num_data', 'mean_between', 'mean_within', 'subject')
%step_between(s,p)= sub_mean_between;
%step_within(s,p) = sub_mean_within;
step = step+1;
    end
    
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
    
    
    
    
   



