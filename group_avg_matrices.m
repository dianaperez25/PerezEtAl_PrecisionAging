%% MAKE GROUP AVG MATRICES


clear all
%% ------------------------------------------------------------------------------------------------
%% PATHS
%% ------------------------------------------------------------------------------------------------

data_dir = '/Users/dianaperez/Desktop/FC_Parcels_333/';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

%% ------------------------------------------------------------------------------------------------
%% OPTIONS
%% ------------------------------------------------------------------------------------------------

datasets = {'Lifespan-NU', 'Lifespan-FSU', 'iNet-FSU', 'iNet-NU', };% 
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 1361; % if this is commented out or set to 0, then the script will calculate it
exclude_subs = {'LS46', 'INET108'};

%% ------------------------------------------------------------------------------------------------
%% DATA MATCHING
%% ------------------------------------------------------------------------------------------------

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

%% ------------------------------------------------------------------------------------------------
%% AVERAGING MATRICES
%% ------------------------------------------------------------------------------------------------

all_datasets = [];
maskmat = ones(333);
maskmat = logical(triu(maskmat, 1));

atlas = 'Parcels333';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 

for d = 1:numel(datasets)
    disp(['dataset: ' datasets{d}])
    
    % sets number of subjects and sessions depending on the dataset being
    [subject, sessions, N] = get_subjects(datasets{d}, exclude_subs);

    % main loop; starts analysis
    dataset_mat = [];
    
    for s = 1:numel(subject)
            
        sub_mat = [];             
            
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

                    %disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(matched_data,1), size(matched_data,2)))

                    % ... calculate the correlation matrix...
                    corrmat_matched_data = paircorr_mod(matched_data');
                    %sorted_data = corrmat_matched_data(atlas_params.sorti,atlas_params.sorti);
                    % ... make it linear and store it in a variable...
                    session_mat = single(FisherTransform(corrmat_matched_data));
                    sub_mat(i,:,:,:) = session_mat;
                    %figure_corrmat_GrattonLab(session_mat, atlas_params,-1,1);

                    % ... then onto the next session.
                else
                    disp(sprintf('sub-%s ses-%d data file is missing!', subject{s}, i))
                end
        end
        sub_mean = squeeze(mean(sub_mat,1));
        dataset_mat(s,:,:) = sub_mean;
        %figure_corrmat_GrattonLab(sub_mean, atlas_params,-1,1);
        %saveas(gcf,sprintf('%s/%s_sessionAvg_parcel_corrmat.jpg', output_dir, subject{s}),'jpg'); close all
    end
    dataset_mean = squeeze(mean(dataset_mat,1));
    all_datasets(d,:,:,:) =  dataset_mean;
    figure_corrmat_GrattonLab(dataset_mean, atlas_params,-1,1);        
    saveas(gcf,sprintf('%s/%s_groupAvg_parcel_corrmat.jpg',output_dir, datasets{d}), 'jpg');
    all_datasets_corrlin(d,:) = single(FisherTransform(dataset_mean(maskmat))); close all
end

% then, calculate the correlation/similarity across all of those linear matrices
simmat = corr(all_datasets_corrlin');

%% MAKE FIGURE
figure('Position',[1 1 1000 800]);
heatmap(simmat, 'XData', ["Lifespan-NU", "iNET-NU", "Lifespan-FSU", "iNET-FSU"],'YData', ["Lifespan-NU", "iNET-NU", "Lifespan-FSU", "iNET-FSU"]); 
load better_jet_colormap.mat
colormap(better_jet_colormap_diff);
ax = gca;
ax.FontSize = 20;
if match_data
    saveas(gcf,[output_dir 'groupAvg_similarityMat_matchedData_' num2str(amt_data) '.jpg'],'jpg');
else
    saveas(gcf,[output_dir 'groupAvg_similarityMat_unMatchedData.jpg'],'jpg');
end
close('all');


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
    
    
    
    
   


% make diff mats
LSNU=squeeze(all_datasets(1,:,:,:));
LSFSU=squeeze(all_datasets(2,:,:,:));
iNetNU=squeeze(all_datasets(3,:,:,:));
iNetFSU=squeeze(all_datasets(4,:,:,:));

LSdiff = LSFSU - LSNU;
figure_corrmat_GrattonLab(LSdiff, atlas_params,-1,1)
clim([-.2 .2])
saveas(gcf,sprintf('%s/LSDiff_groupAvg_parcel_corrmat.jpg',output_dir), 'jpg')

iNetdiff = iNetFSU - iNetNU;
figure_corrmat_GrattonLab(iNetdiff, atlas_params,-1,1)
clim([-.2 .2])
saveas(gcf,sprintf('%s/iNetDiff_groupAvg_parcel_corrmat.jpg',output_dir), 'jpg')

NUdiff = iNetNU - LSNU;
figure_corrmat_GrattonLab(NUdiff, atlas_params,-1,1)
clim([-.2 .2])
saveas(gcf,sprintf('%s/NUDiff_groupAvg_parcel_corrmat.jpg',output_dir), 'jpg')

FSUdiff = iNetFSU - LSFSU;
figure_corrmat_GrattonLab(FSUdiff, atlas_params,-1,1)
clim([-.2 .2])
saveas(gcf,sprintf('%s/FSUDiff_groupAvg_parcel_corrmat.jpg',output_dir), 'jpg')
