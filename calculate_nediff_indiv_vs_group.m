clear all

addpath(genpath('/Users/dianaperez/Documents/Resources/'))

template_match_dir = ['/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/template_match/'];
group_atlas_fname = '/Users/dianaperez/Documents/Resources/WashU120_groupNetworks.dtseries.nii';
surf_area_fname = '/Users/dianaperez/Documents/Resources/surf_areas_verts.dtseries.nii';
out_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FC_Metrics/';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

datasets = {'Lifespan-NU', 'iNet-NU', 'Lifespan-FSU', 'iNet-FSU'};% 
atlas = 'Parcels333';
exclude_subs = {'INET001', 'INET061', 'LS47', 'LS46', 'INET108'};


% load the group average networks
group_avg = ft_read_cifti_mod(group_atlas_fname);
group_avg = group_avg.data(1:59412);
nets = unique(group_avg); nets(nets==0) = [];

% load the surface area file
surf_area = ft_read_cifti_mod(surf_area_fname);
surf_area = surf_area.data;

% for each dataset
for d = 1:numel(datasets)
    [subjects, num_ses, N] = get_subjects(datasets{d}, exclude_subs)

% for each subject
    for sub = 1:N
        out_data = [];

        indiv_net_fname = sprintf('%s/sub-%s_dice_WTA_map_kden0.05.dtseries.nii', template_match_dir, subjects{sub});
        indiv_nets = ft_read_cifti_mod(indiv_net_fname);
        indiv_nets = indiv_nets.data;

        for n = 1:length(nets)
            %index network vertices in individual
            indiv_net_verts = find(indiv_nets==nets(n));
            indiv_this_net = zeros(59412,1);
            indiv_this_net(indiv_net_verts) = 1;
            
            %index network vertices in group average
            group_verts = find(group_avg==nets(n));
            group_net = zeros(59412,1);
            group_net(group_verts) = 1;
            
            %calculate surface area of individual network
            out_data.indiv(2,n) = length(indiv_net_verts);
            out_data.indiv(1,n) = sum(surf_area(indiv_net_verts));
            
            %calculate surface area of individual network
            out_data.group(1,n) = length(group_verts);
            out_data.group(2,n) = sum(surf_area(group_verts));
            
            %calculate dice overlap bewteen individual and group average
            %nets
            out_data.dice_overlap(n) = dice_coefficient_mod(group_net, indiv_this_net);
        end
        outfile = sprintf('%s/sub-%s_FCmetrics.mat', out_dir, subjects{sub});
        save(outfile, 'out_data', '-v7.3')
    end
end



function [subject, sessions, N] = get_subjects(dataset, exclude_subs)

if strcmpi(dataset, 'Lifespan-NU')
    subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'}; %
    sessions = 5;
elseif strcmpi(dataset, 'iNet-NU')
    subject = {'INET002', 'INET003', 'INET005', 'INET006','INET010',...
        'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073',...
        'INET061', 'INET001', 'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
     'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
    'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
    'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
    'INET057', 'INET058', 'INET059', 'INET060', 'INET062', 'INET063',...
    'INET065', 'INET067'}; 
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
    'LS76', 'LS77', 'LS79', 'LS85', 'LS89', 'LS94', 'LS108'}; %
    sessions = 5;
end
subject(:,find(contains(subject, exclude_subs))) = [];
N = numel(subject);

end