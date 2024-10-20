%% Make overlap maps for OA's and for YA's of individualized networks for the networks that differed from Group Avg
% language network (7), FP (3), and SMd (10) differed in size (FP and lang smaller diff
% in OA and SMd bigger diff in OA). 
% visual(2), CO(9), SMl (11) and SMd (10), and Aud (12) differed in spatial arrangement: all
% more similar to typical distribution in YA
network = 12; 
out_str = 'Aud';
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/template_match/';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

datasets = {'Lifespan-NU', 'Lifespan-FSU', 'iNet-FSU', 'iNet-NU'};% 
exclude_subs = {'LS46', 'INET108'};

% load a template
template = ft_read_cifti_mod('/Users/dianaperez/Desktop/template.dtseries.nii');

for d = 1:numel(datasets)
    [subjects, sessions, N] = get_subjects(datasets{d}, exclude_subs);
    all_subs = [];
    for sub = 1:N
        % first load the individualized map for each subject
        indivnets_fname = sprintf('%s/sub-%s_dice_WTA_map_kden0.05.dtseries.nii', data_dir, subjects{sub});
        indiv_nets = ft_read_cifti_mod(indivnets_fname);

        %then binarize the network map to only have ones for the network we
        %want
        bin_netmap = zeros(size(indiv_nets.data));
        bin_netmap = logical(indiv_nets.data==network);
        all_subs(:,sub) = bin_netmap;
    end
    overlap_map = template;
    overlap_map.data = sum(all_subs,2)/N;
    out_fname = sprintf('%s/%s_%s_indivnet_overlap_map.dtseries.nii', output_dir, datasets{d}, out_str);
    ft_write_cifti_mod(out_fname, overlap_map);
end
        