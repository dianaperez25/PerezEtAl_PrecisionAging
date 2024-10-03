function calculate_FC_metrics(subject, indiv_network, dconn_fname, outdir)

% load surface area file
surf_area = ft_read_cifti_mod('/projects/p32293/needed_files/surf_areas_verts.dtseries.nii');
% load group average
group_avg = ft_read_cifti_mod('/projects/p32293/needed_files/WashU120_groupNetworks.dtseries.nii');
group_avg = group_avg.data(1:59412);
%initialize variables to store information
%wFC = []; bFC = []; seg = []; %within-network and between-network FC and segregation index
%surf_area = []; dice_overlap = [];    
out_data = [];
%netmap_fname = sprintf('/scratch/dcr8536/postFCproc_CIFTI/template_matching/sub-%s_dice_WTA_map_kden0.05.dtseries.nii', subjects{sub});
net_map = ft_read_cifti_mod(indiv_network);
dconn = ft_read_cifti_mod(dconn_fname);
dconn = dconn.data(1:59412, 1:59412);
nets = unique(group_avg); nets(nets==0) = [];
%allwFC = []; allbFC = []; allwFC_2=[]; allbFC_2=[];
for n = 1:length(nets)
    %index network vertices in individual
    indiv_verts = find(net_map.data==nets(n));
    indiv_net = zeros(59412,1);
    indiv_net(indiv_verts) = 1;
    
    %index network vertices in group average
    group_verts = find(group_avg==nets(n));
    group_net = zeros(59412,1);
    group_net(group_verts) = 1;
    
    %calculate surface area of individual network
    out_data.indiv(1,n) = sum(surf_area.data(indiv_verts));
    
    %calculate surface area of individual network
    out_data.group(1,n) = sum(surf_area.data(group_verts));
    
    %calculate dice overlap bewteen individual FPN and group average FPN
    out_data.dice_overlap(n) = dice_coefficient_mod(group_net, indiv_net);

    %% individual network
    %index network functional connectivity
    indiv_fc = dconn(indiv_verts,:); %individual FPN
    
    %calculate within-network connectivity
    tmp_wFC = indiv_fc(:,indiv_verts);
    %allwFC = [allwFC; tmp_wFC(:)];
    indiv_wFC = mean(mean(tmp_wFC(tmp_wFC>0), 'omitnan'));
    out_data.indiv(2,n) = indiv_wFC;
    
    %calculate between-network connectivity
    net_mask = ones(size(indiv_fc));
    net_mask(:,indiv_verts) = 0;
    tmp_bFC = indiv_fc(logical(net_mask));
    %allbFC = [allbFC; tmp_bFC(:)];
    indiv_bFC = mean(tmp_bFC(tmp_bFC>0), 'omitnan');
    out_data.indiv(3,n) = indiv_bFC;
    clear net_mask indiv_fc tmp_wFC tmp_bFC
    
    %calculate network segregation index
    indiv_seg = (indiv_wFC - indiv_bFC)/indiv_wFC;
    out_data.indiv(4,n) = indiv_seg;
    clear indiv_seg indiv_wFC indiv_bFC
    
    %% group network
    group_fc = dconn(group_verts,:); %group network
    
    %calculate within-network connectivity
    tmp_wFC = group_fc(:,group_verts);
    %allwFC_2 = [allwFC_2; tmp_wFC(:)];
    group_wFC = mean(mean(tmp_wFC(tmp_wFC>0), 'omitnan'));
    out_data.group(2,n) = group_wFC;

    %calculate between-network connectivity
    net_mask = ones(size(group_fc));
    net_mask(:,group_verts) = 0;
    tmp_bFC = group_fc(logical(net_mask));
    %allbFC_2 = [allbFC_2; tmp_bFC(:)];
    group_bFC = mean(tmp_bFC(tmp_bFC>0), 'omitnan');
    out_data.group(3,n) = group_bFC;
    clear net_mask group_fc tmp_wFC tmp_bFC

    %calculate network segregation index
    group_seg = (group_wFC - group_bFC)/group_wFC;
    out_data.group(4,n) = group_seg;
    clear group_seg group_wFC group_bFC    
end
%out_data.indiv(2,end+1)=mean(allwFC(allwFC>0), 'omitnan');
%out_data.group(2,end+1)=mean(allwFC_2(allwFC_2>0), 'omitnan');
%out_data.indiv(3,end)=mean(allbFC(allbFC>0), 'omitnan');
%out_data.group(3,end)=mean(allbFC_2(allbFC_2>0), 'omitnan');
%out_data.indiv(4,end)=(mean(allwFC(allwFC>0), 'omitnan')-mean(allbFC(allbFC>0), 'omitnan'))/(mean(allwFC(allwFC>0), 'omitnan'));
%out_data.group(4,end)=(mean(allwFC_2(allwFC_2>0), 'omitnan')-mean(allbFC_2(allbFC_2>0), 'omitnan'))/(mean(allwFC_2(allwFC_2>0), 'omitnan'));
outfile = sprintf('%s/sub-%s_FCmetrics.mat', outdir, subject);
save(outfile, 'out_data', '-v7.3')
end