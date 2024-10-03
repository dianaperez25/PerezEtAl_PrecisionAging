%% Figures for papers


data_dir = '/Users/dianaperez/Desktop/';

%% Average Number of Parcels
figure;
%load the file
fname = sprintf('%s/Lifespan-NU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = num_parcs;
type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])

fname = sprintf('%s/iNet-NU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = num_parcs;
type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])

fname = sprintf('%s/Lifespan-FSU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = num_parcs;
type = 'left';
pos = 2.5;
col = [0 0.4470 0.7410];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])

fname = sprintf('%s/iNet-FSU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = num_parcs;
type = 'right';
pos = 2.5;
col = [0.6350 0.0780 0.1840];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])


fname = sprintf('%s/GroupAvg_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

yline(num_parcs, '-', 'Group Avg', 'LineWidth',3)

%% Average Parcel Size in Number of Vertices
figure;
%load the file
fname = sprintf('%s/Lifespan-NU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,1);
type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])

fname = sprintf('%s/iNet-NU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,1);
type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])

fname = sprintf('%s/Lifespan-FSU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,1);
type = 'left';
pos = 2.5;
col = [0 0.4470 0.7410];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])

fname = sprintf('%s/iNet-FSU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,1);
type = 'right';
pos = 2.5;
col = [0.6350 0.0780 0.1840];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 0.5, col, type, [], [])


fname = sprintf('%s/GroupAvg_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

yline(avg_parc_size(1), '-', 'Group Avg', 'LineWidth',3)

%% Average Parcel Size in mm^2
figure;
%load the file
fname = sprintf('%s/Lifespan-NU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,2);
type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(x, pos, 1, col, type, [], [])

fname = sprintf('%s/iNet-NU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,2);
type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 1, col, type, [], [])

fname = sprintf('%s/Lifespan-FSU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,2);
type = 'left';
pos = 2.5;
col = [0 0.4470 0.7410];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 1, col, type, [], [])

fname = sprintf('%s/iNet-FSU_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

x = avg_parc_size(:,2);
type = 'right';
pos = 2.5;
col = [0.6350 0.0780 0.1840];

[h, mu, sigma, q, notch] = al_goodplot(x, pos, 1, col, type, [], [])


fname = sprintf('%s/GroupAvg_descriptiveIndivParcelInfo.mat', data_dir);
load(fname);

yline(avg_parc_size(2), '-', 'Group Avg', 'LineWidth',3)

%% COMPARING INDIV VS GROUP NETS
load('/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FC_Metrics/Lifespan-NU_group_vs_indiv_compare_size_spatial.mat')
Lifespan_NU = out_data;
load('/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FC_Metrics/Lifespan-FSU_group_vs_indiv_compare_size_spatial.mat')
Lifespan_FSU = out_data;
load('/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FC_Metrics/iNet-NU_group_vs_indiv_compare_size_spatial.mat')
iNet_NU = out_data;
load('/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FC_Metrics/iNet-FSU_group_vs_indiv_compare_size_spatial.mat')
iNet_FSU = out_data;
clear out_data

% size difference group - indiv in mm2 averaged across nets
Lifespan_NU_size_diff = mean((Lifespan_NU.group.size_mm2 - Lifespan_NU.indiv.size_mm2)/Lifespan_NU.group.size_mm2, 2);
iNet_NU_size_diff = mean((iNet_NU.group.size_mm2 - iNet_NU.indiv.size_mm2)/iNet_NU.group.size_mm2,2);
Lifespan_FSU_size_diff = mean((Lifespan_FSU.group.size_mm2 - Lifespan_FSU.indiv.size_mm2)/Lifespan_FSU.group.size_mm2,2);
iNet_FSU_size_diff = mean((iNet_FSU.group.size_mm2 - iNet_FSU.indiv.size_mm2)/iNet_FSU.group.size_mm2,2);

figure; 

type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_NU_size_diff(:), pos,.5, col, type, [], [])

type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_NU_size_diff(:), pos, .5, col, type, [], [])

figure;
type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_FSU_size_diff(:), pos, 1, col, type, [], [])

type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_FSU_size_diff(:), pos, 1, col, type, [], [])

%%

% size difference group - indiv in mm2 looking at each network 
Lifespan_NU_size_diff = (Lifespan_NU.group.size_mm2 - Lifespan_NU.indiv.size_mm2)./Lifespan_NU.group.size_mm2;
iNet_NU_size_diff = (iNet_NU.group.size_mm2 - iNet_NU.indiv.size_mm2)./iNet_NU.group.size_mm2;
Lifespan_FSU_size_diff = (Lifespan_FSU.group.size_mm2 - Lifespan_FSU.indiv.size_mm2)./Lifespan_FSU.group.size_mm2;
iNet_FSU_size_diff = (iNet_FSU.group.size_mm2 - iNet_FSU.indiv.size_mm2)./iNet_FSU.group.size_mm2;

Lifespan_NU_size_diff = (Lifespan_NU.group.size_mm2 - Lifespan_NU.indiv.size_mm2);
iNet_NU_size_diff = (iNet_NU.group.size_mm2 - iNet_NU.indiv.size_mm2);
Lifespan_FSU_size_diff = (Lifespan_FSU.group.size_mm2 - Lifespan_FSU.indiv.size_mm2);
iNet_FSU_size_diff = (iNet_FSU.group.size_mm2 - iNet_FSU.indiv.size_mm2);

figure; 

for net = 1:11
    type = 'left';
    pos = net;
    col = [0 0.4470 0.7410];
    [h, mu, sigma, q, notch] = al_goodplot(Lifespan_NU_size_diff(:,net), pos, .5, col, type, [], 1);

    type = 'right';
    pos = net;
    col = [0.6350 0.0780 0.1840];
    [h, mu, sigma, q, notch] = al_goodplot(iNet_NU_size_diff(:,net), pos, .5, col, type, [], 1);
end

figure;
type = 'left';
pos = 1:11;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_FSU_size_diff(:,1:11), pos, .5, col, type, [], 1.5)

type = 'right';
pos = 1:11;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_FSU_size_diff(:,1:11), pos, .5, col, type, [], 1.5)


%% DICE COEFFICIENT

Lifespan_NU_dice = mean(Lifespan_NU.dice_overlap,2);
iNet_NU_dice = mean(iNet_NU.dice_overlap,2);
Lifespan_FSU_dice = mean(Lifespan_FSU.dice_overlap,2);
iNet_FSU_dice = mean(iNet_FSU.dice_overlap,2);
figure; 

type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_NU_dice(:), pos, 1, col, type, [], [])

type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_NU_dice(:), pos, 1, col, type, [], [])

figure;
type = 'left';
pos = 1;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_FSU_dice(:), pos, 1, col, type, [], [])

type = 'right';
pos = 1;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_FSU_dice(:), pos, 1, col, type, [], [])

% by network

Lifespan_NU_dice = Lifespan_NU.dice_overlap;
iNet_NU_dice = iNet_NU.dice_overlap;
Lifespan_FSU_dice = Lifespan_FSU.dice_overlap;
iNet_FSU_dice = iNet_FSU.dice_overlap;
figure; 

type = 'left';
pos = 1:11;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_NU_dice(:,1:11), pos, 1, col, type, [], 2)

type = 'right';
pos = 1:11;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_NU_dice(:,1:11), pos, 1, col, type, [], 2)

figure; 
type = 'left';
pos = 1:11;
col = [0 0.4470 0.7410];
[h, mu, sigma, q, notch] = al_goodplot(Lifespan_FSU_dice(:,1:11), pos, 1, col, type, [], 1.5)

type = 'right';
pos = 1:11;
col = [0.6350 0.0780 0.1840];
[h, mu, sigma, q, notch] = al_goodplot(iNet_FSU_dice(:,1:11), pos, 1, col, type, [], 1.5)