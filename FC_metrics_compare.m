%% analyze FC metrics across young and older adults
% networks:
% 1- DMN, 2 - visual, 3 - FP, 5 - DAN
% 7 - VAN/lang, 8 - salience, 9 - CON, 10 - SMd
% 11 - SMl, 12 - auditory, 13 - Tpole, 14 - MTL
% 15 - PMN, 16 - PON
% cognitive control: FP, DAN, salience, CON, VAN/lang?
% somatomotor: visual, SMd, SMl, auditory
% memory/default: DMN, Tpole, MTL, PMN, PON?
% Chan et al definitions:
% sensorimotor: Visual, SMl, SMd, Auditory (indices: 2, 8, 9, 10)
% Association: DMN, FP, VAN, CON, DAN, Salience (indices: 1, 3, 4, 5, 6, 7)


FC_metrics_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FCmetrics';
dataset = {'Lifespan-NU', 'iNet-FSU', 'Lifespan-FSU', 'iNet-NU'};
exclude_subs = {'LS46', 'INET108', 'LS108', 'INET057'};


for d = 1:numel(dataset)
    [subjects, sessions, N] = get_subjects(dataset{d}, exclude_subs);

    for s = 1:numel(subjects)
        load(sprintf('%s/sub-%s_FCmetrics.mat', FC_metrics_dir, subjects{s}));
        %out_data = net_data;
        %% calculate weights for each network
        total_surf_indiv = sum(out_data.indiv(1,:));
        total_surf_group = sum(out_data.group(1,:));
        weights_indiv = out_data.indiv(1,:)./total_surf_indiv;
        weights_group = out_data.group(1,:)./total_surf_group;
        tmp = []; tmp_SM = []; tmp_assoc = [];
        for n = 1:14
        %% calculate wFC for all networks combined
            tmp(1,n) = out_data.indiv(2,n)*weights_indiv(n);
            tmp(2,n) = out_data.group(2,n)*weights_group(n);
        %% calculate bFC for all networks combined
            tmp(3,n) = out_data.indiv(3,n)*weights_indiv(n);
            tmp(4,n) = out_data.group(3,n)*weights_group(n);
        %% calculate segregation index for all networks combined
            tmp(5,n) = out_data.indiv(4,n)*weights_indiv(n);
            tmp(6,n) = out_data.group(4,n)*weights_group(n);
            if n == 2 || n == 8 || n == 9 || n == 10
                tmp_SM(1,end+1) = out_data.indiv(2,n);
                tmp_SM(2,end) = out_data.group(2,n);
                tmp_SM(3,end) = out_data.indiv(3,n);
                tmp_SM(4,end) = out_data.group(3,n);
                tmp_SM(5,end) = out_data.indiv(4,n);
                tmp_SM(6,end) = out_data.group(4,n);
            elseif n == 1 || n == 3 || n == 4 || n == 5 || n ==6 || n ==7
                tmp_assoc(1,end+1) = out_data.indiv(2,n);
                tmp_assoc(2,end) = out_data.group(2,n);
                tmp_assoc(3,end) = out_data.indiv(3,n);
                tmp_assoc(4,end) = out_data.group(3,n);
                tmp_assoc(5,end) = out_data.indiv(4,n);
                tmp_assoc(6,end) = out_data.group(4,n);
            end
        end

        
        if strcmpi(dataset{d}, 'iNet-NU')
            iNet_NU.dice(s,:) = out_data.dice_overlap;
            iNet_NU.size(s,:) = out_data.indiv(1,:);
            iNet_NU.size_diff(s,:) = out_data.indiv(1,:) - out_data.group(1,:);
            iNet_NU.allnets.summary.wFC(s,1) = sum(tmp(1,:)); % individualized nets
            iNet_NU.allnets.summary.wFC(s,2) = sum(tmp(2,:)); % group nets
            iNet_NU.allnets.bynet.wFC.indiv(s,:) = out_data.indiv(2,:);
            iNet_NU.allnets.bynet.wFC.group(s,:) = out_data.group(2,:);
            iNet_NU.allnets.bynet.bFC.indiv(s,:) = out_data.indiv(3,:);
            iNet_NU.allnets.bynet.bFC.group(s,:) = out_data.group(3,:);
            iNet_NU.allnets.summary.bFC(s,1) = sum(tmp(3,:)); % individualized nets
            iNet_NU.allnets.summary.bFC(s,2) = sum(tmp(4,:)); % group nets 
            iNet_NU.allnets.summary.seg(s,1) = (iNet_NU.allnets.summary.wFC(s,1) - iNet_NU.allnets.summary.bFC(s,1))/iNet_NU.allnets.summary.wFC(s,1); % individualized nets
            iNet_NU.allnets.summary.seg(s,2) = (iNet_NU.allnets.summary.wFC(s,2) - iNet_NU.allnets.summary.bFC(s,2))/iNet_NU.allnets.summary.wFC(s,2); % group nets
            iNet_NU.SM.summary.wFC(s,1) = mean(tmp_SM(1,:)); iNet_NU.assoc.wFC(s,1) = mean(tmp_assoc(1,:)); % individualized nets
            iNet_NU.SM.wFC(s,2) = mean(tmp_SM(2,:)); iNet_NU.assoc.summary.wFC(s,2) = mean(tmp_assoc(2,:)); % group nets
            iNet_NU.SM.summary.bFC(s,1) = mean(tmp_SM(3,:)); iNet_NU.assoc.summary.bFC(s,1) = mean(tmp_assoc(3,:)); % individualized nets
            iNet_NU.SM.summary.bFC(s,2) = mean(tmp_SM(4,:)); iNet_NU.assoc.summary.bFC(s,2) = mean(tmp_assoc(4,:)); % group nets
            iNet_NU.SM.summary.seg(s,1) = mean(tmp_SM(5,:)); iNet_NU.assoc.summary.seg(s,1) = mean(tmp_assoc(5,:)); % individualized nets
            iNet_NU.SM.summary.seg(s,2) = mean(tmp_SM(6,:)); iNet_NU.assoc.summary.seg(s,2) = mean(tmp_assoc(6,:)); % group nets
        elseif strcmpi(dataset{d}, 'Lifespan-NU')
            Lifespan_NU.dice(s,:) = out_data.dice_overlap;
            Lifespan_NU.size(s,:) = out_data.indiv(1,:);
            Lifespan_NU.size_diff(s,:) = out_data.indiv(1,:) - out_data.group(1,:);
            Lifespan_NU.allnets.summary.wFC(s,1) = sum(tmp(1,:));
            Lifespan_NU.allnets.summary.wFC(s,2) = sum(tmp(2,:));
            Lifespan_NU.allnets.summary.bFC(s,1) = sum(tmp(3,:));
            Lifespan_NU.allnets.summary.bFC(s,2) = sum(tmp(4,:));
            Lifespan_NU.allnets.summary.seg(s,1) = (Lifespan_NU.allnets.summary.wFC(s,1) - Lifespan_NU.allnets.summary.bFC(s,1))/Lifespan_NU.allnets.summary.wFC(s,1);
            Lifespan_NU.allnets.summary.seg(s,2) = (Lifespan_NU.allnets.summary.wFC(s,2) - Lifespan_NU.allnets.summary.bFC(s,2))/Lifespan_NU.allnets.summary.wFC(s,2);
            Lifespan_NU.allnets.bynet.wFC.indiv(s,:) = out_data.indiv(2,:);
            Lifespan_NU.allnets.bynet.wFC.group(s,:) = out_data.group(2,:);
            Lifespan_NU.allnets.bynet.bFC.indiv(s,:) = out_data.indiv(3,:);
            Lifespan_NU.allnets.bynet.bFC.group(s,:) = out_data.group(3,:);
            Lifespan_NU.SM.summary.wFC(s,1) = mean(tmp_SM(1,:)); Lifespan_NU.assoc.summary.wFC(s,1) = mean(tmp_assoc(1,:));
            Lifespan_NU.SM.summary.wFC(s,2) = mean(tmp_SM(2,:)); Lifespan_NU.assoc.summary.wFC(s,2) = mean(tmp_assoc(2,:));
            Lifespan_NU.SM.summary.bFC(s,1) = mean(tmp_SM(3,:)); Lifespan_NU.assoc.summary.bFC(s,1) = mean(tmp_assoc(3,:));
            Lifespan_NU.SM.summary.bFC(s,2) = mean(tmp_SM(4,:)); Lifespan_NU.assoc.summary.bFC(s,2) = mean(tmp_assoc(4,:));
            Lifespan_NU.SM.summary.seg(s,1) = mean(tmp_SM(5,:)); Lifespan_NU.assoc.summary.seg(s,1) = mean(tmp_assoc(5,:));
            Lifespan_NU.SM.summary.seg(s,2) = mean(tmp_SM(6,:)); Lifespan_NU.assoc.summary.seg(s,2) = mean(tmp_assoc(6,:));
        elseif contains(dataset{d}, 'iNet-FSU')
            iNet_FSU.dice(s,:) = out_data.dice_overlap;
            iNet_FSU.size(s,:) = out_data.indiv(1,:);
            iNet_FSU.size_diff(s,:) = out_data.indiv(1,:) - out_data.group(1,:);
            iNet_FSU.allnets.summary.wFC(s,1) = sum(tmp(1,:)); % individualized nets
            iNet_FSU.allnets.summary.wFC(s,2) = sum(tmp(2,:)); % group nets
            iNet_FSU.allnets.summary.bFC(s,1) = sum(tmp(3,:)); % individualized nets
            iNet_FSU.allnets.summary.bFC(s,2) = sum(tmp(4,:)); % group nets 
            iNet_FSU.allnets.summary.seg(s,1) = (iNet_FSU.allnets.summary.wFC(s,1) - iNet_FSU.allnets.summary.bFC(s,1))/iNet_FSU.allnets.summary.wFC(s,1); % individualized nets
            iNet_FSU.allnets.summary.seg(s,2) = (iNet_FSU.allnets.summary.wFC(s,2) - iNet_FSU.allnets.summary.bFC(s,2))/iNet_FSU.allnets.summary.wFC(s,2); % group nets
            iNet_FSU.allnets.bynet.wFC.indiv(s,:) = out_data.indiv(2,:);
            iNet_FSU.allnets.bynet.wFC.group(s,:) = out_data.group(2,:);
            iNet_FSU.allnets.bynet.bFC.indiv(s,:) = out_data.indiv(3,:);
            iNet_FSU.allnets.bynet.bFC.group(s,:) = out_data.group(3,:);
            iNet_FSU.SM.summary.wFC(s,1) = mean(tmp_SM(1,:)); iNet_FSU.assoc.summary.wFC(s,1) = mean(tmp_assoc(1,:)); % individualized nets
            iNet_FSU.SM.summary.wFC(s,2) = mean(tmp_SM(2,:)); iNet_FSU.assoc.summary.wFC(s,2) = mean(tmp_assoc(2,:)); % group nets
            iNet_FSU.SM.summary.bFC(s,1) = mean(tmp_SM(3,:)); iNet_FSU.assoc.summary.bFC(s,1) = mean(tmp_assoc(3,:)); % individualized nets
            iNet_FSU.SM.summary.bFC(s,2) = mean(tmp_SM(4,:)); iNet_FSU.assoc.summary.bFC(s,2) = mean(tmp_assoc(4,:)); % group nets
            iNet_FSU.SM.summary.seg(s,1) = mean(tmp_SM(5,:)); iNet_FSU.assoc.summary.seg(s,1) = mean(tmp_assoc(5,:)); % individualized nets
            iNet_FSU.SM.summary.seg(s,2) = mean(tmp_SM(6,:)); iNet_FSU.assoc.summary.seg(s,2) = mean(tmp_assoc(6,:)); % group nets
        elseif strcmpi(dataset{d}, 'Lifespan-FSU')
            Lifespan_FSU.dice(s,:) = out_data.dice_overlap;
            Lifespan_FSU.size(s,:) = out_data.indiv(1,:);
            Lifespan_FSU.size_diff(s,:) = out_data.indiv(1,:) - out_data.group(1,:);
            Lifespan_FSU.allnets.summary.wFC(s,1) = sum(tmp(1,:));
            Lifespan_FSU.allnets.summary.wFC(s,2) = sum(tmp(2,:));
            Lifespan_FSU.allnets.summary.bFC(s,1) = sum(tmp(3,:));
            Lifespan_FSU.allnets.summary.bFC(s,2) = sum(tmp(4,:));
            Lifespan_FSU.allnets.summary.seg(s,1) = (Lifespan_FSU.allnets.summary.wFC(s,1) - Lifespan_FSU.allnets.summary.bFC(s,1))/Lifespan_FSU.allnets.summary.wFC(s,1);
            Lifespan_FSU.allnets.summary.seg(s,2) = (Lifespan_FSU.allnets.summary.wFC(s,2) - Lifespan_FSU.allnets.summary.bFC(s,2))/Lifespan_FSU.allnets.summary.wFC(s,2);            
            Lifespan_FSU.allnets.bynet.wFC.indiv(s,:) = out_data.indiv(2,:);
            Lifespan_FSU.allnets.bynet.wFC.group(s,:) = out_data.group(2,:);
            Lifespan_FSU.allnets.bynet.bFC.indiv(s,:) = out_data.indiv(3,:);
            Lifespan_FSU.allnets.bynet.bFC.group(s,:) = out_data.group(3,:);
            Lifespan_FSU.SM.summary.wFC(s,1) = mean(tmp_SM(1,:)); Lifespan_FSU.assoc.summary.wFC(s,1) = mean(tmp_assoc(1,:));
            Lifespan_FSU.SM.summary.wFC(s,2) = mean(tmp_SM(2,:)); Lifespan_FSU.assoc.summary.wFC(s,2) = mean(tmp_assoc(2,:));
            Lifespan_FSU.SM.summary.bFC(s,1) = mean(tmp_SM(3,:)); Lifespan_FSU.assoc.summary.bFC(s,1) = mean(tmp_assoc(3,:));
            Lifespan_FSU.SM.summary.bFC(s,2) = mean(tmp_SM(4,:)); Lifespan_FSU.assoc.summary.bFC(s,2) = mean(tmp_assoc(4,:));
            Lifespan_FSU.SM.summary.seg(s,1) = mean(tmp_SM(5,:)); Lifespan_FSU.assoc.summary.seg(s,1) = mean(tmp_assoc(5,:));
            Lifespan_FSU.SM.summary.seg(s,2) = mean(tmp_SM(6,:)); Lifespan_FSU.assoc.summary.seg(s,2) = mean(tmp_assoc(6,:));
        end
        clear tmp
    end
end

%% t-test: comparing age group on difference in size between group and indidualized networks
[H,P,CI,STATS] = ttest2(mean(Lifespan_NU.size_diff,2), mean(iNet_NU.size_diff,2), 'Vartype', 'unequal')
[H,P,CI,STATS] = ttest2(mean(Lifespan_FSU.size_diff,2), mean(iNet_FSU.size_diff,2), 'Vartype', 'unequal')

%% t-test: comparing age group on difference in spatial correspondence between group and indidualized networks
[H,P,CI,STATS] = ttest2(mean(Lifespan_NU.dice,2), mean(iNet_NU.dice,2), 'Vartype', 'unequal')
[H,P,CI,STATS] = ttest2(mean(Lifespan_FSU.dice,2), mean(iNet_FSU.dice,2), 'Vartype', 'unequal')

%% 2x2 ANOVAS: age group(young adults & older adults) x data collection site (NU & FSU)
% set up labels for age_group factor
age_group = [];
age_group(1:8) = 1; % older adults NU
age_group(9:53) = 2; % younger adults NU
age_group(54:75) = 1; % older adults FSU
age_group(76:118) = 2; % younger adults FSU

site = [];
site(1:53) = 3; % data collection site --> NU
site(54:118) = 4; % data collection site --> FSU

factors = {age_group', site'};

% difference in how different the sizes are between individualized and
% group average networks
size_diff = [mean(Lifespan_NU.size_diff,2); mean(iNet_NU.size_diff,2); mean(Lifespan_FSU.size_diff,2); mean(iNet_FSU.size_diff,2)];
%aov_size_diff = anova(factors, size_diff, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site"])

% difference in dice coefficient between individualized and
% group average networks
dice = [mean(Lifespan_NU.dice,2); mean(iNet_NU.dice,2); mean(Lifespan_FSU.dice,2); mean(iNet_FSU.dice,2)];
%aov_dice = anova(factors, dice, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site"])


%% 2x2x2 ANOVAS: age group(young adults & older adults) x network scheme (group avg & individualized) x data collection site (NU & FSU)
% set up labels for age_group factor
% age_group = [];
% age_group(1:8) = 1; % older adults NU
% age_group(9:53) = 2; % younger adults NU
% age_group(54:75) = 1; % older adults FSU
% age_group(76:118) = 2; % younger adults FSU
age_group = [age_group'; age_group'];

% site = [];
% site(1:53) = 3; % data collection site --> NU
% site(54:118) = 4; % data collection site --> FSU

site = [site'; site'];

% set up labels for networks factor
networks = [];
networks(1:118) = 5; % individualized parcels
networks(119:236) = 6; % group parcels

factors = {age_group, site, networks'};

%% comparing within-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
wFC = [Lifespan_NU.allnets.wFC(:,1); iNet_NU.allnets.wFC(:,1); Lifespan_FSU.allnets.wFC(:,1); iNet_FSU.allnets.wFC(:,1); Lifespan_NU.allnets.wFC(:,2); iNet_NU.allnets.wFC(:,2); Lifespan_FSU.allnets.wFC(:,2); iNet_FSU.allnets.wFC(:,2)];
%aov_wFC = anova(factors, wFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site" "Network Parcellation"])

%% comparing between-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
bFC = [Lifespan_NU.allnets.bFC(:,1); iNet_NU.allnets.bFC(:,1); Lifespan_FSU.allnets.bFC(:,1); iNet_FSU.allnets.bFC(:,1); Lifespan_NU.allnets.bFC(:,2); iNet_NU.allnets.bFC(:,2); Lifespan_FSU.allnets.bFC(:,2); iNet_FSU.allnets.bFC(:,2)];
%aov_bFC = anova(factors, bFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site" "Network Parcellation"])

%% comparing segregation index
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
seg = [Lifespan_NU.allnets.seg(:,1); iNet_NU.allnets.seg(:,1); Lifespan_FSU.allnets.seg(:,1); iNet_FSU.allnets.seg(:,1); Lifespan_NU.allnets.seg(:,2); iNet_NU.allnets.seg(:,2); Lifespan_FSU.allnets.seg(:,2); iNet_FSU.allnets.seg(:,2)];
%aov_seg = anova(factors, seg, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site" "Network Parcellation"])

%% 2x2x2x2 ANOVAS: age group(young adults & older adults) x network parcellation (group avg & individualized) x network class (sensorimotor & association) x data collection site (NU & FSU)
% set up labels for age_group factor
age_group = [age_group; age_group];

% set up labels for site factor
site = [site; site];

% set up labels for network parcellation factor
networks = [networks'; networks'];

% set up labels for networks class factor
class = [];
class(1:236) = 7; % sensorimotor systems
class(237:472) = 8; % association systems
factors = {age_group, site, networks, class'};

%% comparing within-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
SM_wFC = [Lifespan_NU.SM.wFC(:,1); iNet_NU.SM.wFC(:,1); Lifespan_FSU.SM.wFC(:,1); iNet_FSU.SM.wFC(:,1); Lifespan_NU.SM.wFC(:,2); iNet_NU.SM.wFC(:,2); Lifespan_FSU.SM.wFC(:,2); iNet_FSU.SM.wFC(:,2)];
assoc_wFC = [Lifespan_NU.assoc.wFC(:,1); iNet_NU.assoc.wFC(:,1); Lifespan_FSU.assoc.wFC(:,1); iNet_FSU.assoc.wFC(:,1); Lifespan_NU.assoc.wFC(:,2); iNet_NU.assoc.wFC(:,2); Lifespan_FSU.assoc.wFC(:,2); iNet_FSU.assoc.wFC(:,2)];
wFC = [SM_wFC; assoc_wFC];
%aov_wFC = anova(factors, wFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site" "Network Parcellation" "Network Class"])

%% comparing between-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
SM_bFC = [Lifespan_NU.SM.bFC(:,1); iNet_NU.SM.bFC(:,1); Lifespan_FSU.SM.bFC(:,1); iNet_FSU.SM.bFC(:,1); Lifespan_NU.SM.bFC(:,2); iNet_NU.SM.bFC(:,2); Lifespan_FSU.SM.bFC(:,2); iNet_FSU.SM.bFC(:,2)];
assoc_bFC = [Lifespan_NU.assoc.bFC(:,1); iNet_NU.assoc.bFC(:,1); Lifespan_FSU.assoc.bFC(:,1); iNet_FSU.assoc.bFC(:,1); Lifespan_NU.assoc.bFC(:,2); iNet_NU.assoc.bFC(:,2); Lifespan_FSU.assoc.bFC(:,2); iNet_FSU.assoc.bFC(:,2)];
bFC = [SM_bFC; assoc_bFC];
%aov_bFC = anova(factors, bFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site" "Network Parcellation" "Network Class"])

%% comparing segregation index
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
SM_seg = [Lifespan_NU.assoc.seg(:,1); iNet_NU.assoc.seg(:,1); Lifespan_FSU.assoc.seg(:,1); iNet_FSU.assoc.seg(:,1); Lifespan_NU.assoc.seg(:,2); iNet_NU.assoc.seg(:,2); Lifespan_FSU.assoc.seg(:,2); iNet_FSU.assoc.seg(:,2)];
assoc_seg = [Lifespan_NU.assoc.seg(:,1); iNet_NU.assoc.seg(:,1); Lifespan_FSU.assoc.seg(:,1); iNet_FSU.assoc.seg(:,1); Lifespan_NU.assoc.seg(:,2); iNet_NU.assoc.seg(:,2); Lifespan_FSU.assoc.seg(:,2); iNet_FSU.assoc.seg(:,2)];
seg = [SM_seg; assoc_seg];
%aov_seg = anova(factors, seg, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Site" "Network Parcellation" "Network Class"])



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