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


FC_metrics_dir = '/projects/p32293/Lifespan/FCmetrics';
dataset = {'iNet', 'Lifespan'};

iNet_wFC = []; iNet_bFC = []; iNet_seg = []; iNet_dice = []; iNet_size_diff = []; iNet_SM = []; iNet_assoc = [];
Lifespan_wFC = []; Lifespan_bFC = []; Lifespan_seg = []; Lifespan_dice = []; Lifespan_size_diff = []; Lifespan_SM = []; Lifespan_assoc = [];

for d = 1:numel(dataset)
    if strcmp(dataset{d}, 'iNet')
         subjects = {'INET002', 'INET003', 'INET005', 'INET006','INET010',...
        'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
        'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
        'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
        'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
        'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
        'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'};%,...'INET001',
    elseif strcmp(dataset{d}, 'Lifespan')
        subjects={'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
    end

    for s = 1:numel(subjects)
        load(sprintf('%s/sub-%s_FCmetrics.mat', FC_metrics_dir, subjects{s}));
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
        %indiv_wFC(d,s) = sum(tmp(1,:));
        %group_wFC(d,s) = sum(tmp(2,:));
%         indiv_bFC(d,s) = sum(tmp(3,:));
%         group_bFC(d,s) = sum(tmp(4,:));
%         indiv_seg(d,s) = (indiv_wFC(s) - indiv_bFC(s))/indiv_wFC(s);
%         group_seg(d,s) = (group_wFC(s) - group_bFC(s))/group_wFC(s);
        
        if strcmp(dataset{d}, 'iNet')
            iNet_dice(s,:) = out_data.dice_overlap;
            iNet_size_diff(s,:) = out_data.indiv(1,:) - out_data.group(1,:);
            iNet_wFC(s,1) = sum(tmp(1,:)); % individualized nets
            iNet_wFC(s,2) = sum(tmp(2,:)); % group nets
            iNet_bFC(s,1) = sum(tmp(3,:)); % individualized nets
            iNet_bFC(s,2) = sum(tmp(4,:)); % group nets 
            iNet_seg(s,1) = (iNet_wFC(s,1) - iNet_bFC(s,1))/iNet_wFC(s,1); % individualized nets
            iNet_seg(s,2) = (iNet_wFC(s,2) - iNet_bFC(s,2))/iNet_wFC(s,2); % group nets
            iNet_SM.wFC(s,1) = mean(tmp_SM(1,:)); iNet_assoc.wFC(s,1) = mean(tmp_assoc(1,:)); % individualized nets
            iNet_SM.wFC(s,2) = mean(tmp_SM(2,:)); iNet_assoc.wFC(s,2) = mean(tmp_assoc(2,:)); % group nets
            iNet_SM.bFC(s,1) = mean(tmp_SM(3,:)); iNet_assoc.bFC(s,1) = mean(tmp_assoc(3,:)); % individualized nets
            iNet_SM.bFC(s,2) = mean(tmp_SM(4,:)); iNet_assoc.bFC(s,2) = mean(tmp_assoc(4,:)); % group nets
            iNet_SM.seg(s,1) = mean(tmp_SM(5,:)); iNet_assoc.seg(s,1) = mean(tmp_assoc(5,:)); % individualized nets
            iNet_SM.seg(s,2) = mean(tmp_SM(6,:)); iNet_assoc.seg(s,2) = mean(tmp_assoc(6,:)); % group nets
        elseif strcmp(dataset{d}, 'Lifespan')
            Lifespan_dice(s,:) = out_data.dice_overlap;
            Lifespan_size_diff(s,:) = out_data.indiv(1,:) - out_data.group(1,:);
            Lifespan_wFC(s,1) = sum(tmp(1,:));
            Lifespan_wFC(s,2) = sum(tmp(2,:));
            Lifespan_bFC(s,1) = sum(tmp(3,:));
            Lifespan_bFC(s,2) = sum(tmp(4,:));
            Lifespan_seg(s,1) = (Lifespan_wFC(s,1) - Lifespan_bFC(s,1))/Lifespan_wFC(s,1);
            Lifespan_seg(s,2) = (Lifespan_wFC(s,2) - Lifespan_bFC(s,2))/Lifespan_wFC(s,2);
            Lifespan_SM.wFC(s,1) = mean(tmp_SM(1,:)); Lifespan_assoc.wFC(s,1) = mean(tmp_assoc(1,:));
            Lifespan_SM.wFC(s,2) = mean(tmp_SM(2,:)); Lifespan_assoc.wFC(s,2) = mean(tmp_assoc(2,:));
            Lifespan_SM.bFC(s,1) = mean(tmp_SM(3,:)); Lifespan_assoc.bFC(s,1) = mean(tmp_assoc(3,:));
            Lifespan_SM.bFC(s,2) = mean(tmp_SM(4,:)); Lifespan_assoc.bFC(s,2) = mean(tmp_assoc(4,:));
            Lifespan_SM.seg(s,1) = mean(tmp_SM(5,:)); Lifespan_assoc.seg(s,1) = mean(tmp_assoc(5,:));
            Lifespan_SM.seg(s,2) = mean(tmp_SM(6,:)); Lifespan_assoc.seg(s,2) = mean(tmp_assoc(6,:));
        end
        clear tmp
    end
end

%% t-test: comparing age group on difference in size between group and indidualized networks
[H,P,CI,STATS] = ttest2(mean(Lifespan_size_diff,2), mean(iNet_size_diff,2), 'Vartype', 'unequal');
%% t-test: comparing age group on difference in spatial correspondence between group and indidualized networks
[H,P,CI,STATS] = ttest2(mean(Lifespan_dice,2), mean(iNet_dice,2), 'Vartype', 'unequal');

%% 2x2 ANOVAS: age group(young adults & older adults) x network scheme (group avg & individualized)
% set up labels for age_group factor
age_group = [];
age_group(1:8) = 1; % older adults
age_group(9:55) = 2; % younger adults
age_group = [age_group'; age_group'];

% set up labels for networks factor
networks = [];
networks(1:55) = 3; % individualized parcels
networks(56:110) = 4; % group parcels

factors = {age_group, networks};

%% comparing within-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
wFC = [Lifespan_wFC(:,1); iNet_wFC(:,1); Lifespan_wFC(:,2); iNet_wFC(:,2)];
aov_wFC = anova(factors, wFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group", "Network Parcellation"])

%% comparing between-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
bFC = [Lifespan_bFC(:,1); iNet_bFC(:,1); Lifespan_bFC(:,2); iNet_bFC(:,2)];
aov_bFC = anova(factors, bFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Network Parcellation"])

%% comparing segregation index
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
seg = [Lifespan_seg(:,1); iNet_seg(:,1); Lifespan_seg(:,2); iNet_seg(:,2)];
aov_seg = anova(factors, seg, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Network Parcellation"])

%% 2x2x@ ANOVAS: age group(young adults & older adults) x network parcellation (group avg & individualized) x network class (sensorimotor & association)
% set up labels for age_group factor
age_group = [];
age_group(1:8) = 1; % older adults
age_group(9:55) = 2; % younger adults
age_group = [age_group'; age_group'; age_group'; age_group'];

% set up labels for network parcellation factor
networks = [];
networks(1:55) = 3; % individualized parcels
networks(56:110) = 4; % group parcels
networks = [networks'; networks'];

% set up labels for networks class factor
class = [];
class(1:110) = 5; % sensorimotor systems
class(111:220) = 6; % association systems
factors = {age_group(:), networks(:), class'};

%% comparing within-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
SM_wFC = [Lifespan_SM.wFC(:,1); iNet_SM.wFC(:,1); Lifespan_SM.wFC(:,2); iNet_SM.wFC(:,2)];
assoc_wFC = [Lifespan_assoc.wFC(:,1); iNet_assoc.wFC(:,1); Lifespan_assoc.wFC(:,2); iNet_assoc.wFC(:,2)];
wFC = [SM_wFC; assoc_wFC];
aov_wFC = anova(factors, wFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Network Parcellation" "Network Class"])

%% comparing between-net FC
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
SM_bFC = [Lifespan_SM.bFC(:,1); iNet_SM.bFC(:,1); Lifespan_SM.bFC(:,2); iNet_SM.bFC(:,2)];
assoc_bFC = [Lifespan_assoc.bFC(:,1); iNet_assoc.bFC(:,1); Lifespan_assoc.bFC(:,2); iNet_assoc.bFC(:,2)];
bFC = [SM_bFC; assoc_bFC];
aov_bFC = anova(factors, bFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Network Parcellation" "Network Class"])

%% comparing segregation index
% prepare data: first Lifespan individual nets, second iNet individual
% nets, third Lifespan group nets, fourth iNet group nets
SM_seg = [Lifespan_SM.seg(:,1); iNet_SM.seg(:,1); Lifespan_SM.seg(:,2); iNet_SM.seg(:,2)];
assoc_seg = [Lifespan_assoc.seg(:,1); iNet_assoc.seg(:,1); Lifespan_assoc.seg(:,2); iNet_assoc.seg(:,2)];
seg = [SM_seg; assoc_seg];
aov_seg = anova(factors, seg, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Network Parcellation" "Network Class"])
