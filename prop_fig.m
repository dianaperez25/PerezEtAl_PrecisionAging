%% how many are there
dataset = {'Lifespan-NU', 'iNet-FSU', 'Lifespan-FSU', 'iNet-NU'};
exclude_subs = {'LS46', 'INET108', 'LS108', 'INET057'};
FC_metrics_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/FCmetrics';
outdir = '/Users/diana/Desktop/';
rgb_colors = [1 0 0;
                  0 0 .6;
                  1 1 0;
                  0 .8 0;
                  0 .6 .6;
                  0 0 0;
                  .3 0 .6;
                  .2 1 1;
                  1 .5 0;
                  .6 .2 1;
                  .2 1 .2;
                  0 .2 .4;
                  0 0 1;
                  .8 .8 .6];

for d = 1:numel(dataset)
    close all;
    net_surf_area = [];
    [subjects, sessions, N] = get_subjects(dataset{d}, exclude_subs);

    for s=1:length(subjects)
        load(sprintf('%s/sub-%s_FCmetrics.mat', FC_metrics_dir, subjects{s}));
        %% calculate weights for each network
        indiv_nets_size = out_data.indiv(1,:);
        group_nets_size = out_data.group(1,:);
        total_surf_indiv = sum(indiv_nets_size);
        total_surf_group = sum(group_nets_size);
%         start=1+sum(subjects_counts(1:s-1));
%         fin=start+subjects_counts(s)-1;
        for i=1:14
            net_surf_area(s,i) = 100*(indiv_nets_size(i)/total_surf_indiv);
            %subject_hubs{s}=idx(start:fin);
        end
        if s == length(subjects)
            subjects{end+1} = 'Group';
            for i=1:14
            net_surf_area(s+1,i) = 100*(group_nets_size(i)/total_surf_group);
            %subject_hubs{s}=idx(start:fin);
            end
        end
    end
    figure(1)
    b=barh(categorical(subjects),net_surf_area, 'stacked','FaceColor','flat');
    for i = 1:14
        b(i).CData = rgb_colors(i,:);
    end
    xlim([0, 100]);
    xlabel('Percent Surface Area');
    saveas(figure(1),[outdir dataset{d} '_network_surface_area_proportions.jpg'], 'jpg')
end