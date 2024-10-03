%% how many are there
for s=1:length(subjects)
    start=1+sum(subjects_counts(1:s-1));
    fin=start+subjects_counts(s)-1;
    for i=1:length(cluster_labels)
        pct_hub_cluster_counts(s,i) = 100*length(find(idx(start:fin)==i))/length(idx(start:fin));
        subject_hubs{s}=idx(start:fin);
    end
end
figure(1)
b=barh(categorical(subjects),pct_hub_cluster_counts, 'stacked','FaceColor','flat');
for i = 1:size(cluster_labels,2)
    b(i).CData = colors(i,:);
end
xlim([0, 100]);
xlabel('Percent of Hubs');
saveas(figure(1),[outdir '/hub_subtypes_subject_counts.jpg'], 'jpg')