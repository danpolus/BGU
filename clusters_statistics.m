%
% compute statistics based in clustering:
% P(cluster), P(conditon), P(cluster|conditon), P(conditon|cluster)
%
function Stats = clusters_statistics(Clusters_vec, Id_vec, MultiFileAchVecs, TestingSet, iTau, fig_name, plot_flg)

nof_ach = length(Id_vec);

Stats.CondIds = TestingSet.CondIds;
Stats.P_cond = zeros(1,length(TestingSet.CondIds));
Stats.P_clst = [];
Stats.P_clstGINVcond = [];
Stats.P_clondGINVclst = [];

%dataset nof epoch based
for iFile = 1:length(MultiFileAchVecs)
    cond_idx = find(contains(Stats.CondIds, MultiFileAchVecs{iFile}(iTau).file_id(1:17)));
    Stats.P_cond(cond_idx) = Stats.P_cond(cond_idx) + length(MultiFileAchVecs{iFile}(iTau).epochs_vecs) - length(TestingSet.EpochIds{cond_idx});
end
Stats.P_cond = Stats.P_cond/sum(Stats.P_cond);
Stats.P_cond_uni = ones(1,length(MultiFileAchVecs))/length(MultiFileAchVecs); %uniform

clstVScond = zeros(length(Stats.CondIds),max(Clusters_vec));
for iAch=1:nof_ach
    cond_idx = find(contains(Stats.CondIds, Id_vec{iAch}(1:17)));
    clstVScond(cond_idx,Clusters_vec(iAch)) = clstVScond(cond_idx,Clusters_vec(iAch)) + 1;
end

% Stats.P_cond = sum(clstVScond,2)'/nof_ach; %avalanche based
Stats.P_clst = sum(clstVScond,1)/nof_ach;
Stats.P_clondGINVclst = clstVScond ./ sum(clstVScond,1);
Stats.P_clstGINVcond = clstVScond ./ sum(clstVScond,2);

Stats.P_cond(isnan(Stats.P_cond)) = 0;
Stats.P_clst(isnan(Stats.P_clst)) = 0;
Stats.P_clondGINVclst(isnan(Stats.P_clondGINVclst)) = 0;
Stats.P_clstGINVcond(isnan(Stats.P_clstGINVcond)) = 0;

if plot_flg
    figure('Name',fig_name);
    for i=1:size(clstVScond,1)  
        pie_labels = string(1:size(clstVScond,2));
        if length(pie_labels) <= 1
            pie_labels = {pie_labels};
        end
        if sum(clstVScond(i,:)) > 0
            subplot(ceil(size(clstVScond,1)/2),2,i); pie(clstVScond(i,:),pie_labels);
            title(Stats.CondIds(i));
        end
    end
end
