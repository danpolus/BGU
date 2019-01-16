%
%Step 4
%
% https://www.mathworks.com/help/stats/hierarchical-clustering.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_avalanches()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select avalanche similarity matrix files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end


for iFile = 1:length(files)
    tic
    AvalancheAnalysisData = load([fp files{iFile}]);
    
    for iTau = 1:length(AvalancheAnalysisData.SimilarityMat)
        
        if isempty(AvalancheAnalysisData.SimilarityMat(iTau).tau)
            continue;
        end
        
        ClusteringData(iTau).ClustersFull = [];
        ClusteringData(iTau).ClustersLen = [];
        ClusteringData(iTau).StatsFull = [];
        ClusteringData(iTau).StatsLen = [];
        ClusteringData(iTau).tau = AvalancheAnalysisData.SimilarityMat(iTau).tau;
        
        fig_name = ['Full Matrix - \tau '  num2str(AvalancheAnalysisData.SimilarityMat(iTau).tau)   '  ' files{iFile}];
        ClusteringData(iTau).ClustersFull = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatFull, fig_name ,1);
        ClusteringData(iTau).StatsFull = clusters_statistics(ClusteringData(iTau).ClustersFull, AvalancheAnalysisData.SimilarityMat(iTau).Id, fig_name, 1);
        
        concat_max_nof_clusters_Len = 0;
        for iLen = 1:length(AvalancheAnalysisData.SimilarityMat(iTau).MatLen)
            fig_name = ['Avalanche Length: ' num2str(iLen) ' \tau '  num2str(AvalancheAnalysisData.SimilarityMat(iTau).tau)   '  ' files{iFile}];
            ClusteringData(iTau).ClustersLen{iLen} = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatLen{iLen}, fig_name, 1);
            ClusteringData(iTau).StatsLen{iLen} = clusters_statistics(ClusteringData(iTau).ClustersLen{iLen}, AvalancheAnalysisData.SimilarityMat(iTau).IdLen{iLen}, fig_name, 0);
            concat_max_nof_clusters_Len = max([concat_max_nof_clusters_Len; ClusteringData(iTau).ClustersLen{iLen}]);
        end
        
        concat_cluster_id_prefix_Len = 10^ceil(log10(concat_max_nof_clusters_Len));
        concat_ClustersLen = [];
        concat_IdLen = [];
        for iLen = 1:length(AvalancheAnalysisData.SimilarityMat(iTau).MatLen)
            concat_ClustersLen = [concat_ClustersLen;  (concat_cluster_id_prefix_Len*iLen)+ClusteringData(iTau).ClustersLen{iLen}];
            concat_IdLen = [concat_IdLen  AvalancheAnalysisData.SimilarityMat(iTau).IdLen{iLen}];
        end
        fig_name = ['Concatinated Different Length - \tau '  num2str(AvalancheAnalysisData.SimilarityMat(iTau).tau)   '  ' files{iFile}];
        ClusteringData(iTau).StatsConcatLen = clusters_statistics(concat_ClustersLen, concat_IdLen, fig_name, 1);
        
    end %for iTau
    
    save([fp files{iFile}(1:end-4) '_clusters.mat'],'ClusteringData');
    toc
end

% tau = 1;  av_len = 3;
% figure; imagesc(AvalancheAnalysisData.SimilarityMat(tau).MatFull); title(['all avalanches similarity,  \tau =' num2str(AvalancheAnalysisData.SimilarityMat(tau).tau)]);
% figure; imagesc(AvalancheAnalysisData.SimilarityMat(tau).MatLen{av_len}); title([num2str(av_len) ' bins length avalanches similarity,  \tau =' num2str(AvalancheAnalysisData.SimilarityMat(tau).tau)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_opt = similarity_mat_2_clusters(sim_M, fig_name, plot_flg)

max_nof_clusters = Inf; % 15 Inf
similarity_dist_threshold = 0.5; %for cases when there is same distance all over the matrix
contrast_tolerance = 0.05; %allow contrast to be 5% (of dynamic range) less to nminimize nof clusters

if numel(sim_M) == 0
    T_opt = [];
elseif numel(sim_M) == 1
    T_opt = 1;
else
    N = size(sim_M,1);
    sim_M = min(1,sim_M);
    sim_M = ~eye(N,N).*sim_M; %put zeros on diagonal
    dist_v = squareform(~eye(N,N).*(1-sim_M),'tovector'); %convert to ditance, put zeros on diagonal, convert to vector
    cluster_tree = linkage(dist_v,'average');
    
    %divide to clusters
    clusters_distance = unique(cluster_tree(:,3));
    if length(clusters_distance) == 1
        if clusters_distance <= similarity_dist_threshold
            T_opt = ones(N,1);
        else
            T_opt = (1:N)';
        end
        return;
    end
    
    NofsClusters = []; Cutoffs = []; Contrasts = [];
    for iCut=2:length(clusters_distance)
        cutoff = (clusters_distance(iCut)+clusters_distance(iCut-1))/2;
        T = cluster(cluster_tree,'cutoff',cutoff,'criterion','distance');
        D_in = 0; D_in_cnt = 0;
        for iClust = unique(T)'
            D_in = D_in + sum(sum(sim_M(T==iClust, T==iClust)));
            D_in_cnt = D_in_cnt + sum(T==iClust)^2;
        end
        D_out = (sum(sum(sim_M)) - D_in)/(N*N - D_in_cnt);
        D_in = D_in/D_in_cnt;
        C = (D_in-D_out)/(D_in+D_out);
        
        NofsClusters = [NofsClusters max(T)];
        Cutoffs = [Cutoffs cutoff];
        Contrasts = [Contrasts C];
    end
    %     cutoff_opt = 0.9;
    %     T_opt = cluster(cluster_tree,'cutoff',cutoff_opt,'criterion','distance');
    %     inconsistency_cutoff = 1.2; inconsistency_depth = 2;
    %     T_ = cluster(cluster_tree,'cutoff',inconsistency_cutoff,'depth',inconsistency_depth);
    %     nof_clusters = 2;
    %     T_ = cluster(cluster_tree,'maxclust',nof_clusters);
    
    if length(NofsClusters)>1
        [NofsClusters, sort_inces] = sort(NofsClusters);
        opt_inx = find(Contrasts(sort_inces) > max(Contrasts) - (max(Contrasts)-min(Contrasts))*contrast_tolerance, 1);
        opt_inx = min(opt_inx, find(NofsClusters <= max_nof_clusters, 1, 'last'));
    else
        sort_inces = 1; opt_inx = 1;
    end
    
    %     NofsClusters(opt_inx)
    %     Cutoffs(sort_inces(opt_inx))
    %     Contrasts(sort_inces(opt_inx))
    T_opt = cluster(cluster_tree,'cutoff',Cutoffs(sort_inces(opt_inx)),'criterion','distance');
    
    if plot_flg
        figure('Name',fig_name);
        subplot(1,3,1);imagesc(sim_M);title('original');
        subplot(1,3,2);[~,T_dend,perm_dend] = dendrogram(cluster_tree, NofsClusters(opt_inx), 'ColorThreshold',Cutoffs(sort_inces(opt_inx)), 'Orientation','left');
        title(['cutoff distance = ' num2str(Cutoffs(sort_inces(opt_inx))) '  contrast = ' num2str(Contrasts(sort_inces(opt_inx))) '  nof clusters = ' num2str(NofsClusters(opt_inx))]);
        T_perm = [];
        for iClust = perm_dend
            T_perm = [T_perm; find(T_dend==iClust)];
        end
        subplot(1,3,3);imagesc(sim_M(T_perm,T_perm));title('clustered');
        
        %     mtd = 'average';  %'average' 'weighted'
        %     cluster_tree = linkage(dist_v,mtd);
        %     c = cophenet(cluster_tree,dist_v)
        %     inconsistency_M = inconsistent(cluster_tree,inconsistency_depth);
        %     leafOrder = optimalleaforder(cluster_tree,dist_v, 'Criteria','adjacent'); %'adjacent' 'group'
        %     cutoff = median([cluster_tree(end-2,3) cluster_tree(end-1,3)]); %'default'
        %     dendogram_leafs_num = 10; % 0  30
        %     dendrogram(cluster_tree,dendogram_leafs_num, 'ColorThreshold',cutoff, 'Reorder',leafOrder, 'CheckCrossing',true); title(['method - ' mtd]);
        
        %             figure;plot(NofsClusters,Contrasts(sort_inces), NofsClusters(opt_inx),Contrasts(sort_inces(opt_inx)),'xk');
    end
end


function Stats = clusters_statistics(Clusters_vec, Id_vec, fig_name, plot_flg)

subj_inxs = [8 9];

nof_ach = length(Id_vec);

Stats.CondIds = {};
Stats.P_clst = [];
Stats.P_cond = [];
Stats.P_clstGINVcond = [];
Stats.P_clondGINVclst = [];

for iAch=1:nof_ach
    if sum(strcmp(Stats.CondIds, Id_vec{iAch}(1:17))) == 0
        Stats.CondIds = [Stats.CondIds {Id_vec{iAch}(1:17)}];
    end
end

clstVScond = zeros(length(Stats.CondIds),max(Clusters_vec));
for iAch=1:nof_ach
    cond_idx = find(contains(Stats.CondIds, Id_vec{iAch}(1:17)));
    clstVScond(cond_idx,Clusters_vec(iAch)) = clstVScond(cond_idx,Clusters_vec(iAch)) + 1;
end

Stats.P_clst = sum(clstVScond,1)/nof_ach;
Stats.P_cond = sum(clstVScond,2)/nof_ach;
Stats.P_clondGINVclst = clstVScond ./ sum(clstVScond,1);
Stats.P_clstGINVcond = clstVScond ./ sum(clstVScond,2);

if plot_flg
    figure('Name',fig_name);
    for i=1:size(clstVScond,1)  
        pie_labels = string(1:size(clstVScond,2));
        if length(pie_labels) <= 1
            pie_labels = {pie_labels};
        end
        subplot(ceil(size(clstVScond,1)/2),2,i); pie(clstVScond(i,:),pie_labels);
        title(Stats.CondIds(i));
    end
end
