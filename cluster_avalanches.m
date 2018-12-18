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

fig = figure;

for iFile = 1:length(files)
    AvalancheAnalysisData = load([fp files{iFile}]);
    
    for iTau = 1:length(AvalancheAnalysisData.SimilarityMat)
        
        ClusteringData(iTau).ClustersFull = [];
        ClusteringData(iTau).ClustersLen = [];
        ClusteringData(iTau).tau = AvalancheAnalysisData.SimilarityMat(iTau).tau;
        
        ClusteringData(iTau).ClustersFull = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatFull,1);
        fig.Name = 'Full Matrix';
        for iLen = 1:length(AvalancheAnalysisData.SimilarityMat(iTau).MatLen)
            ClusteringData(iTau).ClustersLen{iLen} = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatLen{iLen},1);
            fig.Name = ['Avalanche Length: ' num2str(iLen)];
        end
        
    end %for iTau
    
    save([fp files{iFile}(1:end-4) '_clusters.mat'],'ClusteringData');
end

% tau = 1;  av_len = 3;
% figure; imagesc(AvalancheAnalysisData.SimilarityMat(tau).MatFull); title(['all avalanches similarity,  \tau =' num2str(AvalancheAnalysisData.SimilarityMat(tau).tau)]);
% figure; imagesc(AvalancheAnalysisData.SimilarityMat(tau).MatLen{av_len}); title([num2str(av_len) ' bins length avalanches similarity,  \tau =' num2str(AvalancheAnalysisData.SimilarityMat(tau).tau)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_opt = similarity_mat_2_clusters(sim_M, plot_flg)

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
        if clusters_distance >= similarity_dist_threshold
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
    else
        sort_inces = 1; opt_inx = 1;
    end       
%     NofsClusters(opt_inx)
%     Cutoffs(sort_inces(opt_inx))
%     Contrasts(sort_inces(opt_inx))
    T_opt = cluster(cluster_tree,'cutoff',Cutoffs(sort_inces(opt_inx)),'criterion','distance');
     
    if plot_flg
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
