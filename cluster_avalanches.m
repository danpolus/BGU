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
    AvalancheAnalysisData = load([fp files{iFile}]);
    
    for iTau = 1:length(AvalancheAnalysisData.SimilarityMat)
        
        ClusteringData(iTau).ClustersFull = [];
        ClusteringData(iTau).ClustersLen = [];
        ClusteringData(iTau).tau = AvalancheAnalysisData.SimilarityMat(iTau).tau;
        
        ClusteringData(iTau).ClustersFull = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatFull);
        for iLen = 1:length(AvalancheAnalysisData.SimilarityMat(iTau).MatLen)
            ClusteringData(iTau).ClustersLen{iLen} = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatLen{iLen});
        end
        
    end %for iTau
    
    save([fp files{iFile}(1:end-4) '_clusters.mat'],'ClusteringData');
end

% tau = 1;  av_len = 3;
% figure; imagesc(AvalancheAnalysisData.SimilarityMat(tau).MatFull); title(['all avalanches similarity,  \tau =' num2str(AvalancheAnalysisData.SimilarityMat(tau).tau)]);
% figure; imagesc(AvalancheAnalysisData.SimilarityMat(tau).MatLen{av_len}); title([num2str(av_len) ' bins length avalanches similarity,  \tau =' num2str(AvalancheAnalysisData.SimilarityMat(tau).tau)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_opt = similarity_mat_2_clusters(sim_M)

if numel(sim_M) == 0
    T_opt = [];
elseif numel(sim_M) == 1
    T_opt = 1;
else
    sim_M = min(1,sim_M);
    dist_v = squareform(~eye(size(sim_M)).*(1-sim_M),'tovector'); %convert to ditance, put zeros on diagonal, convert to vector
    cluster_tree = linkage(dist_v,'average');
    
    %divide to clusters
    sim_M_zerodiag = ~eye(size(sim_M)).*sim_M; %put zeros on diagonal
    clusters_distance = unique(cluster_tree(:,3));
    max_contrast = 0;
    cutoff_opt = 0;
    T_opt = 1:size(sim_M,1);
    for iCut=2:length(clusters_distance)
        cutoff = (clusters_distance(iCut)+clusters_distance(iCut-1))/2;
        T = cluster(cluster_tree,'cutoff',cutoff,'criterion','distance');
        D_in = 0; D_in_cnt = 0;
        for iClust = unique(T)'
            D_in = D_in + sum(sum(sim_M_zerodiag(T==iClust, T==iClust)));
            D_in_cnt = D_in_cnt + sum(T==iClust)^2;
        end
        D_out = (sum(sum(sim_M_zerodiag)) - D_in)/(numel(sim_M_zerodiag) - D_in_cnt);
        D_in = D_in/D_in_cnt;
        C = (D_in-D_out)/(D_in+D_out);
        if C > max_contrast
            max_contrast = C;
            cutoff_opt = cutoff;
            T_opt = T;
        end
    end
%     cutoff_opt = 0.9;
%     T_opt = cluster(cluster_tree,'cutoff',cutoff_opt,'criterion','distance');    
    
    %     inconsistency_cutoff = 1.2; inconsistency_depth = 2;
    %     T_ = cluster(cluster_tree,'cutoff',inconsistency_cutoff,'depth',inconsistency_depth);
    %     nof_clusters = 2;
    %     T_ = cluster(cluster_tree,'maxclust',nof_clusters);
     
    %debug plots
    subplot(1,3,1);imagesc(sim_M);title('original');
    subplot(1,3,2);[~,T_dend,perm_dend] = dendrogram(cluster_tree,max(T_opt), 'ColorThreshold',cutoff_opt, 'Orientation','left');title('dendogram with optimal cutoff');
    T_perm = [];
    for iClust = perm_dend
        T_perm = [T_perm; find(T_dend==iClust)];
    end
    subplot(1,3,3);imagesc(sim_M(T_perm,T_perm));title('after hierarchical clustering');
 
    %     mtd = 'average';  %'average' 'weighted'
    %     cluster_tree = linkage(dist_v,mtd);
    %     c = cophenet(cluster_tree,dist_v)
    %     inconsistency_M = inconsistent(cluster_tree,inconsistency_depth);
    %     leafOrder = optimalleaforder(cluster_tree,dist_v, 'Criteria','adjacent'); %'adjacent' 'group'
    %     cutoff = median([cluster_tree(end-2,3) cluster_tree(end-1,3)]); %'default'
    %     dendogram_leafs_num = 10; % 0  30
    %     dendrogram(cluster_tree,dendogram_leafs_num, 'ColorThreshold',cutoff, 'Reorder',leafOrder, 'CheckCrossing',true); title(['method - ' mtd]);
end
