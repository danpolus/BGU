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
        
        if isempty(AvalancheAnalysisData.SimilarityMat(iTau).tau)
            continue;
        end
        
        ClusteringData(iTau).ClustersFull = [];
        ClusteringData(iTau).ClustersLen = [];
        ClusteringData(iTau).AffinityFull = [];
        ClusteringData(iTau).AffinityLen = [];
        ClusteringData(iTau).tau = AvalancheAnalysisData.SimilarityMat(iTau).tau;
        
        fig_name = ['Full Matrix - \tau '  num2str(AvalancheAnalysisData.SimilarityMat(iTau).tau)   '  ' files{iFile}];
        ClusteringData(iTau).ClustersFull = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatFull, fig_name ,1);
        ClusteringData(iTau).AffinityFull = avalanche_affinity(ClusteringData(iTau).ClustersFull, AvalancheAnalysisData.SimilarityMat(iTau).Id, fig_name, 1);
        for iLen = 1:length(AvalancheAnalysisData.SimilarityMat(iTau).MatLen)
            fig_name = ['Avalanche Length: ' num2str(iLen) ' \tau '  num2str(AvalancheAnalysisData.SimilarityMat(iTau).tau)   '  ' files{iFile}];
            ClusteringData(iTau).ClustersLen{iLen} = similarity_mat_2_clusters(AvalancheAnalysisData.SimilarityMat(iTau).MatLen{iLen}, fig_name, 1);
            ClusteringData(iTau).AffinityLen{iLen} = avalanche_affinity(ClusteringData(iTau).ClustersLen{iLen}, AvalancheAnalysisData.SimilarityMat(iTau).IdLen{iLen}, fig_name, 1);
        end
        
    end %for iTau
    
    save([fp files{iFile}(1:end-4) '_clusters.mat'],'ClusteringData');
    
    %plot pie chart of composing avalanches: per scenario per subject per condition per word
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


function Affinity = avalanche_affinity(Clusters_vec, Id_vec, fig_name, plot_flg)

subj_inxs = [8 9];

Affinity = [];
for iAch=1:length(Id_vec)
    subj_idx = subj_inxs(isstrprop(Id_vec{iAch}(subj_inxs),'digit'));
    Affinity{str2num(Id_vec{iAch}(4)),str2num(Id_vec{iAch}(subj_idx)),str2num(Id_vec{iAch}(13)),str2num(Id_vec{iAch}(17))} = [];
end

for iAch=1:length(Id_vec)
    subj_idx = subj_inxs(isstrprop(Id_vec{iAch}(subj_inxs),'digit'));
    Affinity{str2num(Id_vec{iAch}(4)),str2num(Id_vec{iAch}(subj_idx)),str2num(Id_vec{iAch}(13)),str2num(Id_vec{iAch}(17))} = ...
        [Affinity{str2num(Id_vec{iAch}(4)),str2num(Id_vec{iAch}(subj_idx)),str2num(Id_vec{iAch}(13)),str2num(Id_vec{iAch}(17))}  Clusters_vec(iAch)];
end

if plot_flg
    
    plot_id = [];
    for k = 1:size(Affinity,1)
        for l = 1:size(Affinity,2)
            for m = 1:size(Affinity,3)
                for n = 1:size(Affinity,4)
                    if ~isempty(Affinity{k,l,m,n})
                        plot_id = [plot_id; [k,l,m,n]];
                    end
                end
            end
        end
    end
    
    if ~isempty(plot_id)
        figure('Name',fig_name);
        for i=1:size(plot_id,1)
            uniqe_ach = unique(Affinity{plot_id(i,1),plot_id(i,2),plot_id(i,3),plot_id(i,4)});
            avalanch_type_cnt = [];
            for iAch=1:length(uniqe_ach)
                avalanch_type_cnt(iAch) = sum(Affinity{plot_id(i,1),plot_id(i,2),plot_id(i,3),plot_id(i,4)} == uniqe_ach(iAch));
            end
            pie_labels = string(uniqe_ach);
            if length(pie_labels) <= 1
                pie_labels = {pie_labels};
            end
            subplot(ceil(size(plot_id,1)/2),2,i);pie(avalanch_type_cnt,pie_labels);
            title(['scn: ' num2str(plot_id(i,1)) '  subj: ' num2str(plot_id(i,2)) '  cond: ' num2str(plot_id(i,3)) '  word: ' num2str(plot_id(i,4))]);
        end
    end
    
end
