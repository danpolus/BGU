%
%Step 6 : cluster_avalanches - hierarcial clustering; generate statistics
%https://www.mathworks.com/help/stats/hierarchical-clustering.html
%
%  inputs:
% MultiFileAchVecs - avalanche vectors from all sets
% SimilarityMat - similarity between avalanches distances matrix
% TrainingSets - avalanches (epochs) for training
% TrainValidTest - (just saves it)
% save_str - add this sting to filename
% saveFlg
% plotFlg
%
%  outputs:
% ClusteringDataSets - clusters and statistics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClusteringDataSets = cluster_avalanches(MultiFileAchVecs, SimilarityMat, TrainingSets, TrainValidTest, save_str, saveFlg, plotFlg)

params_t = global_params();

optimal_tau_t = MultiFileAchVecs{1}(~cellfun(@isempty,{MultiFileAchVecs{1}.is_optimal_tau}));
fileInfo = optimal_tau_t(([MultiFileAchVecs{1}.is_optimal_tau] == 1)).dataInfo.FileInfo;

if saveFlg
    output_fp = [fileInfo.base_fp '4 clusters\'];
    mkdir(output_fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

percent_waitbar = 0;
f_waitbar = waitbar(percent_waitbar, ['cluster avalanches ' num2str(100*percent_waitbar) '%'], 'Name',fileInfo.orig_fn );
        
ClusteringDataSets = [];
for iCrossValid = 1:length(TrainingSets)
    
    TrainSimilarityMat = get_training_similarity_mat(MultiFileAchVecs, SimilarityMat, TrainingSets{iCrossValid});       
    
    for iTau = 1:length(TrainSimilarityMat)
        
        if isempty(TrainSimilarityMat(iTau).tau)
            continue;
        end
        
        ClusteringData(iTau).Id = [];
        ClusteringData(iTau).Clusters = [];
        ClusteringData(iTau).Stats = [];
        ClusteringData(iTau).tau = TrainSimilarityMat(iTau).tau;
        ClusteringData(iTau).is_optimal_tau = TrainSimilarityMat(iTau).is_optimal_tau;
        
        nofMat = length(TrainSimilarityMat(iTau).Mat);
        concat_max_nof_clusters_Len = 0;
        for iLen = 1:nofMat
            percent_waitbar = (iCrossValid - 1 +iLen/(nofMat+1))/length(TrainingSets);
            waitbar(percent_waitbar,f_waitbar,['cross-valid ' num2str(iCrossValid) '  cluster avalanches \tau=' num2str(TrainSimilarityMat(iTau).tau) '   ' num2str(100*percent_waitbar) '%']);
            
            if iLen < nofMat
                fig_name = ['Avalanche Length: ' num2str(iLen) ' \tau '  num2str(ClusteringData(iTau).tau)   '  ' fileInfo.orig_fn];
            else
                fig_name = ['Full Matrix - \tau '  num2str(ClusteringData(iTau).tau)   '  ' fileInfo.orig_fn];
            end
            ClusteringData(iTau).Id{iLen} = TrainSimilarityMat(iTau).Id{iLen};
            ClusteringData(iTau).Clusters{iLen} = similarity_mat_2_clusters(TrainSimilarityMat(iTau).Mat{iLen}, fig_name, params_t, plotFlg);
            ClusteringData(iTau).Stats{iLen} = clusters_statistics(ClusteringData(iTau).Clusters{iLen}, ClusteringData(iTau).Id{iLen}, TrainingSets{iCrossValid}(iTau), fig_name, 0);
            if iLen < nofMat
                concat_max_nof_clusters_Len = max([concat_max_nof_clusters_Len; ClusteringData(iTau).Clusters{iLen}]);
            end
        end
        
        concat_cluster_id_prefix_Len = 10^ceil(log10(concat_max_nof_clusters_Len));
        ClusteringData(iTau).Id{nofMat+1} = [];
        ClusteringData(iTau).Clusters{nofMat+1} = [];
        for iLen = 1:nofMat-1
            ClusteringData(iTau).Clusters{nofMat+1} = [ClusteringData(iTau).Clusters{nofMat+1};  (concat_cluster_id_prefix_Len*iLen)+ClusteringData(iTau).Clusters{iLen}];
            ClusteringData(iTau).Id{nofMat+1} = [ClusteringData(iTau).Id{nofMat+1}  ClusteringData(iTau).Id{iLen}];
        end
        fig_name = ['Concatinated Different Length - \tau '  num2str(ClusteringData(iTau).tau)   '  ' fileInfo.orig_fn];
        ClusteringData(iTau).Stats{nofMat+1} = clusters_statistics(ClusteringData(iTau).Clusters{nofMat+1}, ClusteringData(iTau).Id{nofMat+1}, TrainingSets{iCrossValid}(iTau), fig_name, plotFlg);

    end %for iTau
    ClusteringDataSets{iCrossValid} = ClusteringData;
    
end %for iCrossValid

close(f_waitbar);

if saveFlg
    save([output_fp fileInfo.orig_fn '_' save_str 'Clusters.mat'],'ClusteringDataSets','MultiFileAchVecs','TrainValidTest');
end

% tau = 1;  av_len = 3;
% figure; imagesc(SimilarityMat(tau).Mat{end}); title(['all avalanches similarity,  \tau =' num2str(SimilarityMat(tau).tau)]);
% figure; imagesc(SimilarityMat(tau).Mat{av_len}); title([num2str(av_len) ' bins length avalanches similarity,  \tau =' num2str(SimilarityMat(tau).tau)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TrainSimilarityMat = get_training_similarity_mat(MultiFileAchVecs, SimilarityMat, TrainingSet)

TrainSimilarityMat = SimilarityMat;

for iTau = 1:length(SimilarityMat)
    if isempty(SimilarityMat(iTau).tau)
        continue;
    end
    
    nofMat = length(SimilarityMat(iTau).Mat);
    
    TrainingIds{nofMat} = [];
    for iCond = 1:length(TrainingSet(iTau).CondIds)
        for iEpoch = 1:length(TrainingSet(iTau).EpochIds{iCond})
            fln = str2num(TrainingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
            epc = str2num(TrainingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
            ach_vectors_t = MultiFileAchVecs{fln}(iTau).epochs_vecs{epc};
            for iVec = 1:length(ach_vectors_t)
                TrainingIds{ach_vectors_t(iVec).length_bins} = [TrainingIds{ach_vectors_t(iVec).length_bins} {ach_vectors_t(iVec).id}];
                TrainingIds{nofMat} = [TrainingIds{nofMat} {ach_vectors_t(iVec).id}];          
            end
        end
    end   
    
    for iLen = 1:nofMat
        if isempty(TrainingIds{iLen}) 
            TrainSimilarityMat(iTau).Mat{iLen} = [];
            TrainSimilarityMat(iTau).Id{iLen} = [];
        else
            ach_idx = contains(SimilarityMat(iTau).Id{iLen}, TrainingIds{iLen});
            TrainSimilarityMat(iTau).Mat{iLen} = SimilarityMat(iTau).Mat{iLen}(ach_idx,ach_idx);
            TrainSimilarityMat(iTau).Id{iLen} = SimilarityMat(iTau).Id{iLen}(ach_idx);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_opt = similarity_mat_2_clusters(sim_M, fig_name, params_t, plotFlg)

similarity_dist_threshold = 0.5; %for cases when there is same distance all over the matrix

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
            T_opt = ones(N,1); %single cluster
        else
            T_opt = (1:N)'; %all separate clusters
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
        
        if C >= params_t.minimal_contrast
            NofsClusters = [NofsClusters max(T)];
            Cutoffs = [Cutoffs cutoff];
            Contrasts = [Contrasts C];
        end
    end
    %     cutoff_opt = 0.9;
    %     T_opt = cluster(cluster_tree,'cutoff',cutoff_opt,'criterion','distance');
    %     inconsistency_cutoff = 1.2; inconsistency_depth = 2;
    %     T_ = cluster(cluster_tree,'cutoff',inconsistency_cutoff,'depth',inconsistency_depth);
    %     nof_clusters = 2;
    %     T_ = cluster(cluster_tree,'maxclust',nof_clusters);
    
    if length(NofsClusters)>1
        [NofsClusters, sort_inces] = sort(NofsClusters);
        Contrasts = Contrasts(sort_inces);
        Cutoffs = Cutoffs(sort_inces);
        NofsClusters_scaled = log10(NofsClusters)/log10(params_t.reasonable_nof_clusters);
        %NofsClusters_scaled = (NofsClusters_scaled - log10(2)/log10(params_t.reasonable_nof_clusters)) / (1-log10(2)/log10(params_t.reasonable_nof_clusters)); %scale log[2:reasonable_nof_clusters]
        if min(NofsClusters_scaled)<1
            NofsClusters_scaled = (NofsClusters_scaled - min(NofsClusters_scaled)) / (1-min(NofsClusters_scaled)); %scale log[min dynamic range:reasonable_nof_clusters]
        else
            error('min(NofsClusters_scaled)>=1');
        end
        %Contrasts_scaled = (Contrasts-params_t.minimal_contrast) / (1-params_t.minimal_contrast); %mininal contrast scaling
        Contrasts_scaled = (Contrasts-min(Contrasts)) / (max(Contrasts)-min(Contrasts)); %dynamic range scaling
        if isinf(params_t.max_nof_clusters)
            [~, opt_inx] = min(sqrt((0-NofsClusters_scaled).^2 + (1-Contrasts_scaled).^2)); % find minimal distance to (0,1)
        else
            [~, opt_inx] = max(Contrasts(NofsClusters <= params_t.max_nof_clusters));
            if isempty(opt_inx)
                opt_inx = 1;
            end
        end
    elseif length(NofsClusters) == 1
        sort_inces = 1; opt_inx = 1;
    else %clustering failed due to low contrast
        if max(clusters_distance) <= similarity_dist_threshold
            T_opt = ones(N,1); %single cluster
        else
            T_opt = (1:N)'; %all separate clusters
        end
        return;
    end
    
    T_opt = cluster(cluster_tree,'cutoff',Cutoffs(opt_inx),'criterion','distance');
    
    if plotFlg
        figure('Name',fig_name);
        subplot(2,2,1);imagesc(sim_M);title('original');colorbar;
        subplot(2,2,3);plot(log10(NofsClusters),Contrasts, log10(NofsClusters(opt_inx)),Contrasts(opt_inx), 'xr');xlabel('log10(nof clusters)');ylabel('contrast');title('Contrast');
        subplot(2,2,4);[~,T_dend,perm_dend] = dendrogram(cluster_tree, NofsClusters(opt_inx), 'ColorThreshold',Cutoffs(opt_inx), 'Orientation','left');
        title(['cutoff distance = ' num2str(Cutoffs(opt_inx)) '  contrast = ' num2str(Contrasts(opt_inx)) '  nof clusters = ' num2str(NofsClusters(opt_inx))]);
        T_perm = [];
        for iClust = perm_dend
            T_perm = [T_perm; find(T_dend==iClust)];
        end
        subplot(2,2,2);imagesc(sim_M(T_perm,T_perm));title('clustered');colorbar;
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute statistics based in clustering:  P(cluster), P(conditon), P(cluster|conditon), P(conditon|cluster)
function Stats = clusters_statistics(Clusters_vec, Id_vec, TrainingSet, fig_name, plotFlg)

nof_ach = length(Id_vec);

Stats.CondIds = TrainingSet.CondIds;
Stats.P_cond = zeros(1,length(TrainingSet.CondIds));
Stats.P_clst = [];
Stats.P_clstGINVcond = [];
Stats.P_clondGINVclst = [];

%dataset nof epoch based
for iCond= 1:length(Stats.CondIds)
    Stats.P_cond(iCond) = length(TrainingSet.EpochIds{iCond});
end
Stats.P_cond = Stats.P_cond/sum(Stats.P_cond);
% Stats.P_cond = ones(1,length(Stats.CondIds))/length(Stats.CondIds); %uniform

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

if plotFlg
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
