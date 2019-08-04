%
%Step 6 : get_testing_clusters - find clusters for training and testing set vectors
%
%  inputs:
% ClusteringData - clusters and statistics
% MultiFileAchVecs - avalanche vectors from all sets
% TestingSet - avalanches (epochs) for testing
% saveFlg
%
%  outputs:
% TestingSetClusters - testing set clusters
% TrainingSetClusters - training set clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TestingSetClusters, TrainingSetClusters] = get_testing_clusters(ClusteringData, MultiFileAchVecs, TestingSet, saveFlg)

params_t = global_params();

optimal_tau_t = MultiFileAchVecs{1}(~cellfun(@isempty,{MultiFileAchVecs{1}.is_optimal_tau}));
fileInfo = optimal_tau_t(([MultiFileAchVecs{1}.is_optimal_tau] == 1)).dataInfo.FileInfo;

if saveFlg
    output_fp = [fileInfo.base_fp '5 testing\'];
    mkdir(output_fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TestingSetClusters = [];
TrainingSetClusters = [];
for iTau = 1:length(ClusteringData)
    
    if isempty(ClusteringData(iTau).tau)
        continue;
    end
    
    percent_waitbar = 0;
    f_waitbar = waitbar(percent_waitbar, ['get testing set clusters \tau=' num2str(TestingSet(iTau).tau) '   ' num2str(100*percent_waitbar) '%'], 'Name',fileInfo.orig_fn );
    
    %testing set
    TestingSetClusters(iTau).CondClst = [];
    TestingSetClusters(iTau).CondIds = TestingSet(iTau).CondIds;
    TestingSetClusters(iTau).tau = TestingSet(iTau).tau;
    TestingSetClusters(iTau).is_optimal_tau = TestingSet(iTau).is_optimal_tau;
    TestingSetClusters(iTau).fileInfo.base_fp = fileInfo.base_fp;
    TestingSetClusters(iTau).fileInfo.orig_fn = fileInfo.orig_fn;
    
    testEpochMat = zeros(length(TestingSetClusters(iTau).CondIds),params_t.max_nof_epochs);
    for iCond = 1:length(TestingSet(iTau).CondIds)
        %test epochs vectors
        for iEpoch = 1:length(TestingSet(iTau).EpochIds{iCond})
            percent_waitbar = (iCond - 1 + iEpoch/length(TestingSet(iTau).EpochIds{iCond}))/length(TestingSet(iTau).CondIds);
            waitbar(percent_waitbar,f_waitbar,['get testing set clusters \tau=' num2str(TestingSet(iTau).tau) '   ' num2str(100*percent_waitbar) '%']);
            
            fln = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
            epc = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
            ach_vectors_t = MultiFileAchVecs{fln}(iTau).epochs_vecs{epc};
            for iVec = 1:length(ach_vectors_t)
                [TestingSetClusters(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec,:), TestingSetClusters(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_sim(iVec,:)] = ...
                    find_vector_cluster(ach_vectors_t(iVec).vec, ach_vectors_t(iVec).length_bins, MultiFileAchVecs, ClusteringData, params_t.similarity_method, iTau, 0);
            end
            testEpochMat(fln,epc) = 1;
        end
    end
    
    %traning set
    TrainingSetClusters(iTau).CondClst = [];
    TrainingSetClusters(iTau).CondIds = TestingSet(iTau).CondIds;
    TrainingSetClusters(iTau).tau = TestingSet(iTau).tau;
    TrainingSetClusters(iTau).is_optimal_tau = TestingSet(iTau).is_optimal_tau;
    TrainingSetClusters(iTau).fileInfo.base_fp = fileInfo.base_fp;
    TrainingSetClusters(iTau).fileInfo.orig_fn = fileInfo.orig_fn;
    
    nofLen = length(ClusteringData(iTau).Clusters);
    iFullEpoch = zeros(1,length(TrainingSetClusters(iTau).CondIds));
    for iFile = 1:length(MultiFileAchVecs)
        cond_idx = find(contains(TrainingSetClusters(iTau).CondIds, MultiFileAchVecs{iFile}(iTau).file_id(1:17)));        
        for iEpoch = randperm(length(MultiFileAchVecs{iFile}(iTau).epochs_vecs))
            if testEpochMat(iFile,iEpoch) == 1 || isempty(MultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch})
                continue;
            end
            iFullEpoch(cond_idx) = iFullEpoch(cond_idx) + 1;
            avch_length_bins = [MultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch}.length_bins];
            for iAvalanche=1:length(avch_length_bins)
                id = [MultiFileAchVecs{iFile}(iTau).file_id 'epc' num2str(iEpoch,'%03d') 'ach' num2str(iAvalanche,'%04d')];
                cluster_num = zeros(1,nofLen);
                cluster_sim = zeros(1,nofLen);                
                ach_idx = contains(ClusteringData(iTau).Id{avch_length_bins(iAvalanche)}, id);
                cluster_num(avch_length_bins(iAvalanche)) = ClusteringData(iTau).Clusters{avch_length_bins(iAvalanche)}(ach_idx);
                ach_idx = contains(ClusteringData(iTau).Id{end-1}, id);
                cluster_num(nofLen-1) = ClusteringData(iTau).Clusters{nofLen-1}(ach_idx);
                ach_idx = contains(ClusteringData(iTau).Id{nofLen}, id);
                cluster_num(nofLen) = ClusteringData(iTau).Clusters{nofLen}(ach_idx);
                cluster_sim(cluster_num > 0) = 1;  
                TrainingSetClusters(iTau).CondClst(cond_idx).EpochClst(iFullEpoch(cond_idx)).cluster_num(iAvalanche,:) = cluster_num;
                TrainingSetClusters(iTau).CondClst(cond_idx).EpochClst(iFullEpoch(cond_idx)).cluster_sim(iAvalanche,:) = cluster_sim;
            end         
        end
    end
    
    close(f_waitbar);
end

if saveFlg
    save([output_fp fileInfo.orig_fn '_testclust.mat'],'TestingSetClusters','TrainingSetClusters','ClusteringData');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_vector_cluster(v, avch_length_bins, MultiFileAchVecs, ClusteringData, similarity_method, iTau, plotFlg)

nofLen = length(ClusteringData(iTau).Clusters);
cluster_num = zeros(1,nofLen);
cluster_sim = zeros(1,nofLen);

lenInx = [avch_length_bins  nofLen-1  nofLen];
for iLen = lenInx
    
    clusters = unique(ClusteringData(iTau).Clusters{iLen});
    if isempty(clusters)
        continue;
    end
    clusters_similarity = zeros(1,length(clusters));
    
    if plotFlg && iLen ~= nofLen-1 %not full matrix
        plot_hist = zeros(length(clusters),length(v));
    end
    
    for iClust = 1:length(clusters)
        ClusterAchIds = ClusteringData(iTau).Id{iLen}(clusters(iClust) == ClusteringData(iTau).Clusters{iLen});
        if iLen == nofLen %concat
            v_compare = MultiFileAchVecs{str2num(ClusterAchIds{1}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{1}(27:29))}(str2num(ClusterAchIds{1}(33:36))).vec;
            if  avch_length_bins ~= length(v_compare)/MultiFileAchVecs{1}(iTau).nof_channels
                continue;
            end
        end
        
        %fln = str2num(ClusterAchIds{iAch}(21:23));
        %epc = str2num(ClusterAchIds{iAch}(27:29));
        %ach = str2num(ClusterAchIds{iAch}(33:36));

        if strcmp(similarity_method, 'jaccard') && iLen ~= nofLen-1 %not full matrix : faster computation
            Vcompare = [];
            for iAch = 1:length(ClusterAchIds)
                Vcompare = [Vcompare;  MultiFileAchVecs{str2num(ClusterAchIds{iAch}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{iAch}(27:29))}(str2num(ClusterAchIds{iAch}(33:36))).vec];
            end            
            sim_M = squareform(1 - pdist([v; Vcompare],'jaccard'), 'tomatrix');
            %clusters_similarity(iClust) = sum(sim_M(1,:)); %average
            clusters_similarity(iClust) = max(sim_M(1,:)); %nearest neighbour
        else
            for iAch = 1:length(ClusterAchIds)
                v_compare = MultiFileAchVecs{str2num(ClusterAchIds{iAch}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{iAch}(27:29))}(str2num(ClusterAchIds{iAch}(33:36))).vec;
                if iLen == nofLen-1 %full
                    %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method); %average
                    clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method)); %nearest neighbour
                else
                    %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_vectors_similarity(v, v_compare, similarity_method); %average
                    clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_vectors_similarity(v, v_compare, similarity_method)); %nearest neighbour
                    
                    if plotFlg
                        plot_hist(iClust,:) = plot_hist(iClust,:) + v_compare;
                    end
                end
                if clusters_similarity(iClust) == 1 % >0.7 %nearest neighbour
                    %figure;stem([v_compare; v]');legend('nearest cluster vector','vector');
                    break;
                end
            end
        end
        %clusters_similarity(iClust) = clusters_similarity(iClust)/length(ClusterAchIds); %average
        
        if clusters_similarity(iClust) == 1 %for iClust
            break;
        end
    end %for iClust
    
    [cluster_sim(iLen),max_inx] = max(clusters_similarity);
    cluster_num(iLen) = clusters(max_inx);
    
    if plotFlg && iLen ~= nofLen-1 %not full matrix
        figure;stem([plot_hist(max_inx,:)/max(plot_hist(max_inx,:))/0.5; v]');legend('cluster histogram','vector');title(['cluster: ' num2str(cluster_num(iLen)) '   similarity: ' num2str(cluster_sim(iLen))]);
    end
    
end
