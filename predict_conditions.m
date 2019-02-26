%
%Step 5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predict_conditions()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select clustering results files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

for iFile = 1:length(files)
    load([fp files{iFile}],'ClusteringData');
    load([fp [files{iFile}(1:end-13) '.mat']],'MultiFileAchVecs','SimilarityMat','TestingSet');
    
    PredictionResults = [];
    for iTau = 1:length(ClusteringData)
        
        if isempty(ClusteringData(iTau).tau)
            continue;
        end
        
        PredictionResults(iTau).CondIds = TestingSet(iTau).CondIds;
        for iCond=1:length(TestingSet(iTau).CondIds)
            PredictionResults(iTau).TotalAccum(iCond).conditions_prob_log10 = NaN;
            PredictionResults(iTau).TotalAccum(iCond).ach_cnt = 0;
            PredictionResults(iTau).ThreshAccum(iCond).conditions_prob_log10 = NaN;
            PredictionResults(iTau).ThreshAccum(iCond).ach_cnt = 0;   
            PredictionResults(iTau).EpochAccum(iCond).conditions_prob_log10 = NaN;
            PredictionResults(iTau).EpochAccum(iCond).ach_cnt = 0;             
            
            for iEpoch=1:length(TestingSet(iTau).EpochIds{iCond})

                
                %fln = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
                %epc = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
                ach_vectors_t = MultiFileAchVecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6))}(iTau).epochs_vecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12))};
                
                for iVec=1:length(ach_vectors_t)
                    avch_length_bins = length(ach_vectors_t(iVec).vec)/MultiFileAchVecs{1}(iTau).nof_channels;
                    [cluster_num, cluster_sim] = find_vector_cluster(ach_vectors_t(iVec).vec, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, params_t.similarity_method, 'full', iTau);
                    if cluster_sim >= params_t.minimal_similarity_threshold
                        PredictionResults(iTau).TotalAccum(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).StatsFull, cluster_num, cluster_sim, PredictionResults(iTau).TotalAccum(iCond).conditions_prob_log10(:,end));
                        PredictionResults(iTau).TotalAccum(iCond).ach_cnt(end) = PredictionResults(iTau).TotalAccum(iCond).ach_cnt(end) + 1;
                        PredictionResults(iTau).ThreshAccum(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).StatsFull, cluster_num, cluster_sim, PredictionResults(iTau).ThreshAccum(iCond).conditions_prob_log10(:,end));
                        PredictionResults(iTau).ThreshAccum(iCond).ach_cnt(end) = PredictionResults(iTau).ThreshAccum(iCond).ach_cnt(end) + 1;
                        PredictionResults(iTau).EpochAccum(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).StatsFull, cluster_num, cluster_sim, PredictionResults(iTau).EpochAccum(iCond).conditions_prob_log10(:,end));
                        PredictionResults(iTau).EpochAccum(iCond).ach_cnt(end) = PredictionResults(iTau).EpochAccum(iCond).ach_cnt(end) + 1;                        
                    end
                    
                    
%                     [cluster_num, cluster_sim] = find_vector_cluster(ach_vectors_t(iVec).vec, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, params_t.similarity_method, 'len', iTau);
%                     conditions_prob_log10 = bayessian_accumulative_predictor(ClusteringData(iTau).StatsLen{avch_length_bins}, cluster_num, conditions_prob_log10, accum_reset_flg));
%                     [cluster_num, cluster_sim] = find_vector_cluster(ach_vectors_t(iVec).vec, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, params_t.similarity_method, 'concat', iTau);
%                     conditions_prob_log10 = bayessian_accumulative_predictor(ClusteringData(iTau).StatsConcatLen, cluster_num, conditions_prob_log10, accum_reset_flg));                    

                    if sum(PredictionResults(iTau).ThreshAccum(iCond).conditions_prob_log10(:,end) >= params_t.condition_descision_threshold)>0
                        PredictionResults(iTau).ThreshAccum(iCond).conditions_prob_log10(:,end+1) = NaN;
                        PredictionResults(iTau).ThreshAccum(iCond).ach_cnt(end+1) = 0;
                    end
                    
                end
                
                
                PredictionResults(iTau).EpochAccum(iCond).conditions_prob_log10(:,end+1) = NaN;
                PredictionResults(iTau).EpochAccum(iCond).ach_cnt(end+1) = 0;

                
            end
        end
        
        
    end
    
    save([fp files{iFile}(1:end-4) '_prediction.mat'],'PredictionResults');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conditions_prob_log10 = bayessian_accumulative_predictor(Stats, cluster_num, cluster_sim, conditions_prob_log10)

%exponential primacy weighting?
%static thresholds? distribution depandent?
%dynamic thresholds? accumulator dynamic range and std depandent?
%threshold collapses over time? as function of nof samples?

if isnan(conditions_prob_log10(1)) %init after reset 
    conditions_prob_log10 = log10(Stats.P_cond);
end
conditions_prob_log10 = conditions_prob_log10 + log10(cluster_sim*Stats.P_clstGINVcond(:,cluster_num)/Stats.P_clst(cluster_num));
% conditions_prob_log10 = conditions_prob_log10 + log10(cluster_sim*Stats.P_clstGINVcond(:,cluster_num));

% %test these:
% conditions_prob = Stats.P_clondGINVclst(:,cluster_num); %straight forward
% conditions_prob = Stats.P_clstGINVcond(:,cluster_num).*Stats.P_cond/Stats.P_clst(cluster_num); %bayes
% conditions_prob = Stats.P_clstGINVcond(:,cluster_num).*Stats.P_cond/sum(Stats.P_clstGINVcond(:,cluster_num).*Stats.P_cond); %full bayes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_vector_cluster(v, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, similarity_method, compare_cluster_type, iTau)

switch compare_cluster_type
    case 'full'
        clusters_v = ClusteringData(iTau).ClustersFull;
        Id_list = SimilarityMat(iTau).Id;
    case 'len'
        clusters_v = ClusteringData(iTau).ClustersLen{avch_length_bins};
        Id_list = SimilarityMat(iTau).IdLen{avch_length_bins};
    case 'concat'
        clusters_v = ClusteringData(iTau).ClustersConcatLen;
        Id_list = ClusteringData(iTau).IdLConcatLen;
end

clusters = unique(clusters_v);
clusters_similarity = zeros(1,length(clusters));
for iClust = 1:length(clusters)
    ClusterAchIds = Id_list(clusters(iClust) == clusters_v);
    if strcmp(compare_cluster_type, 'concat')
        v_compare = MultiFileAchVecs{str2num(ClusterAchIds{1}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{1}(27:29))}(str2num(ClusterAchIds{1}(33:36))).vec;
        if  avch_length_bins ~= length(v_compare)/MultiFileAchVecs{1}(iTau).nof_channels
            continue;
        end
    end
    for iAch = 1:length(ClusterAchIds)
        %fln = str2num(ClusterAchIds{iAch}(21:23));
        %epc = str2num(ClusterAchIds{iAch}(27:29));
        %ach = str2num(ClusterAchIds{iAch}(33:36));
        v_compare = MultiFileAchVecs{str2num(ClusterAchIds{iAch}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{iAch}(27:29))}(str2num(ClusterAchIds{iAch}(33:36))).vec;
        if strcmp(compare_cluster_type, 'full')
            %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method); %average
            clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method)); %nearest neighbour
        else
            %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_vectors_similarity(v, v_compare, similarity_method); %average
            clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_vectors_similarity(v, v_compare, similarity_method)); %nearest neighbour
        end
    end
    %clusters_similarity(iClust) = clusters_similarity(iClust)/length(ClusterAchIds); %average
end
[cluster_sim,max_inx] = max(clusters_similarity);
cluster_num = clusters(max_inx);
