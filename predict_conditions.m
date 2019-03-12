%
%Step 5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predict_conditions()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();

compare_modes = {'Full', 'Len', 'ConcatLen'};%{'Len'};
accumulator_types = {'TotalAccum', 'ThreshAccum', 'EpochAccum', 'SampLimitAccum'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select clustering results files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

for iFile = 1:length(files)
    load([fp files{iFile}],'ClusteringData');
    load([fp [files{iFile}(1:end-13) '.mat']],'MultiFileAchVecs','SimilarityMat','TestingSet');
    
    tic
    
    PredictionResults = [];
    
    for iTau = 1:length(ClusteringData)
        if isempty(ClusteringData(iTau).tau)
            continue;
        end
        
%         %to display ROC, put breakpoint after this command
%         load([fp files{iFile}(1:end-4) '_prediction.mat'],'PredictionResults');
        
        PredictionResults(iTau).CondIds = TestingSet(iTau).CondIds;
        nof_cond = length(TestingSet(iTau).CondIds);
        
        for iMode = 1:length(compare_modes)
            for iCond = 1:nof_cond
                
                %init
                if ~strcmp(compare_modes{iMode}, 'Len')
                    for iAccum = 1:length(accumulator_types)
                        PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10 = nan(nof_cond,1);
                        PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt = 0;
                    end
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        for iAccum = 1:length(accumulator_types)
                            PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10 = nan(nof_cond,1);
                            PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).ach_cnt = 0;
                        end
                    end
                end

                %test epochs vectors
                for iEpoch = 1:length(TestingSet(iTau).EpochIds{iCond})
                    %fln = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
                    %epc = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
                    ach_vectors_t = MultiFileAchVecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6))}(iTau).epochs_vecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12))};
                    
                    for iVec = 1:length(ach_vectors_t)
                        avch_length_bins = length(ach_vectors_t(iVec).vec)/MultiFileAchVecs{1}(iTau).nof_channels;
                        [cluster_num, cluster_sim] = find_vector_cluster(ach_vectors_t(iVec).vec, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, params_t.similarity_method, compare_modes{iMode}, iTau);
                        if cluster_sim >= params_t.minimal_similarity_threshold
                            if ~strcmp(compare_modes{iMode}, 'Len')
                                for iAccum = 1:length(accumulator_types)
                                    PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).(['Stats' compare_modes{iMode}]),...
                                        cluster_num, cluster_sim, PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end));
                                    PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) = PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) + 1;
                                end
                            else
                                for iAccum = 1:length(accumulator_types)
                                    PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).(['Stats' compare_modes{iMode}]){avch_length_bins},...
                                        cluster_num, cluster_sim, PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end));
                                    PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).ach_cnt(end) = PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).ach_cnt(end) + 1;
                                end                                
                            end
                        end
                        
                        %next ThreshAccum, SampLimitAccum
                        for iAccum = 1:length(accumulator_types)
                            if ~strcmp(compare_modes{iMode}, 'Len')
                                if (strcmp(accumulator_types{iAccum},'ThreshAccum') && sum(PredictionResults(iTau).(compare_modes{iMode}).ThreshAccum(iCond).conditions_prob_log10(:,end) >= params_t.condition_descision_threshold) > 0) || ...
                                        (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(compare_modes{iMode}).SampLimitAccum(iCond).ach_cnt(end) >= params_t.accum_sample_limit)
                                    PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                                    PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end+1) = 0;
                                end
                            else
                                if (strcmp(accumulator_types{iAccum},'ThreshAccum') && sum(PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.ThreshAccum(iCond).conditions_prob_log10(:,end) >= params_t.condition_descision_threshold) > 0) || ...
                                        (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.SampLimitAccum(iCond).ach_cnt(end) >= params_t.accum_sample_limit)
                                    PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                                    PredictionResults(iTau).(compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).ach_cnt(end+1) = 0;
                                end
                            end
                        end
                          
                    end %for iVec
                     
                    %next EpochAccum
                    if iEpoch < length(TestingSet(iTau).EpochIds{iCond})
                        if ~strcmp(compare_modes{iMode}, 'Len')
                            PredictionResults(iTau).(compare_modes{iMode}).EpochAccum(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                            PredictionResults(iTau).(compare_modes{iMode}).EpochAccum(iCond).ach_cnt(end+1) = 0;
                        else
                            for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                                PredictionResults(iTau).(compare_modes{iMode}){iLen}.EpochAccum(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                                PredictionResults(iTau).(compare_modes{iMode}){iLen}.EpochAccum(iCond).ach_cnt(:,end+1) = 0;
                            end
                        end
                    end
                     
                end %for iEpoch
                
                %remove last element if it's not based on suficient data
                if ~strcmp(compare_modes{iMode}, 'Len')
                    for iAccum = 1:length(accumulator_types)
                        if size(PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10, 2) > 1
                            max_p = max(PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end),[],1);
                            if (strcmp(accumulator_types{iAccum},'ThreshAccum') && (max_p<params_t.condition_descision_threshold || isnan(max_p))) ||...
                                    (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) < params_t.accum_sample_limit)
                                PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = [];
                                PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) = [];
                            end
                        end
                    end
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        for iAccum = 1:length(accumulator_types)
                            if size(PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10, 2) > 1
                                max_p = max(PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end),[],1);
                                if ((strcmp(accumulator_types{iAccum},'ThreshAccum') && (max_p<params_t.condition_descision_threshold || isnan(max_p))) ||...
                                        (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).ach_cnt(end)  < params_t.accum_sample_limit))
                                    PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = [];
                                    PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).ach_cnt(end) = [];
                                end
                            end
                        end
                    end
                end
                
            end %for iCond
        end %for iMode
        
        %ROC & display
        for iMode = 1:length(compare_modes)
            disp(['----------- Mode:  ' compare_modes{iMode} ' -------------------']);
            disp(['conditions: ' PredictionResults(iTau).CondIds]);
            if ~strcmp(compare_modes{iMode}, 'Len')
                for iAccum = 1:length(accumulator_types)
                    PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum}) = calc_detector_performance(PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum}));
                    disp([accumulator_types{iAccum} ':   TP = ' num2str(PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.tp_total) ...
                        ' FA = ' num2str(PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.fa_total) ' detection matrix (tested VS detected) = ']); 
                      disp(num2str(PredictionResults(iTau).(compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.p_det_m));
                end
            else
                for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                    for iAccum = 1:length(accumulator_types)
                        PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum}) = calc_detector_performance(PredictionResults(iTau).(compare_modes{iMode}){iLen}.(accumulator_types{iAccum}));
                    end 
                end
            end
        end
        
    end %for iTau
    
    toc
    
    save([fp files{iFile}(1:end-4) '_prediction.mat'],'PredictionResults');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conditions_prob_log10 = bayessian_accumulative_predictor(Stats, cluster_num, cluster_sim, conditions_prob_log10)

%exponential primacy weighting? other weighting? w1*log(P(Ci|A1))?+ w2*log(P(Ci|A2))? ...+ ...wk*log(P(Ci|Ak))
%static thresholds? distribution depandent?
%dynamic thresholds? accumulator dynamic range and std depandent?
%threshold collapses over time? as function of nof samples?
if cluster_sim <= 0
    return;
end

if isnan(conditions_prob_log10(1)) %init after reset
    conditions_prob_log10 = log10(Stats.P_cond)';
    conditions_prob_log10(isinf(conditions_prob_log10)) = 0;
end
if size(Stats.P_clstGINVcond,2) >= cluster_num && size(Stats.P_clst,2) >= cluster_num
    log_p = log10(cluster_sim*Stats.P_clstGINVcond(:,cluster_num)/Stats.P_clst(cluster_num));
    log_p(isinf(log_p)) = 0;
    conditions_prob_log10 = conditions_prob_log10 + log_p;
end
% conditions_prob_log10 = conditions_prob_log10 + log10(cluster_sim*Stats.P_clstGINVcond(:,cluster_num)).*(Stats.P_clstGINVcond(:,cluster_num) > 0);

% %test these:
% if ~isempty(Stats.P_clondGINVclst)
% conditions_prob = Stats.P_clondGINVclst(:,cluster_num); %straight forward
% conditions_prob = Stats.P_clstGINVcond(:,cluster_num).*Stats.P_cond/Stats.P_clst(cluster_num); %bayes
% conditions_prob = Stats.P_clstGINVcond(:,cluster_num).*Stats.P_cond/sum(Stats.P_clstGINVcond(:,cluster_num).*Stats.P_cond); %full bayes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_vector_cluster(v, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, similarity_method, compare_mode, iTau)

switch compare_mode
    case 'Full'
        clusters_v = ClusteringData(iTau).ClustersFull;
        Id_list = SimilarityMat(iTau).Id;
    case 'Len'
        clusters_v = ClusteringData(iTau).ClustersLen{avch_length_bins};
        Id_list = SimilarityMat(iTau).IdLen{avch_length_bins};
    case 'ConcatLen'
        clusters_v = ClusteringData(iTau).ClustersConcatLen;
        Id_list = ClusteringData(iTau).IdLConcatLen;
end

clusters = unique(clusters_v);
clusters_similarity = zeros(1,length(clusters));
for iClust = 1:length(clusters)
    ClusterAchIds = Id_list(clusters(iClust) == clusters_v);
    if strcmp(compare_mode, 'ConcatLen')
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
        if strcmp(compare_mode, 'Full')
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function accum_t_v = calc_detector_performance(accum_t_v)
    
    nof_cond = length(accum_t_v);

    %count detections
    for iCond = 1:nof_cond
        [max_p, max_i] = max(accum_t_v(iCond).conditions_prob_log10,[],1);
        for iCondDet = 1:nof_cond
            accum_t_v(iCond).roc.nof_det(iCondDet) = sum(max_i==iCondDet & ~isnan(max_p));
        end
        accum_t_v(iCond).roc.p_det = accum_t_v(iCond).roc.nof_det/size(accum_t_v(iCond).conditions_prob_log10,2);
    end
    %calc true-positive, false-alarm
    tp_total = 0; fa_total = 0;
    p_det_m = [];
    for iCond = 1:nof_cond
        accum_t_v(iCond).roc.tp = accum_t_v(iCond).roc.p_det(iCond);
        accum_t_v(iCond).roc.fa = 0;
        fa_sum = 0; fa_cnt = 0;
        for iCondDet = 1:nof_cond
            if iCondDet~= iCond
                fa_sum = fa_sum + accum_t_v(iCondDet).roc.nof_det(iCond);
                fa_cnt = fa_cnt + size(accum_t_v(iCondDet).conditions_prob_log10,2);
            end
        end
        if fa_cnt > 0
            accum_t_v(iCond).roc.fa = fa_sum/fa_cnt;
        end
        tp_total = tp_total + accum_t_v(iCond).roc.tp; 
        fa_total = fa_total + accum_t_v(iCond).roc.fa;
        p_det_m = [p_det_m; accum_t_v(iCond).roc.p_det];
    end
    tp_total = tp_total/nof_cond;
    fa_total = fa_total/nof_cond;
    for iCond = 1:nof_cond
        accum_t_v(iCond).roc.tp_total = tp_total;
        accum_t_v(iCond).roc.fa_total = fa_total;
        accum_t_v(iCond).roc.p_det_m = p_det_m;
    end
