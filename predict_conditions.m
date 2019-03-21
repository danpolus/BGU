%
%Step 6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predict_conditions()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();

accumulator_types = {'TotalAccum', 'ThreshAccum', 'EpochAccum', 'SampLimitAccum'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select testing clusters files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

for iFile = 1:length(files)
    load([fp files{iFile}],'TestingSetClusters','ClusteringData');
    
    PredictionResults = [];
    
    for iTau = 1:length(TestingSetClusters)
        if isempty(TestingSetClusters(iTau).tau)
            continue;
        end
        
        PredictionResults(iTau).tau = TestingSetClusters(iTau).tau;
        PredictionResults(iTau).is_optimal_tau = TestingSetClusters(iTau).is_optimal_tau;
        PredictionResults(iTau).CondIds = TestingSetClusters(iTau).CondIds;
        nof_cond = length(TestingSetClusters(iTau).CondIds);
        
        for iMode = 1:length(params_t.compare_modes)
            for iCond = 1:nof_cond
                
                %init
                if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                    for iAccum = 1:length(accumulator_types)
                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10 = nan(nof_cond,1);
                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt = 0;
                    end
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        for iAccum = 1:length(accumulator_types)
                            PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10 = nan(nof_cond,1);
                            PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).ach_cnt = 0;
                        end
                    end
                end

                %test epochs clusters
                for iEpoch = 1:length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                    for iVec = 1:length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num)
                        avch_length_bins = TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).avch_length_bins(iVec);

                        if TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec) >= params_t.minimal_similarity_threshold
                            if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                                for iAccum = 1:length(accumulator_types)
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]),...
                                        TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec) , TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec),...
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end));
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) = PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) + 1;
                                end
                            else
                                for iAccum = 1:length(accumulator_types)
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = bayessian_accumulative_predictor(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins},...
                                        TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec) , TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec),...
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end));
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).ach_cnt(end) = PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).ach_cnt(end) + 1;
                                end                                
                            end
                        end
                        
                        %next ThreshAccum, SampLimitAccum
                        for iAccum = 1:length(accumulator_types)
                            if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                                if (strcmp(accumulator_types{iAccum},'ThreshAccum') && sum(PredictionResults(iTau).(params_t.compare_modes{iMode}).ThreshAccum(iCond).conditions_prob_log10(:,end) >= params_t.condition_descision_threshold) > 0) || ...
                                        (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(params_t.compare_modes{iMode}).SampLimitAccum(iCond).ach_cnt(end) >= params_t.accum_sample_limit)
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end+1) = 0;
                                end
                            else
                                if (strcmp(accumulator_types{iAccum},'ThreshAccum') && sum(PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.ThreshAccum(iCond).conditions_prob_log10(:,end) >= params_t.condition_descision_threshold) > 0) || ...
                                        (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.SampLimitAccum(iCond).ach_cnt(end) >= params_t.accum_sample_limit)
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond).ach_cnt(end+1) = 0;
                                end
                            end
                        end
                          
                    end %for iVec
                     
                    %next EpochAccum
                    if iEpoch < length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                        if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                            PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                            PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond).ach_cnt(end+1) = 0;
                        else
                            for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                                PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond).conditions_prob_log10(:,end+1) = nan(nof_cond,1);
                                PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond).ach_cnt(:,end+1) = 0;
                            end
                        end
                    end
                     
                end %for iEpoch
                
                %remove last element if it's not based on suficient data
                if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                    for iAccum = 1:length(accumulator_types)
                        if size(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10, 2) > 1
                            max_p = max(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end),[],1);
                            if (strcmp(accumulator_types{iAccum},'ThreshAccum') && (max_p<params_t.condition_descision_threshold || isnan(max_p))) ||...
                                    (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) < params_t.accum_sample_limit)
                                PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = [];
                                PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).ach_cnt(end) = [];
                            end
                        end
                    end
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        for iAccum = 1:length(accumulator_types)
                            if size(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10, 2) > 1
                                max_p = max(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end),[],1);
                                if ((strcmp(accumulator_types{iAccum},'ThreshAccum') && (max_p<params_t.condition_descision_threshold || isnan(max_p))) ||...
                                        (strcmp(accumulator_types{iAccum},'SampLimitAccum') && PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).ach_cnt(end)  < params_t.accum_sample_limit))
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).conditions_prob_log10(:,end) = [];
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond).ach_cnt(end) = [];
                                end
                            end
                        end
                    end
                end
                
            end %for iCond
        end %for iMode
        
        %ROC & display
        for iMode = 1:length(params_t.compare_modes)
            disp(['----------- Mode:  ' params_t.compare_modes{iMode} ' -------------------']);
            if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                for iAccum = 1:length(accumulator_types)
                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}) = calc_detector_performance(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}),...
                        ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]));                    
                    disp_str_roc = [];
                    for iCond = 1:length(PredictionResults(iTau).CondIds)
                        disp_str_roc = [disp_str_roc PredictionResults(iTau).CondIds{iCond} ': tp=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.tp,'%1.2f')...
                            ' fa=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.fa,'%1.2f') '   '];
                    end
                    disp([accumulator_types{iAccum} '  TP=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.tp_total,'%1.2f') '  -    ' disp_str_roc]); 
                    disp('tested (row) VS detected (columns) conditions probability:'); 
                    disp(num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.p_det_m));
                end
            else
                for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                    for iAccum = 1:length(accumulator_types)
                        PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}) = calc_detector_performance(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}),...
                            ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins});
                    end 
                end
            end
        end
        
    end %for iTau
    
    save([fp files{iFile}(1:end-14) '_testpred.mat'],'PredictionResults');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conditions_prob_log10 = bayessian_accumulative_predictor(Stats, cluster_num, cluster_sim, conditions_prob_log10)

% %test predictor
% cluster_sim = 1; nof_clust = 4;
% % TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec) = iCond; % cluster_num = ceil(4*rand); % cluster_num = 3;
% Stats.P_cond = ones(1,4)/4; % Stats.P_cond = [0.00001 0.99999 0.00001 0.00001];
% Stats.P_clst = ones(1,nof_clust)/nof_clust; % Stats.P_clst = [0.05 0 0.95 0];
% Stats.P_clstGINVcond = ones(4,nof_clust)/nof_clust; %Stats.P_clstGINVcond = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% Stats.P_clondGINVclst = ones(4,nof_clust)/4;

if isnan(conditions_prob_log10(1)) %init after reset
    conditions_prob_log10 = log10(Stats.P_cond)';
end
% log_p = log10(Stats.P_clstGINVcond(:,cluster_num)/Stats.P_clst(cluster_num));
log_p = log10(Stats.P_clstGINVcond(:,cluster_num)/(Stats.P_cond * Stats.P_clstGINVcond(:,cluster_num))); %use full bayes because P_cond based on nof_epochs but P_clst based on nof_avalanches
log_p(isinf(log_p) | isnan(log_p)) = 0;
% cluster_sim = (cluster_sim-params_t.minimal_similarity_threshold)/(1-params_t.minimal_similarity_threshold);% normalize?
conditions_prob_log10 = conditions_prob_log10 + cluster_sim*log_p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function accum_t_v = calc_detector_performance(accum_t_v, Stats)
    
    nof_cond = length(accum_t_v);

    %count detections
    p_det_m = []; %detection matrix
    for iCond = 1:nof_cond
        [max_p, max_i] = max(accum_t_v(iCond).conditions_prob_log10,[],1);
        for iCondDet = 1:nof_cond
            accum_t_v(iCond).roc.nof_det(iCondDet) = sum(max_i==iCondDet & ~isnan(max_p));
        end
        accum_t_v(iCond).roc.p_det = accum_t_v(iCond).roc.nof_det/size(accum_t_v(iCond).conditions_prob_log10,2);
        p_det_m = [p_det_m; accum_t_v(iCond).roc.p_det];
    end
    %calc true-positive, false-alarm
    tp_total = Stats.P_cond * diag(p_det_m);
    fa_total = 1 - tp_total; % sum(Stats.P_cond * p_det_m) - tp_total; % can be calculated this way as well
    cond_idx = 1:nof_cond;
    for iCond = 1:nof_cond
        accum_t_v(iCond).roc.tp = p_det_m(iCond,iCond);
        normalized_p_cond = Stats.P_cond(cond_idx(cond_idx~=iCond))/sum(Stats.P_cond(cond_idx(cond_idx~=iCond))); %same as divide by (1-Stats.P_cond(iCond))
        accum_t_v(iCond).roc.fa = normalized_p_cond * p_det_m(cond_idx(cond_idx~=iCond),iCond);
        
        %fa_total = fa_total + accum_t_v(iCond).roc.fa *(1-Stats.P_cond(iCond)); % can be calculated this way as well
        
        accum_t_v(iCond).roc.tp_total = tp_total;
        accum_t_v(iCond).roc.fa_total = fa_total;
        accum_t_v(iCond).roc.p_det_m = p_det_m;         
    end
