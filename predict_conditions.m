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
                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond) = predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), []);
                    end
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        for iAccum = 1:length(accumulator_types)
                            PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond) =  predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){iLen}, []);
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
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond) = predictor_accumulate(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]),...
                                        TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec), TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec),...
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond));
                                    %next ThreshAccum, SampLimitAccum
                                    if strcmp(accumulator_types{iAccum},'ThreshAccum')
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond) = predictor_threshold_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond), params_t, 10);
                                    end
                                    if strcmp(accumulator_types{iAccum},'SampLimitAccum')
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond) = predictor_counter_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond), params_t, 10);
                                    end                                    
                                end
                            else
                                for iAccum = 1:length(accumulator_types)
                                    PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond) = predictor_accumulate(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins},...
                                        TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec), TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec),...
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond));
                                    %next ThreshAccum, SampLimitAccum
                                    if strcmp(accumulator_types{iAccum},'ThreshAccum')
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond) = predictor_threshold_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins}, PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond), params_t, 10);
                                    end
                                    if strcmp(accumulator_types{iAccum},'SampLimitAccum')
                                        PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond) = predictor_counter_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins}, PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum})(iCond), params_t, 10);
                                    end   
                                end                                
                            end
                        end   
                    end %for iVec
                     
                    %next EpochAccum
                    if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                        PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond));
                        if iEpoch < length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                            PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond) = predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond));
                        end
                    else
                        for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                            PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond));
                            if iEpoch < length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                                PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond) = predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins}, PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond));
                            end
                        end
                    end
                     
                end %for iEpoch
                
                %decide TotalAccum
                if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                    PredictionResults(iTau).(params_t.compare_modes{iMode}).TotalAccum(iCond) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}).TotalAccum(iCond));
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.TotalAccum(iCond) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.TotalAccum(iCond));
                    end
                end

            end %for iCond
        end %for iMode
        
        %ROC & display
        for iMode = 1:length(params_t.compare_modes)
            disp(['----------- Mode:  ' params_t.compare_modes{iMode} ' -------------------']);
            if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                for iAccum = 1:length(accumulator_types)
                    PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}) = calc_predictor_performance(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}),...
                        ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]));                    
                    disp_str_roc = [];
                    for iCond = 1:length(PredictionResults(iTau).CondIds)
                        disp_str_roc = [disp_str_roc PredictionResults(iTau).CondIds{iCond}...
                            ': tp=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.tp,'%1.2f')...
                            ' fa=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.fa,'%1.2f') ...
                            ' tp_sal=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.tp_av_salience,'%1.2f') ...
                            ' fa_sal=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.fa_av_salience,'%1.2f') ...
                            ' tp_cnt=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.tp_av_cnt,'%1.2f') ...
                            ' fa_cnt=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond).roc.fa_av_cnt,'%1.2f') '   '];
                    end
                    disp([accumulator_types{iAccum}...
                        '  TP=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.tp_total,'%1.2f')...
                        '  TP_sal=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.tp_av_salience_total,'%1.2f')...
                        '  FA_sal=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.fa_av_salience_total,'%1.2f')...
                        '  TP_cnt=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.tp_av_cnt_total,'%1.2f')...
                        '  FA_cnt=' num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.fa_av_cnt_total,'%1.2f')  '  -    ' disp_str_roc]); 
                    disp('tested (row) VS detected (columns) conditions probability:'); 
                    disp(num2str(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(1).roc.p_det_m));
                end
            else
                for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                    for iAccum = 1:length(accumulator_types)
                        PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}) = calc_predictor_performance(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}),...
                            ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins});
                    end 
                end
            end
        end
        
    end %for iTau
    
    save([fp files{iFile}(1:end-14) '_testpred.mat'],'PredictionResults');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_init_next(Stats, predictor)

if isempty(predictor)
    predictor.conditions_prob_log10 = [];
    predictor.ach_cnt = [];
    predictor.decision_cond = [];
    predictor.decision_salience = [];
end

predictor.conditions_prob_log10(:,end+1) = log10(Stats.P_cond)';
predictor.ach_cnt(end+1) = 0;
predictor.decision_cond(end+1) = 0;
predictor.decision_salience(end+1) = NaN;
predictor.step_accum = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_accumulate(Stats, cluster_num, cluster_sim, predictor)

% %test predictor
% cluster_sim = 1; nof_clust = 4;
% % TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec) = iCond; % cluster_num = ceil(4*rand); % cluster_num = 3;
% Stats.P_cond = ones(1,4)/4; % Stats.P_cond = [0.00001 0.99999 0.00001 0.00001];
% Stats.P_clst = ones(1,nof_clust)/nof_clust; % Stats.P_clst = [0.05 0 0.95 0];
% Stats.P_clstGINVcond = ones(4,nof_clust)/nof_clust; %Stats.P_clstGINVcond = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% Stats.P_clondGINVclst = ones(4,nof_clust)/4;

minimal_p_clst = 0.001; %to avoid zero probability and -Inf log

% log_p = log10(Stats.P_clstGINVcond(:,cluster_num)/Stats.P_clst(cluster_num));
log_p = log10(Stats.P_clstGINVcond(:,cluster_num)/(Stats.P_cond * Stats.P_clstGINVcond(:,cluster_num))); %use full bayes because P_cond based on nof_epochs but P_clst based on nof_avalanches
log_p(log_p == -Inf) = log10(minimal_p_clst);
log_p(log_p == Inf | isnan(log_p)) = 0; %log_p == -Inf is allowed
% cluster_sim = (cluster_sim-params_t.minimal_similarity_threshold)/(1-params_t.minimal_similarity_threshold);% normalize?
predictor.conditions_prob_log10(:,end) = predictor.conditions_prob_log10(:,end) + cluster_sim*log_p;
predictor.ach_cnt(end) = predictor.ach_cnt(end) + 1;

predictor.step_accum = [predictor.step_accum cluster_sim*log_p];
step_av = median(abs(reshape(predictor.step_accum,1,[])));
% step_av = mean(abs(reshape(predictor.step_accum,1,[])));
% step_std = std(abs(predictor.conditions_prob_log10(:,end))/predictor.ach_cnt(end));
conditions_prob_log10_sorted = sort(predictor.conditions_prob_log10(:,end),'descend');
if step_av > 0 % step_std > 0 %
% if abs(conditions_prob_log10_sorted(1))+abs(conditions_prob_log10_sorted(2)) >0
    predictor.decision_salience(end) = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/step_av;
    %predictor.decision_salience(end) = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/step_std;
    %predictor.decision_salience(end) = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/(abs(conditions_prob_log10_sorted(1))+abs(conditions_prob_log10_sorted(2)));
else
    predictor.decision_salience(end) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_threshold_decide(Stats, predictor, params_t, iThs)

above_chance_thresh = log10(1/size(predictor.conditions_prob_log10,1));

if predictor.ach_cnt(end) > 0 && max(predictor.conditions_prob_log10(:,end)) > above_chance_thresh && ...
        predictor.decision_salience(end) >= params_t.condition_descision_step_av_threshold*(1-predictor.ach_cnt(end)/params_t.max_cnt)%collapsing threshold
        %predictor.decision_salience(end) >= params_t.condition_descision_step_std_threshold%*(1-predictor.ach_cnt(end)/params_t.max_cnt)%collapsing threshold
        %predictor.decision_salience(end) >= params_t.condition_descision_contrast_threshold%*(1-predictor.ach_cnt(end)/params_t.max_cnt)%collapsing threshold
        %max(predictor.conditions_prob_log10(:,end)) > params_t.condition_descision_static_threshold(iThs)
    predictor = predictor_decide_last(predictor);
    predictor = predictor_init_next(Stats, predictor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_counter_decide(Stats, predictor, params_t, iCnt)

if predictor.ach_cnt(end) >= params_t.condition_descision_static_counter_limit(iCnt)
    predictor = predictor_decide_last(predictor);
    predictor = predictor_init_next(Stats, predictor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_decide_last(predictor)
if predictor.ach_cnt(end) > 0
    [~,predictor.decision_cond(end)] = max(predictor.conditions_prob_log10(:,end));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cond_predictors = calc_predictor_performance(cond_predictors, Stats)
    
    nof_cond = length(cond_predictors);

    %count detections
    p_det_m = []; %detection matrix
    tp_av_salience_total = []; fa_av_salience_total = []; tp_av_cnt_total = []; fa_av_cnt_total = [];
    for iCond = 1:nof_cond %can't be 0
        nof_det = [];
        for iCondDet = 1:nof_cond %can't be 0
            nof_det(iCondDet) = sum(cond_predictors(iCond).decision_cond == iCondDet);
        end
        cond_predictors(iCond).roc.p_det = nof_det/sum(cond_predictors(iCond).decision_cond > 0);
        p_det_m = [p_det_m; cond_predictors(iCond).roc.p_det];
        
        cond_predictors(iCond).roc.tp_av_salience = mean(cond_predictors(iCond).decision_salience(cond_predictors(iCond).decision_cond == iCond & cond_predictors(iCond).decision_cond > 0),'omitnan');
        tp_av_salience_total = [tp_av_salience_total  cond_predictors(iCond).decision_salience(cond_predictors(iCond).decision_cond == iCond & cond_predictors(iCond).decision_cond > 0)];
        cond_predictors(iCond).roc.fa_av_salience = mean(cond_predictors(iCond).decision_salience(cond_predictors(iCond).decision_cond ~= iCond & cond_predictors(iCond).decision_cond > 0),'omitnan');
        fa_av_salience_total = [fa_av_salience_total  cond_predictors(iCond).decision_salience(cond_predictors(iCond).decision_cond ~= iCond & cond_predictors(iCond).decision_cond > 0)];
        cond_predictors(iCond).roc.tp_av_cnt = mean(cond_predictors(iCond).ach_cnt(cond_predictors(iCond).decision_cond  == iCond & cond_predictors(iCond).decision_cond > 0));
        tp_av_cnt_total = [tp_av_cnt_total  cond_predictors(iCond).ach_cnt(cond_predictors(iCond).decision_cond  == iCond & cond_predictors(iCond).decision_cond > 0)];
        cond_predictors(iCond).roc.fa_av_cnt = mean(cond_predictors(iCond).ach_cnt(cond_predictors(iCond).decision_cond  ~= iCond & cond_predictors(iCond).decision_cond > 0));
        fa_av_cnt_total = [fa_av_cnt_total  cond_predictors(iCond).ach_cnt(cond_predictors(iCond).decision_cond  ~= iCond & cond_predictors(iCond).decision_cond > 0)];
    end
    %calc true-positive, false-alarm
    tp_total = Stats.P_cond * diag(p_det_m);
    fa_total = 1 - tp_total; % sum(Stats.P_cond * p_det_m) - tp_total; % can be calculated this way as well
    tp_av_salience_total = mean(tp_av_salience_total,'omitnan');
    fa_av_salience_total = mean(fa_av_salience_total,'omitnan');
    tp_av_cnt_total = mean(tp_av_cnt_total); 
    fa_av_cnt_total = mean(fa_av_cnt_total);
    cond_idx = 1:nof_cond;
    for iCond = 1:nof_cond
        cond_predictors(iCond).roc.tp = p_det_m(iCond,iCond);
        normalized_p_cond = Stats.P_cond(cond_idx(cond_idx~=iCond))/sum(Stats.P_cond(cond_idx(cond_idx~=iCond))); %same as divide by (1-Stats.P_cond(iCond))
        cond_predictors(iCond).roc.fa = normalized_p_cond * p_det_m(cond_idx(cond_idx~=iCond),iCond);
        
        %fa_total = fa_total + cond_predictors(iCond).roc.fa *(1-Stats.P_cond(iCond)); % can be calculated this way as well
        
        cond_predictors(iCond).roc.tp_total = tp_total;
        cond_predictors(iCond).roc.fa_total = fa_total;
        cond_predictors(iCond).roc.tp_av_salience_total = tp_av_salience_total;
        cond_predictors(iCond).roc.fa_av_salience_total = fa_av_salience_total;
        cond_predictors(iCond).roc.tp_av_cnt_total = tp_av_cnt_total;
        cond_predictors(iCond).roc.fa_av_cnt_total = fa_av_cnt_total;
        cond_predictors(iCond).roc.p_det_m = p_det_m;         
    end
