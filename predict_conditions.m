%
%Step 7
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predict_conditions()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();

accumulator_types = {'TotalAccum', 'ThreshAccum', 'EpochAccum', 'SampLimitAccum'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fn, fp] = uigetfile([fp '*.mat'], 'Select testing clusters file');

load([fp fn],'TestingSetClusters','ClusteringData');

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
                    switch accumulator_types{iAccum}
                        case 'ThreshAccum'
                            nof_ths = length(params_t.condition_descision_threshold);
                        case 'SampLimitAccum'
                            nof_ths = length(params_t.condition_counter_limit);
                        otherwise
                            nof_ths = 1;
                    end
                    for iThs = 1:nof_ths
                        PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum})(iCond,iThs) = predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), []);
                    end
                end
            else
                for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                    for iAccum = 1:length(accumulator_types)
                        switch accumulator_types{iAccum}
                            case 'ThreshAccum'
                                nof_ths = length(params_t.condition_descision_threshold);
                            case 'SampLimitAccum'
                                nof_ths = length(params_t.condition_counter_limit);
                            otherwise
                                nof_ths = 1;
                        end
                        for iThs = 1:nof_ths
                            PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum})(iCond,iThs) =  predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){iLen}, []);
                        end
                    end
                end
            end

            %test epochs clusters
            for iEpoch = 1:length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                for iVec = 1:length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num)
                    avch_length_bins = TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).avch_length_bins(iVec);
                    max_cnt = length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num); %params_t.max_cnt

                    if TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec) >= params_t.minimal_similarity_threshold
                        if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                            for iAccum = 1:length(accumulator_types)
                                cond_predictors = PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum});
                                for iThs = 1:size(cond_predictors,2)
                                    cond_predictors(iCond,iThs) = predictor_accumulate(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec),...
                                        TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec), cond_predictors(iCond,iThs));
                                    %next ThreshAccum, SampLimitAccum
                                    if strcmp(accumulator_types{iAccum},'ThreshAccum')
                                        cond_predictors(iCond,iThs) = predictor_threshold_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), cond_predictors(iCond,iThs), params_t.condition_descision_threshold(iThs), max_cnt);
                                    end
                                    if strcmp(accumulator_types{iAccum},'SampLimitAccum')
                                        cond_predictors(iCond,iThs) = predictor_counter_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), cond_predictors(iCond,iThs), params_t.condition_counter_limit(iThs));
                                    end
                                end
                                PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}) = cond_predictors;
                            end
                        else
                            for iAccum = 1:length(accumulator_types)
                                cond_predictors = PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum});
                                for iThs = 1:size(cond_predictors,2)
                                    cond_predictors(iCond,iThs) = predictor_accumulate(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins}, TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec),...
                                        TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec), cond_predictors(iCond,iThs));
                                    %next ThreshAccum, SampLimitAccum
                                    if strcmp(accumulator_types{iAccum},'ThreshAccum')
                                        cond_predictors(iCond,iThs) = predictor_threshold_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins}, cond_predictors(iCond,iThs), params_t.condition_descision_threshold(iThs), max_cnt);
                                    end
                                    if strcmp(accumulator_types{iAccum},'SampLimitAccum')
                                        cond_predictors(iCond,iThs) = predictor_counter_decide(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){avch_length_bins}, cond_predictors(iCond,iThs), params_t.condition_counter_limit(iThs));
                                    end
                                end
                                PredictionResults(iTau).(params_t.compare_modes{iMode}){avch_length_bins}.(accumulator_types{iAccum}) = cond_predictors;
                            end
                        end
                    end
                end %for iVec

                %next EpochAccum
                if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                    PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond,1) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond,1),[]);
                    if iEpoch < length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                        PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond,1) = predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]), PredictionResults(iTau).(params_t.compare_modes{iMode}).EpochAccum(iCond,1));
                    end
                else
                    for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                        PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond,1) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond,1),[]);
                        if iEpoch < length(TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst)
                            PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond,1) = predictor_init_next(ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){iLen}, PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.EpochAccum(iCond,1));
                        end
                    end
                end

            end %for iEpoch

            %decide TotalAccum
            if ~strcmp(params_t.compare_modes{iMode}, 'Len')
                plot_title = ['total, mode: ' params_t.compare_modes{iMode} ', condition: ' num2str(iCond)];
                PredictionResults(iTau).(params_t.compare_modes{iMode}).TotalAccum(iCond,1) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}).TotalAccum(iCond,1),plot_title);
            else
                for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                    PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.TotalAccum(iCond,1) = predictor_decide_last(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.TotalAccum(iCond,1),[]);
                end
            end

        end %for iCond
    end %for iMode

    %ROC & display
    for iMode = 1:length(params_t.compare_modes)
        if ~strcmp(params_t.compare_modes{iMode}, 'Len')
            for iAccum = 1:length(accumulator_types)
                PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}) = calc_predictor_performance(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}),...
                    ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]));
                disp_str = ['Mode: ' params_t.compare_modes{iMode} ', ' accumulator_types{iAccum}];
                diplay_predictor_results(PredictionResults(iTau).(params_t.compare_modes{iMode}).(accumulator_types{iAccum}), PredictionResults(iTau).CondIds, params_t,disp_str);
            end
        else
            for iLen = 1:length(ClusteringData(iTau).ClustersLen)
                for iAccum = 1:length(accumulator_types)
                    PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}) = calc_predictor_performance(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}),...
                        ClusteringData(iTau).(['Stats' params_t.compare_modes{iMode}]){iLen});
                    %disp_str = ['Mode: ' params_t.compare_modes{iMode} ', ' accumulator_types{iAccum} 'Avalanche Length: ' num2str(iLen,'%d')];
                    %diplay_predictor_results(PredictionResults(iTau).(params_t.compare_modes{iMode}){iLen}.(accumulator_types{iAccum}), PredictionResults(iTau).CondIds, params_t, disp_str);
                end
            end
        end
    end

end %for iTau

save([fp fn(1:end-4) '_pred.mat'],'PredictionResults');
    



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
conditions_prob_log10_sorted = sort(predictor.conditions_prob_log10(:,end),'descend');

% %average-based salience
% step_av = median(abs(reshape(predictor.step_accum,1,[]))); %mean % median
% if step_av > 0
%     predictor.decision_salience(end) = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/step_av;

%contrast-based salience
if abs(conditions_prob_log10_sorted(1))+abs(conditions_prob_log10_sorted(2)) >0
    predictor.decision_salience(end) = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/(abs(conditions_prob_log10_sorted(1))+abs(conditions_prob_log10_sorted(2)));

% %std-based salience: poor results    
% step_std = std(abs(predictor.conditions_prob_log10(:,end))/predictor.ach_cnt(end));
% if step_std > 0
%     predictor.decision_salience(end) = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/step_std;

else
    predictor.decision_salience(end) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_threshold_decide(Stats, predictor, condition_descision_threshold, max_cnt)

% max_cnt = 20;
above_chance_thresh = log10(1/size(predictor.conditions_prob_log10,1));

if predictor.ach_cnt(end) > 0 && max(predictor.conditions_prob_log10(:,end)) > above_chance_thresh
%     if predictor.decision_salience(end) >= condition_descision_threshold   %average-based salience or std-based salience
    if predictor.decision_salience(end) >= condition_descision_threshold * (1-predictor.ach_cnt(end)/max_cnt)   %contrast-based salience. collapsing threshold due to saturation at high cnt 
%     if max(predictor.conditions_prob_log10(:,end)) > condition_descision_threshold   %static threshold
        plot_title = [];%['thresh = ' num2str(condition_descision_threshold * (1-predictor.ach_cnt(end)/max_cnt),'%1.2f')];
        predictor = predictor_decide_last(predictor,plot_title);
        predictor = predictor_init_next(Stats, predictor);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_counter_decide(Stats, predictor, condition_counter_limit)

if predictor.ach_cnt(end) >= condition_counter_limit
    predictor = predictor_decide_last(predictor,[]);
    predictor = predictor_init_next(Stats, predictor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_decide_last(predictor,plot_title)
if predictor.ach_cnt(end) > 0
    [~,predictor.decision_cond(end)] = max(predictor.conditions_prob_log10(:,end));
end
if ~isempty(plot_title)
    figure;plot(predictor.step_accum');xlabel('step');title(['Accumulators - ' plot_title]);legend('cond 1','cond 2','cond 3','cond 4');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cond_predictors = calc_predictor_performance(cond_predictors, Stats)

nof_cond = size(cond_predictors,1);
for iThs = 1:size(cond_predictors,2)
    %count detections
    p_det_m = []; %detection matrix
    tp_av_salience_total = []; fa_av_salience_total = []; tp_av_cnt_total = []; fa_av_cnt_total = [];
    for iCond = 1:nof_cond %can't be 0
        nof_det = [];
        for iCondDet = 1:nof_cond %can't be 0
            nof_det(iCondDet) = sum(cond_predictors(iCond,iThs).decision_cond == iCondDet);
        end
        cond_predictors(iCond,iThs).roc.p_det = nof_det/sum(cond_predictors(iCond,iThs).decision_cond > 0);
        p_det_m = [p_det_m; cond_predictors(iCond,iThs).roc.p_det];
        
        cond_predictors(iCond,iThs).roc.tp_av_salience = mean(cond_predictors(iCond,iThs).decision_salience(cond_predictors(iCond,iThs).decision_cond == iCond & cond_predictors(iCond,iThs).decision_cond > 0),'omitnan');
        tp_av_salience_total = [tp_av_salience_total  cond_predictors(iCond,iThs).decision_salience(cond_predictors(iCond,iThs).decision_cond == iCond & cond_predictors(iCond,iThs).decision_cond > 0)];
        cond_predictors(iCond,iThs).roc.fa_av_salience = mean(cond_predictors(iCond,iThs).decision_salience(cond_predictors(iCond,iThs).decision_cond ~= iCond & cond_predictors(iCond,iThs).decision_cond > 0),'omitnan');
        fa_av_salience_total = [fa_av_salience_total  cond_predictors(iCond,iThs).decision_salience(cond_predictors(iCond,iThs).decision_cond ~= iCond & cond_predictors(iCond,iThs).decision_cond > 0)];
        cond_predictors(iCond,iThs).roc.tp_av_cnt = mean(cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  == iCond & cond_predictors(iCond,iThs).decision_cond > 0));
        tp_av_cnt_total = [tp_av_cnt_total  cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  == iCond & cond_predictors(iCond,iThs).decision_cond > 0)];
        cond_predictors(iCond,iThs).roc.fa_av_cnt = mean(cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  ~= iCond & cond_predictors(iCond,iThs).decision_cond > 0));
        fa_av_cnt_total = [fa_av_cnt_total  cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  ~= iCond & cond_predictors(iCond,iThs).decision_cond > 0)];
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
        cond_predictors(iCond,iThs).roc.tp = p_det_m(iCond,iCond);
        normalized_p_cond = Stats.P_cond(cond_idx(cond_idx~=iCond))/sum(Stats.P_cond(cond_idx(cond_idx~=iCond))); %same as divide by (1-Stats.P_cond(iCond,iThs))
        cond_predictors(iCond,iThs).roc.fa = normalized_p_cond * p_det_m(iCond,cond_idx(cond_idx~=iCond))';
        
        %fa_total = fa_total + cond_predictors(iCond,iThs).roc.fa *(1-Stats.P_cond(iCond,iThs)); % can be calculated this way as well
        
        cond_predictors(iCond,iThs).roc.tp_total = tp_total;
        cond_predictors(iCond,iThs).roc.fa_total = fa_total;
        cond_predictors(iCond,iThs).roc.tp_av_salience_total = tp_av_salience_total;
        cond_predictors(iCond,iThs).roc.fa_av_salience_total = fa_av_salience_total;
        cond_predictors(iCond,iThs).roc.tp_av_cnt_total = tp_av_cnt_total;
        cond_predictors(iCond,iThs).roc.fa_av_cnt_total = fa_av_cnt_total;
        cond_predictors(iCond,iThs).roc.p_det_m = p_det_m;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diplay_predictor_results(cond_predictors, CondIds, params_t, disp_str)

roc_fields = {'tp','fa','tp_av_salience','fa_av_salience','tp_av_cnt','fa_av_cnt'};
total_roc_fields = {'tp_total','fa_total','tp_av_salience_total','fa_av_salience_total','tp_av_cnt_total','fa_av_cnt_total'};

for iRoc = 1:length(roc_fields)
    for iCond = 1:length(CondIds)
        for iThs = 1:size(cond_predictors,2)
            eval([roc_fields{iRoc} '(iCond,iThs) = cond_predictors(iCond,iThs).roc.' roc_fields{iRoc} ';']);
        end
    end
end
for iRoc = 1:length(total_roc_fields)
    for iThs = 1:size(cond_predictors,2)
        eval([total_roc_fields{iRoc} '(iThs) = cond_predictors(1,iThs).roc.' total_roc_fields{iRoc} ';']);
    end
end

ax = [];    
figure('Name',disp_str);
for iCond = 1:length(CondIds)
    ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*0+iCond+1)];scatter(fa(iCond,:),tp(iCond,:),'x');xlabel('fa');ylabel('tp');title([CondIds(iCond) ' - ROC']);
    ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*1+iCond+1)];scatter(tp_av_cnt(iCond,:),tp_av_salience(iCond,:),'x');xlabel('cnt');ylabel('salience');title([CondIds(iCond) ' - TP Counter VS Salience']);
    ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*2+iCond+1)];scatter(fa_av_cnt(iCond,:),fa_av_salience(iCond,:),'x');xlabel('cnt');ylabel('salience');title([CondIds(iCond) ' - FA Counter VS Salience']);    
end
ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*0+1)];scatter(fa_total,tp_total,'x');xlabel('fa');ylabel('tp');title([disp_str ' - Total ROC']);
ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*1+1)];scatter(tp_av_cnt_total,tp_av_salience_total,'x');xlabel('cnt');ylabel('salience');title(['Total TP Counter VS Salience']);
ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*2+1)];scatter(fa_av_cnt_total,fa_av_salience_total,'x');xlabel('cnt');ylabel('salience');title(['Total FA Counter VS Salience']);

for iAx = 1:length(ax)
    x=get(get(ax(iAx),'children'),'Xdata');
    y=get(get(ax(iAx),'children'),'Ydata');
    if contains(disp_str, 'ThreshAccum')
        text(ax(iAx),x,y,num2str(params_t.condition_descision_threshold','%1.2f'));
    end
    if contains(disp_str, 'SampLimitAccum')
        text(ax(iAx),x,y,num2str(params_t.condition_counter_limit','%1.0f'));
    end
end

% disp_str_roc = disp_str;
% for iCond = 1:length(CondIds)
%     disp_str_roc = [disp_str_roc CondIds{iCond}...
%         ': tp=' num2str(cond_predictors(iCond,1).roc.tp,'%1.2f')...
%         ' fa=' num2str(cond_predictors(iCond,1).roc.fa,'%1.2f') ...
%         ' tp_sal=' num2str(cond_predictors(iCond,1).roc.tp_av_salience,'%1.2f') ...
%         ' fa_sal=' num2str(cond_predictors(iCond,1).roc.fa_av_salience,'%1.2f') ...
%         ' tp_cnt=' num2str(cond_predictors(iCond,1).roc.tp_av_cnt,'%1.2f') ...
%         ' fa_cnt=' num2str(cond_predictors(iCond,1).roc.fa_av_cnt,'%1.2f') '   '];
% end
% disp([disp_str ...
%     '  TP=' num2str(cond_predictors(1,1).roc.tp_total,'%1.2f')...
%     '  TP_sal=' num2str(cond_predictors(1,1).roc.tp_av_salience_total,'%1.2f')...
%     '  FA_sal=' num2str(cond_predictors(1,1).roc.fa_av_salience_total,'%1.2f')...
%     '  TP_cnt=' num2str(cond_predictors(1,1).roc.tp_av_cnt_total,'%1.2f')...
%     '  FA_cnt=' num2str(cond_predictors(1,1).roc.fa_av_cnt_total,'%1.2f')  '  -    ' disp_str_roc]);
disp([disp_str '  tested (row) VS detected (columns) conditions probability:']);
disp(num2str(cond_predictors(1,1).roc.p_det_m));
