%
%Step 8 : predict_conditions - accumulative Beyesian predictor
%
%  inputs:
% TestClusterSets - clusters for test
% ClusteringDataSets - clusters and statistics
% save_str - add this sting to filename
% saveFlg
% plotFlg
%
%  outputs:
% PredictionResultSets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PredictionResultSets = predict_conditions(TestClusterSets, ClusteringDataSets, save_str, saveFlg, plotFlg)

params_t = global_params();

optimal_tau_t = TestClusterSets{1}(~cellfun(@isempty,{TestClusterSets{1}.is_optimal_tau}));
fileInfo = optimal_tau_t(([TestClusterSets{1}.is_optimal_tau] == 1)).fileInfo;

if saveFlg
    output_fp = [fileInfo.base_fp '5 testing\'];
    mkdir(output_fp);
end

debug_len_to_plot = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PredictionResultSets = [];

percent_waitbar = 0;
f_waitbar = waitbar(percent_waitbar, ['predict conditions ' num2str(100*percent_waitbar) '%'], 'Name',fileInfo.orig_fn );   

for iCrossValid = 1:length(TestClusterSets)
    for iTau = 1:length(TestClusterSets{iCrossValid})
        if isempty(TestClusterSets{iCrossValid}(iTau).tau)
            continue;
        end
        
        percent_waitbar = iCrossValid/length(TestClusterSets);
        waitbar(percent_waitbar,f_waitbar,['cross-valid ' num2str(iCrossValid) ' predict conditions  ' num2str(100*percent_waitbar) '%']);  
        
        PredictionResultSets{iCrossValid}(iTau).CondPred = [];
        PredictionResultSets{iCrossValid}(iTau).CondIds = TestClusterSets{iCrossValid}(iTau).CondIds;
        PredictionResultSets{iCrossValid}(iTau).tau = TestClusterSets{iCrossValid}(iTau).tau;
        PredictionResultSets{iCrossValid}(iTau).is_optimal_tau = TestClusterSets{iCrossValid}(iTau).is_optimal_tau;
        PredictionResultSets{iCrossValid}(iTau).fileInfo = fileInfo;
        
        nof_cond = length(TestClusterSets{iCrossValid}(iTau).CondIds);
        for iCond = 1:nof_cond
            
            %init
            for iLen = 1:length(ClusteringDataSets{iCrossValid}(iTau).Clusters)
                for iAccum = 1:length(params_t.accTypes)
                    switch params_t.accTypes{iAccum}
                        case 'ThreshAccum'
                            nof_ths = length(params_t.condition_descision_threshold);
                        case 'SampLimitAccum'
                            nof_ths = length(params_t.condition_counter_limit);
                        otherwise
                            nof_ths = 1;
                    end
                    for iThs = 1:nof_ths
                        PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum})(iCond,iThs) =  predictor_init_next(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, []);
                    end
                end
            end
            
            %test epochs clusters
            for iEpoch = 1:length(TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst)
                nof_epoch_vecs = size(TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_num,1);
                for iVec = 1:nof_epoch_vecs
                    avch_length_bins = find(TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec,:));
                    avch_length_bins(TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_sim(iVec,avch_length_bins) <= params_t.minimal_similarity_threshold) = [];
                    for iLen = avch_length_bins
                        for iAccum = 1:length(params_t.accTypes)
                            cond_predictors = PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum});
                            for iThs = 1:size(cond_predictors,2)
                                cond_predictors(iCond,iThs) = predictor_accumulate(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec,iLen),...
                                    TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_sim(iVec,iLen), cond_predictors(iCond,iThs));
                                %next ThreshAccum, SampLimitAccum
                                if strcmp(params_t.accTypes{iAccum},'ThreshAccum')
                                    cond_predictors(iCond,iThs) = predictor_threshold_decide(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, cond_predictors(iCond,iThs), params_t.condition_descision_threshold(iThs), params_t.check_prediction_salience, nof_epoch_vecs);
                                end
                                if strcmp(params_t.accTypes{iAccum},'SampLimitAccum')
                                    cond_predictors(iCond,iThs) = predictor_counter_decide(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, cond_predictors(iCond,iThs), params_t.condition_counter_limit(iThs), params_t.check_prediction_salience);
                                end
                            end
                            PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum}) = cond_predictors;
                        end
                    end
                end %for iVec
                
                %next EpochAccum
                if any(contains(params_t.accTypes,'EpochAccum'))
                    for iLen = 1:length(ClusteringDataSets{iCrossValid}(iTau).Clusters)
                        PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.EpochAccum(iCond,1) = predictor_decide_last(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.EpochAccum(iCond,1), params_t.check_prediction_salience, []);
                        if iEpoch < length(TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst)
                            PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.EpochAccum(iCond,1) = predictor_init_next(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.EpochAccum(iCond,1));
                        end
                    end
                end
            end %for iEpoch
            
            %decide TotalAccum
            if any(contains(params_t.accTypes,'TotalAccum'))
                for iLen = 1:length(ClusteringDataSets{iCrossValid}(iTau).Clusters)
                    plot_title = [];
                    if plotFlg
                        if iLen == length(ClusteringDataSets{iCrossValid}(iTau).Clusters)
                            plot_title = ['total, concat, condition: ' num2str(iCond)];
                        elseif iLen == length(ClusteringDataSets{iCrossValid}(iTau).Clusters)-1
                            plot_title = ['total, full, condition: ' num2str(iCond)];
                        elseif any(iLen == debug_len_to_plot)
                            plot_title = ['total, Avalanche Length: ' num2str(iLen,'%d') ' actual condition: ' num2str(iCond)];
                        end
                    end
                    PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.TotalAccum(iCond,1) = predictor_decide_last(ClusteringDataSets{iCrossValid}(iTau).Stats{iLen}, PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.TotalAccum(iCond,1), params_t.check_prediction_salience, plot_title);
                end
            end
            
        end %for iCond
        
        %ROC & display
        for iLen = 1:length(ClusteringDataSets{iCrossValid}(iTau).Clusters)
            for iAccum = 1:length(params_t.accTypes)
                PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum}) = calc_predictor_performance(PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum}), ClusteringDataSets{iCrossValid}(iTau).Stats{iLen});
                if plotFlg
                    if iLen == length(ClusteringDataSets{iCrossValid}(iTau).Clusters)
                        diplay_predictor_results(PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum}), PredictionResultSets{iCrossValid}(iTau).CondIds, params_t, ['concat, ' params_t.accTypes{iAccum}]);
                    elseif iLen == length(ClusteringDataSets{iCrossValid}(iTau).Clusters)-1
                        diplay_predictor_results(PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum}), PredictionResultSets{iCrossValid}(iTau).CondIds, params_t, ['full, ' params_t.accTypes{iAccum}]);
                    elseif any(iLen == debug_len_to_plot)
                        diplay_predictor_results(PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum}), PredictionResultSets{iCrossValid}(iTau).CondIds, params_t, [params_t.accTypes{iAccum}, ' Avalanche Length: ' num2str(iLen,'%d')]);
                    end
                end
            end
        end
        
    end %for iTau
end %for iCrossValid

close(f_waitbar);

if saveFlg
    save([output_fp fileInfo.orig_fn '_' save_str 'Predict.mat'],'PredictionResultSets','TestClusterSets','ClusteringDataSets');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_init_next(Stats, predictor)

if isempty(predictor)
    predictor.conditions_prob_log10 = [];
    predictor.ach_cnt = [];
    predictor.decision_cond = [];
end

predictor.conditions_prob_log10(:,end+1) = log10(Stats.P_cond)';
predictor.ach_cnt(end+1) = 0;
predictor.decision_cond(end+1) = 0;
predictor.conditions_prob_log10_history = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_accumulate(Stats, cluster_num, cluster_sim, predictor)

cluster_sim = 1;
% %test predictor
% cluster_sim = 1; nof_clust = 4;
% % TestClusterSets{iCrossValid}(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec) = iCond; % cluster_num = ceil(4*rand); % cluster_num = 3;
% Stats.P_cond = ones(1,4)/4; % Stats.P_cond = [0.00001 0.99999 0.00001 0.00001];
% Stats.P_clst = ones(1,nof_clust)/nof_clust; % Stats.P_clst = [0.05 0 0.95 0];
% Stats.P_clstGINVcond = ones(4,nof_clust)/nof_clust; %Stats.P_clstGINVcond = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% Stats.P_clondGINVclst = ones(4,nof_clust)/4;

minimal_p_clst = 0.001; %to avoid zero probability and -Inf log

% log_p = log10(Stats.P_clstGINVcond(:,cluster_num)/Stats.P_clst(cluster_num));
log_p = log10(Stats.P_clstGINVcond(:,cluster_num)/(Stats.P_cond * Stats.P_clstGINVcond(:,cluster_num))); %use full bayes because P_cond based on nof_epochs but P_clst based on nof_avalanches
log_p(log_p == -Inf) = log10(minimal_p_clst);
log_p(log_p == Inf | isnan(log_p)) = 0; %log_p == -Inf is allowed
predictor.conditions_prob_log10(:,end) = predictor.conditions_prob_log10(:,end) + cluster_sim*log_p;
predictor.ach_cnt(end) = predictor.ach_cnt(end) + 1;
predictor.conditions_prob_log10_history = [predictor.conditions_prob_log10_history  predictor.conditions_prob_log10(:,end)];

% %average-based salience
% step_accum = diff(conditions_prob_log10_history);
% step_av = median(abs(reshape(diff(conditions_prob_log10_history),1,[]))); %mean % median
% conditions_prob_log10_sorted = sort(predictor.conditions_prob_log10(:,end),'descend');
% decision_salience = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/step_av;
% %contrast-based salience
% decision_salience = (conditions_prob_log10_sorted(1)-conditions_prob_log10_sorted(2))/(abs(conditions_prob_log10_sorted(1))+abs(conditions_prob_log10_sorted(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_threshold_decide(Stats, predictor, condition_descision_threshold, check_prediction_salience, max_cnt)

if max(predictor.conditions_prob_log10(:,end)) >= condition_descision_threshold % * (1-predictor.ach_cnt(end)/max_cnt)   %collapsing threshold
    plot_title = [];%['thresh = ' num2str(condition_descision_threshold,'%1.2f')];
    predictor = predictor_decide_last(Stats, predictor, check_prediction_salience, plot_title);
    if predictor.decision_cond(end) > 0 %descision was made
        predictor = predictor_init_next(Stats, predictor);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_counter_decide(Stats, predictor, condition_counter_limit, check_prediction_salience)

if predictor.ach_cnt(end) >= condition_counter_limit
    predictor = predictor_decide_last(Stats, predictor, check_prediction_salience, []);
    predictor = predictor_init_next(Stats, predictor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predictor = predictor_decide_last(Stats, predictor, check_prediction_salience, plot_title)

[max_prob, max_inx] = max(predictor.conditions_prob_log10(:,end));    
if predictor.ach_cnt(end) > 0 && max_prob > log10(Stats.P_cond(max_inx)) %above chance
    if ~check_prediction_salience || is_salient(predictor)
        predictor.decision_cond(end) = max_inx;
    end
end
if ~isempty(plot_title)
    figure;plot(predictor.conditions_prob_log10_history');xlabel('step');title(['Accumulators - ' plot_title ' predicted condition: ' num2str(predictor.decision_cond(end))]);
    legend('cond 1','cond 2','cond 3','cond 4', 'cond 5', 'cond 6');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isSalient = is_salient(predictor)

mean_diff_thresh = 0; %should not be negative

isSalient = true;           
[~, max_inx] = max(predictor.conditions_prob_log10(:,end));
for i=1:size(predictor.conditions_prob_log10_history,1)
    if i ~= max_inx
        history_diff = predictor.conditions_prob_log10_history(max_inx,:)-predictor.conditions_prob_log10_history(i,:);
        history_diff_mean = mean(history_diff, 2);
        if history_diff_mean <= mean_diff_thresh
            isSalient = false;
            return;
        end
        last_negative_inx = find(history_diff<=0, 1, 'last'); %find last non-growing point
        if ~isempty(last_negative_inx) && history_diff(end) <= max(abs(history_diff(1:last_negative_inx))) %check if diff > max peak before growing
            isSalient = false;
            return;
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cond_predictors = calc_predictor_performance(cond_predictors, Stats)

nof_cond = size(cond_predictors,1);
for iThs = 1:size(cond_predictors,2)
    %count detections
    p_det_m = []; %detection matrix
    tp_av_cnt_total = []; fa_av_cnt_total = [];
    for iCond = 1:nof_cond %can't be 0
        nof_det = [];
        for iCondDet = 1:nof_cond %can't be 0
            nof_det(iCondDet) = sum(cond_predictors(iCond,iThs).decision_cond == iCondDet);
        end
        cond_predictors(iCond,iThs).roc.p_det = max(0,nof_det/sum(cond_predictors(iCond,iThs).decision_cond > 0));
        if any(isnan(cond_predictors(iCond,iThs).roc.p_det))
            error('p_det  NAN');
        end
        p_det_m = [p_det_m; cond_predictors(iCond,iThs).roc.p_det];
        
        cond_predictors(iCond,iThs).roc.tp_av_cnt = mean(cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  == iCond & cond_predictors(iCond,iThs).decision_cond > 0));
        tp_av_cnt_total = [tp_av_cnt_total  cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  == iCond & cond_predictors(iCond,iThs).decision_cond > 0)];
        cond_predictors(iCond,iThs).roc.fa_av_cnt = mean(cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  ~= iCond & cond_predictors(iCond,iThs).decision_cond > 0));
        fa_av_cnt_total = [fa_av_cnt_total  cond_predictors(iCond,iThs).ach_cnt(cond_predictors(iCond,iThs).decision_cond  ~= iCond & cond_predictors(iCond,iThs).decision_cond > 0)];
    end
    %calc true-positive, false-alarm
    tp_total = Stats.P_cond * diag(p_det_m);
    fa_total = 1 - tp_total; % sum(Stats.P_cond * p_det_m) - tp_total; % can be calculated this way as well
    tp_av_cnt_total = mean(tp_av_cnt_total);
    fa_av_cnt_total = mean(fa_av_cnt_total);
    for iCond = 1:nof_cond
        cond_predictors(iCond,iThs).roc.tp = p_det_m(iCond,iCond);
        normalized_p_cond = Stats.P_cond(1:nof_cond~=iCond)/sum(Stats.P_cond(1:nof_cond~=iCond)); %same as divide by (1-Stats.P_cond(iCond))
        cond_predictors(iCond,iThs).roc.fa = normalized_p_cond * p_det_m(1:nof_cond ~= iCond,iCond);
        %fa_total = fa_total + cond_predictors(iCond,iThs).roc.fa *(1-Stats.P_cond(iCond)); % can be calculated this way as well
        
        cond_predictors(iCond,iThs).roc.tp_total = tp_total;
        cond_predictors(iCond,iThs).roc.fa_total = fa_total;
        cond_predictors(iCond,iThs).roc.tp_av_cnt_total = tp_av_cnt_total;
        cond_predictors(iCond,iThs).roc.fa_av_cnt_total = fa_av_cnt_total;
        cond_predictors(iCond,iThs).roc.p_det_m = p_det_m;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diplay_predictor_results(cond_predictors, CondIds, params_t, disp_str)

roc_fields = {'tp','fa','tp_av_cnt','fa_av_cnt'};
total_roc_fields = {'tp_total','fa_total','tp_av_cnt_total','fa_av_cnt_total'};

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
    ax = [ax subplot(1,length(CondIds)+1,(length(CondIds)+1)*0+iCond+1)];scatter(fa(iCond,:),tp(iCond,:),'m*');xlabel('fa');ylabel('tp');title([CondIds(iCond) ' - ROC']);xlim([-0.1 1.1]);ylim([-0.1 1.1]);grid on;grid minor;
%     ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*2+iCond+1)];scatter(fa_av_cnt(iCond,:),fa_av_salience(iCond,:),'m*');xlabel('cnt');ylabel('salience');title([CondIds(iCond) ' - FA Counter VS Salience']);grid on;grid minor;
end
ax = [ax subplot(1,length(CondIds)+1,(length(CondIds)+1)*0+1)];scatter(fa_total,tp_total,'m*');xlabel('fa');ylabel('tp');title([disp_str ' - Total ROC']);xlim([-0.1 1.1]);ylim([-0.1 1.1]);grid on; grid minor;
% ax = [ax subplot(3,length(CondIds)+1,(length(CondIds)+1)*2+1)];scatter(fa_av_cnt_total,fa_av_salience_total,'m*');xlabel('cnt');ylabel('salience');title('Total FA Counter VS Salience');grid on;grid minor;

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
%         ' tp_cnt=' num2str(cond_predictors(iCond,1).roc.tp_av_cnt,'%1.2f') ...
%         ' fa_cnt=' num2str(cond_predictors(iCond,1).roc.fa_av_cnt,'%1.2f') '   '];
% end
% disp([disp_str ...
%     '  TP=' num2str(cond_predictors(1,1).roc.tp_total,'%1.2f')...
%     '  TP_cnt=' num2str(cond_predictors(1,1).roc.tp_av_cnt_total,'%1.2f')...
%     '  FA_cnt=' num2str(cond_predictors(1,1).roc.fa_av_cnt_total,'%1.2f')  '  -    ' disp_str_roc]);
disp([disp_str '  tested (row) VS detected (columns) conditions probability:']);
disp(num2str(cond_predictors(1,1).roc.p_det_m));
