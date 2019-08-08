clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();

thresh_inx = round(length(params_t.condition_descision_threshold)/2);
cntlmt_inx = round(length(params_t.condition_counter_limit)/2);

debug_len_to_plot = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select results data files', 'MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

plotLenInx = [debug_len_to_plot  max([0 debug_len_to_plot])+1  max([0 debug_len_to_plot])+2];

allSubjectsResults = [];
P_cond = [];
CondIds = [];
for iFile = 1:length(files)
    
    load([fp files{iFile}],'PredictionResultSets','ClusteringDataSets');

    for iCrossValid = 1:length(PredictionResultSets)
        for iTau = 1:length(PredictionResultSets{iCrossValid})
            if isempty(PredictionResultSets{iCrossValid}(iTau).tau)
                continue;
            end 
            
            P_cond(iFile,iCrossValid,:) = ClusteringDataSets{iCrossValid}(iTau).Stats{1}.P_cond;
            if isempty(CondIds)
                for iCond = 1:length(PredictionResultSets{1}(iTau).CondIds)
                    CondIds{iCond} = PredictionResultSets{1}(iTau).CondIds{iCond}(end-7:end);
                end
            end
            
            for iLen = 1:length(PredictionResultSets{iCrossValid}(iTau).CondPred)
                if any(iLen == debug_len_to_plot)
                    plotLen = iLen;
                elseif iLen == length(PredictionResultSets{iCrossValid}(iTau).CondPred)-1
                    plotLen = plotLenInx(end-1);
                elseif iLen == length(PredictionResultSets{iCrossValid}(iTau).CondPred)
                    plotLen = plotLenInx(end);
                else
                    continue;
                end
                for iAccum = 1:length(params_t.accTypes)
                    switch params_t.accTypes{iAccum}
                        case 'ThreshAccum'
                            nof_ths = thresh_inx;
                        case 'SampLimitAccum'
                            nof_ths = cntlmt_inx;
                        otherwise
                            nof_ths = 1;
                    end
                    roc = PredictionResultSets{iCrossValid}(iTau).CondPred{iLen}.(params_t.accTypes{iAccum})(1,nof_ths).roc;
                    allSubjectsResults(plotLen,iAccum).tpConditions(iFile,iCrossValid,:) = diag(roc.p_det_m);
                    %roc.tp_total = squeeze(P_cond(iFile,iCrossValid,:))' * diag(roc.p_det_m);
                end
            end
            
        end
    end
    
end

%%%%%%%

nof_files = length(files);
nof_cond = length(CondIds);

for iLen = plotLenInx
    for iAccum = 1:length(params_t.accTypes)
        if iLen < max(plotLenInx)-1
            len_str = num2str(iLen);
        elseif iLen == max(plotLenInx)-1
            len_str = 'full';
        elseif iLen == max(plotLenInx)
            len_str = 'concat';
        end
        figure('Name',['results per subject:  Len = ' len_str '  ' params_t.accTypes{iAccum}]);
        
        tpConditionsAllFileMat = [];
        P_condAllFileMat = [];
        for iFile = 1:nof_files
            tpValidMat = squeeze(allSubjectsResults(iLen,iAccum).tpConditions(iFile,:,:));
            P_condValidMat = squeeze(P_cond(iFile,:,:));
            kappaMat = 1 - (1 - tpValidMat) ./ (1 - P_condValidMat);
            if size(tpValidMat,2) ~= nof_cond
                tpValidMat = tpValidMat';
                P_condValidMat = P_condValidMat';
                kappaMat = kappaMat';
            end
            tpConditionsAllFileMat = [tpConditionsAllFileMat; tpValidMat];
            P_condAllFileMat = [P_condAllFileMat; P_condValidMat];
        
            crosss_valid_conditions_median = median(tpValidMat,1);
            crosss_valid_conditions_std = std(tpValidMat,0,1);
            crosss_valid_conditions_kappa = median(kappaMat, 1);
            crosss_valid_conditions_kappa_std = std(kappaMat, 0, 1);

            crosss_valid_tp_median(iFile) = median(diag(P_condValidMat*tpValidMat'));
            mu = mean(diag(P_condValidMat*tpValidMat'));
            crosss_valid_tp_std(iFile) = sqrt(mean(diag(P_condValidMat*(tpValidMat' - mu).^2)) * numel(tpValidMat)/(numel(tpValidMat)-1)); %unbiased
            crosss_valid_tp_kappa_median(iFile) = median(diag(P_condValidMat*kappaMat'));
            mu = mean(diag(P_condValidMat*kappaMat'));
            crosss_valid_tp_kappa_std(iFile) = sqrt(mean(diag(P_condValidMat*(kappaMat' - mu).^2)) * numel(kappaMat)/(numel(kappaMat)-1)); %unbiased
            
            subplot(2,nof_files,iFile);bar(1:nof_cond,crosss_valid_conditions_median);hold on;
            errorbar(1:nof_cond,crosss_valid_conditions_median,crosss_valid_conditions_std,'LineStyle','none');hold off;xticklabels(CondIds);title([files{iFile}(1:end-4) ': cross-validation median']);
            subplot(2,nof_files,nof_files+iFile);bar(1:nof_cond,crosss_valid_conditions_kappa);hold on;
            errorbar(1:nof_cond,crosss_valid_conditions_kappa,crosss_valid_conditions_kappa_std,'LineStyle','none');xticklabels(CondIds);title('cross-validation median Kappa');
        end
        
        crosss_subj_conditions_median = median(tpConditionsAllFileMat,1);
        crosss_subj_conditions_std = std(tpConditionsAllFileMat,0,1);
        crosss_subj_conditions_kappa = median(1 - (1 - tpConditionsAllFileMat) ./ (1 - P_condAllFileMat), 1);
        crosss_subj_conditions_kappa_std = std(1 - (1 - tpConditionsAllFileMat) ./ (1 - P_condAllFileMat),0,1);
            
        figure('Name',['inter-condition results:  Len = ' len_str '  ' params_t.accTypes{iAccum}]);
        subplot(1,2,1);bar(1:nof_cond,crosss_subj_conditions_median);hold on;
        errorbar(1:nof_cond,crosss_subj_conditions_median,crosss_subj_conditions_std,'LineStyle','none');hold off;xticklabels(CondIds);title('inter-subject median');
        subplot(1,2,2);bar(1:nof_cond,crosss_subj_conditions_kappa);hold on;
        errorbar(1:nof_cond,crosss_subj_conditions_kappa,crosss_subj_conditions_kappa_std,'LineStyle','none');xticklabels(CondIds);title('inter-subject median Kappa');
        
        figure('Name',['inter-subject results:  Len = ' len_str '  ' params_t.accTypes{iAccum}]);
        subplot(1,2,1);bar(1:nof_files,crosss_valid_tp_median);hold on;
        errorbar(1:nof_files,crosss_valid_tp_median,crosss_valid_tp_std,'LineStyle','none');hold off;xlabel('subjects');title(['detection rate cross-validation median. average: ' num2str(mean(crosss_valid_tp_median))]);
        subplot(1,2,2);bar(1:nof_files,crosss_valid_tp_kappa_median);hold on;
        errorbar(1:nof_files,crosss_valid_tp_kappa_median,crosss_valid_tp_kappa_std,'LineStyle','none');xlabel('subjects');title(['detection rate cross-validation median Kappa. average: ' num2str(mean(crosss_valid_tp_kappa_median))]);
    end
end
