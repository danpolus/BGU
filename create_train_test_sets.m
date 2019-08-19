%
%Step 5.2 : create_train_test_sets
%
%  inputs:
% MultiFileAchVecs - avalanche vectors from all sets
% usedTauInfo - used taus info (optimal or others)
%
%  outputs:
% TrainValidTest - train, test, cross-validation sets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function TrainValidTest = create_train_test_sets(MultiFileAchVecs, usedTauInfo)

params_t = global_params();

if strcmp(params_t.taus_to_use, 'optimal_multi_files')
    usedTauInfo.tau_idxs = usedTauInfo.multi_files_tau_optimal_idx;
end

fileInfo = MultiFileAchVecs{1}(usedTauInfo.multi_files_tau_optimal_idx).dataInfo.FileInfo;

nof_cross_valid = floor(1/params_t.test_set_percent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TrainValidTest = [];

for iTau = usedTauInfo.tau_idxs
    
    %get conditions and their epoch ids
    CondIds = {};
    ConditionsEpochIds = {};
    for iFile = 1:length(MultiFileAchVecs)
        if ~isempty(MultiFileAchVecs{iFile}(iTau).file_id)
            cond_idx = find(contains(CondIds, MultiFileAchVecs{iFile}(iTau).file_id(1:17)));
            if isempty(cond_idx)
                CondIds = [CondIds {MultiFileAchVecs{iFile}(iTau).file_id(1:17)}];
                cond_idx = length(CondIds);
                ConditionsEpochIds{cond_idx} = [];
            end
            for iEpoch=1:length(MultiFileAchVecs{iFile}(iTau).epochs_vecs)
                ConditionsEpochIds{cond_idx} = [ConditionsEpochIds{cond_idx} {[MultiFileAchVecs{iFile}(iTau).file_id(18:end) 'epc' num2str(iEpoch,'%03d')]}];
            end
        end
    end
    
    for iCond = 1:length(CondIds)
        nof_cond_epoch = length(ConditionsEpochIds{iCond});
        cond_epochs_inx = randperm(nof_cond_epoch);
        nof_test_epochs = nof_cond_epoch/nof_cross_valid;
        
        for iCrossValid = 1:nof_cross_valid
            test_epochs_inx = cond_epochs_inx(round((iCrossValid-1)*nof_test_epochs+1):round(iCrossValid*nof_test_epochs));
            if iCrossValid == 1     
                train_epochs_inx = setdiff(cond_epochs_inx,test_epochs_inx,'stable');
                TrainValidTest.TrainingSet(iTau).CondIds = CondIds;
                TrainValidTest.TrainingSet(iTau).EpochIds{iCond} = ConditionsEpochIds{iCond}(train_epochs_inx);
                TrainValidTest.TrainingSet(iTau).tau = usedTauInfo.tau_vec(iTau);
                TrainValidTest.TrainingSet(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
                TrainValidTest.TrainingSet(iTau).fileInfo.base_fp = fileInfo.base_fp;
                TrainValidTest.TrainingSet(iTau).fileInfo.orig_fn = fileInfo.orig_fn;   
                
                TrainValidTest.TestingSet(iTau).CondIds = CondIds;
                TrainValidTest.TestingSet(iTau).EpochIds{iCond}  = ConditionsEpochIds{iCond}(test_epochs_inx);
                TrainValidTest.TestingSet(iTau).tau = usedTauInfo.tau_vec(iTau);
                TrainValidTest.TestingSet(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
                TrainValidTest.TestingSet(iTau).fileInfo.base_fp = fileInfo.base_fp;
                TrainValidTest.TestingSet(iTau).fileInfo.orig_fn = fileInfo.orig_fn;  
            else
                valid_train_epochs_inx = setdiff(train_epochs_inx,test_epochs_inx,'stable');
                TrainValidTest.CrossValid(iCrossValid-1).TrainingSet(iTau).CondIds = CondIds;
                TrainValidTest.CrossValid(iCrossValid-1).TrainingSet(iTau).EpochIds{iCond} = ConditionsEpochIds{iCond}(valid_train_epochs_inx);
                TrainValidTest.CrossValid(iCrossValid-1).TrainingSet(iTau).tau = usedTauInfo.tau_vec(iTau);
                TrainValidTest.CrossValid(iCrossValid-1).TrainingSet(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
                TrainValidTest.CrossValid(iCrossValid-1).TrainingSet(iTau).fileInfo.base_fp = fileInfo.base_fp;
                TrainValidTest.CrossValid(iCrossValid-1).TrainingSet(iTau).fileInfo.orig_fn = fileInfo.orig_fn;  
                
                TrainValidTest.CrossValid(iCrossValid-1).TestingSet(iTau).CondIds = CondIds;
                TrainValidTest.CrossValid(iCrossValid-1).TestingSet(iTau).EpochIds{iCond} = ConditionsEpochIds{iCond}(test_epochs_inx);     
                TrainValidTest.CrossValid(iCrossValid-1).TestingSet(iTau).tau = usedTauInfo.tau_vec(iTau);
                TrainValidTest.CrossValid(iCrossValid-1).TestingSet(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
                TrainValidTest.CrossValid(iCrossValid-1).TestingSet(iTau).fileInfo.base_fp = fileInfo.base_fp;
                TrainValidTest.CrossValid(iCrossValid-1).TestingSet(iTau).fileInfo.orig_fn = fileInfo.orig_fn;                   
            end
        end
    end
    
end %for iTau
