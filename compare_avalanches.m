%
%Step 4
%split to training and testing sets
%generate similarity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compare_avalanches()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fn, fp] = uigetfile([fp '*.mat'], 'Select avalanche vectors file');
load([fp fn],'MultiFileAchVecs','tau');

if strcmp(params_t.taus_to_use, 'optimal_multi_files')
    tau.tau_idxs = tau.multi_files_tau_optimal_idx;
end

SimilarityMat = [];
TestingSet = [];
for iTau = tau.tau_idxs
    tic
    
    %prepare testing and training sets
    TestingSet(iTau).CondIds = {};
    TestingSet(iTau).AchIds = [];
    TestingSet(iTau).EpochIds = [];
    TestingSet(iTau).tau = tau.tau_vec(iTau);
    TestingSet(iTau).is_optimal_tau = iTau == tau.multi_files_tau_optimal_idx;
    ConditionsAchIds = {};
    ConditionsEpochIds = {};
    for iFile = 1:length(MultiFileAchVecs)
        cond_idx = find(contains(TestingSet(iTau).CondIds, MultiFileAchVecs{iFile}(iTau).file_id(1:17)));
        if isempty(cond_idx)
            TestingSet(iTau).CondIds = [TestingSet(iTau).CondIds {MultiFileAchVecs{iFile}(iTau).file_id(1:17)}];
            cond_idx = length(TestingSet(iTau).CondIds);
            ConditionsAchIds{cond_idx} = [];
            ConditionsEpochIds{cond_idx} = [];
        end
        ConditionsAchIds{cond_idx} = [ConditionsAchIds{cond_idx} MultiFileAchVecs{iFile}(iTau).Id];
        for iEpoch=1:length(MultiFileAchVecs{iFile}(iTau).epochs_vecs)
            ConditionsEpochIds{cond_idx} = [ConditionsEpochIds{cond_idx} {[MultiFileAchVecs{iFile}(iTau).file_id(18:end) 'epc' num2str(iEpoch,'%03d')]}];
        end
    end
    
    for cond_idx=1:length(TestingSet(iTau).CondIds)
        test_epochs_inx = randperm(length(ConditionsEpochIds{cond_idx}), round(length(ConditionsEpochIds{cond_idx})*params_t.test_set_percent));
        TestingSet(iTau).EpochIds{cond_idx} = ConditionsEpochIds{cond_idx}(test_epochs_inx);
        test_ach_inx = contains(ConditionsAchIds{cond_idx},TestingSet(iTau).EpochIds{cond_idx});
        TestingSet(iTau).AchIds{cond_idx} = ConditionsAchIds{cond_idx}(test_ach_inx);
        ConditionsAchIds{cond_idx}(test_ach_inx) = [];
    end

    
    %compute similarity matrix
    display(['calc similarity all files,  tau: ' num2str(tau.tau_vec(iTau))]);
    
    SimilarityMat(iTau).MatFull = [];
    SimilarityMat(iTau).MatLen = [];
    SimilarityMat(iTau).tau = tau.tau_vec(iTau);
    SimilarityMat(iTau).is_optimal_tau = iTau == tau.multi_files_tau_optimal_idx;
    SimilarityMat(iTau).Id = [];
    SimilarityMat(iTau).IdLen = [];
    for iFile = 1:length(MultiFileAchVecs)
        cond_idx = find(contains(TestingSet(iTau).CondIds, MultiFileAchVecs{iFile}(iTau).file_id(1:17)));
        SimilarityMat(iTau).Id = [SimilarityMat(iTau).Id  ConditionsAchIds{cond_idx}];
        
        max_avch_length_bins = length(MultiFileAchVecs{iFile}(iTau).IdLen);
        for iLen = 1:max_avch_length_bins
            if length(SimilarityMat(iTau).IdLen ) < max_avch_length_bins %alocate indexes for new avalanche length
                SimilarityMat(iTau).IdLen{max_avch_length_bins} = [];
            end
            if ~isempty(MultiFileAchVecs{iFile}(iTau).IdLen{iLen})
                train_samples_inx_Len = contains(ConditionsAchIds{cond_idx},MultiFileAchVecs{iFile}(iTau).IdLen{iLen});
                SimilarityMat(iTau).IdLen{iLen} = [SimilarityMat(iTau).IdLen{iLen}  ConditionsAchIds{cond_idx}(train_samples_inx_Len)];
            end
        end
    end
    SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, params_t.similarity_method);

    %plot
    figure; 
    subplot(1,2,1); imagesc(SimilarityMat(iTau).MatFull); title(['all avalanches similarity,  \tau =' num2str(SimilarityMat(iTau).tau)]);
    subplot(1,2,2); imagesc(SimilarityMat(iTau).MatLen{3}); title(['3 bins length avalanches similarity,  \tau =' num2str(SimilarityMat(iTau).tau)]);


toc
end %for iTau

save([fp fn '_similarity.mat'],'MultiFileAchVecs','SimilarityMat','TestingSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, similarity_method)

%             fln = str2num(SimilarityMat(iTau).Id{m}(21:23));
%             epc = str2num(SimilarityMat(iTau).Id{m}(27:29));
%             ach = str2num(SimilarityMat(iTau).Id{m}(33:36));

    SimilarityMat(iTau).MatFull = zeros(length(SimilarityMat(iTau).Id));
    for m=1:length(SimilarityMat(iTau).Id)
        for n=1:m-1
            v1 = MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{m}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).Id{m}(27:29))}(str2num(SimilarityMat(iTau).Id{m}(33:36))).vec;
            v2 = MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{n}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).Id{n}(27:29))}(str2num(SimilarityMat(iTau).Id{n}(33:36))).vec;
            SimilarityMat(iTau).MatFull(m,n) = calc_sliding_similarity(v1, v2, MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{1}(21:23))}(iTau).nof_channels, similarity_method);
        end
    end
    SimilarityMat(iTau).MatFull = SimilarityMat(iTau).MatFull + SimilarityMat(iTau).MatFull';
    
    for iLen = 1:length(SimilarityMat(iTau).IdLen)
        SimilarityMat(iTau).MatLen{iLen} = zeros(length(SimilarityMat(iTau).IdLen{iLen}));
        for m=1:length(SimilarityMat(iTau).IdLen{iLen})
            for n=1:m-1
                v1 = MultiFileAchVecs{str2num(SimilarityMat(iTau).IdLen{iLen}{m}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).IdLen{iLen}{m}(27:29))}(str2num(SimilarityMat(iTau).IdLen{iLen}{m}(33:36))).vec;
                v2 = MultiFileAchVecs{str2num(SimilarityMat(iTau).IdLen{iLen}{n}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).IdLen{iLen}{n}(27:29))}(str2num(SimilarityMat(iTau).IdLen{iLen}{n}(33:36))).vec;
                SimilarityMat(iTau).MatLen{iLen}(m,n) = calc_vectors_similarity(v1, v2, similarity_method);
            end
        end
        SimilarityMat(iTau).MatLen{iLen} = SimilarityMat(iTau).MatLen{iLen} + SimilarityMat(iTau).MatLen{iLen}';
    end
    