%
%Step 4 : compare_avalanches - split to training and testing sets; generate similarity matrix
%
%  inputs:
% MultiFileAchVecs - avalanche vectors from all sets
% usedTauInfo - used taus info (optimal or others)
% saveFlg
% plotFlg
%
%  outputs:
% SimilarityMat - similarity between avalanches distances matrix
% TestingSet - avalanches (epochs) for testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SimilarityMat, TestingSet] = compare_avalanches(MultiFileAchVecs, usedTauInfo, saveFlg, plotFlg)

params_t = global_params();

if strcmp(params_t.taus_to_use, 'optimal_multi_files')
    usedTauInfo.tau_idxs = usedTauInfo.multi_files_tau_optimal_idx;
end

fileInfo = MultiFileAchVecs{1}(usedTauInfo.multi_files_tau_optimal_idx).dataInfo.FileInfo;

if saveFlg
    output_fp = [fileInfo.base_fp '3 similarity\'];
    mkdir(output_fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimilarityMat = [];
TestingSet = [];

for iTau = usedTauInfo.tau_idxs
    tic
    
    %prepare testing and training sets
    TestingSet(iTau).CondIds = {};
    TestingSet(iTau).AchIds = [];
    TestingSet(iTau).EpochIds = [];
    TestingSet(iTau).tau = usedTauInfo.tau_vec(iTau);
    TestingSet(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
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
        ConditionsAchIds{cond_idx} = [ConditionsAchIds{cond_idx} MultiFileAchVecs{iFile}(iTau).Id{end}];
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
    SimilarityMat(iTau).Mat = [];
    SimilarityMat(iTau).Id = [];
    SimilarityMat(iTau).tau = usedTauInfo.tau_vec(iTau);
    SimilarityMat(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
    AllIds = [];
    for iFile = 1:length(MultiFileAchVecs)
        cond_idx = find(contains(TestingSet(iTau).CondIds, MultiFileAchVecs{iFile}(iTau).file_id(1:17)));
        AllIds = [AllIds  ConditionsAchIds{cond_idx}];
        
        max_avch_length_bins = length(MultiFileAchVecs{iFile}(iTau).Id)-1; %last one is all avalanche matrix
        for iLen = 1:max_avch_length_bins
            if length(SimilarityMat(iTau).Id ) < max_avch_length_bins %alocate indexes for new avalanche length
                SimilarityMat(iTau).Id{max_avch_length_bins} = [];
            end
            if ~isempty(MultiFileAchVecs{iFile}(iTau).Id{iLen})
                train_samples_inx_Len = contains(ConditionsAchIds{cond_idx}, MultiFileAchVecs{iFile}(iTau).Id{iLen});
                SimilarityMat(iTau).Id{iLen} = [SimilarityMat(iTau).Id{iLen}  ConditionsAchIds{cond_idx}(train_samples_inx_Len)];
            end
        end
    end
    SimilarityMat(iTau).Id{end+1} = AllIds;
    SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, params_t.similarity_method);

    if plotFlg %plot
        figure; 
        subplot(1,2,1); imagesc(SimilarityMat(iTau).Mat{end}); title(['all avalanches similarity,  \tau =' num2str(SimilarityMat(iTau).tau)]);
        subplot(1,2,2); imagesc(SimilarityMat(iTau).Mat{3}); title(['3 bins length avalanches similarity,  \tau =' num2str(SimilarityMat(iTau).tau)]);
    end

    display(['calc similarity matrix: ' fileInfo.orig_fn  ' tau: ' num2str(usedTauInfo.tau_vec(iTau))]);
    toc
end %for iTau

if saveFlg
    save([output_fp fileInfo.orig_fn '_similarity.mat'],'MultiFileAchVecs','SimilarityMat','TestingSet');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, similarity_method)

%             fln = str2num(SimilarityMat(iTau).Id{iLen}{m}(21:23));
%             epc = str2num(SimilarityMat(iTau).Id{iLen}{m}(27:29));
%             ach = str2num(SimilarityMat(iTau).Id{iLen}{m}(33:36)); 

    nextLen = 1;
    
    if strcmp(similarity_method, 'jaccard') %faster computation
        for iLen = 1:length(SimilarityMat(iTau).Id)-1
            V = [];
            if length(SimilarityMat(iTau).Id{iLen}) == 1
                SimilarityMat(iTau).Mat{iLen} = 0;
            else
                for m=1:length(SimilarityMat(iTau).Id{iLen})
                    V = [V; MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{iLen}{m}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).Id{iLen}{m}(27:29))}(str2num(SimilarityMat(iTau).Id{iLen}{m}(33:36))).vec];
                end
                SimilarityMat(iTau).Mat{iLen} = squareform(1 - pdist(V,'jaccard'), 'tomatrix');
            end
        end
        nextLen = length(SimilarityMat(iTau).Id);
    end

    for iLen = nextLen:length(SimilarityMat(iTau).Id)
        SimilarityMat(iTau).Mat{iLen} = zeros(length(SimilarityMat(iTau).Id{iLen}));
        for m=1:length(SimilarityMat(iTau).Id{iLen})
            for n=1:m-1
                v1 = MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{iLen}{m}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).Id{iLen}{m}(27:29))}(str2num(SimilarityMat(iTau).Id{iLen}{m}(33:36))).vec;
                v2 = MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{iLen}{n}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).Id{iLen}{n}(27:29))}(str2num(SimilarityMat(iTau).Id{iLen}{n}(33:36))).vec;
                if iLen < length(SimilarityMat(iTau).Id)
                    SimilarityMat(iTau).Mat{iLen}(m,n) = calc_vectors_similarity(v1, v2, similarity_method);
                else %full matrix
                    fln = str2num(SimilarityMat(iTau).Id{iLen}{m}(21:23));
                    SimilarityMat(iTau).Mat{iLen}(m,n) = calc_sliding_similarity(v1, v2, MultiFileAchVecs{fln}(iTau).nof_channels, similarity_method);
                end
            end
        end
        SimilarityMat(iTau).Mat{iLen} = SimilarityMat(iTau).Mat{iLen} + SimilarityMat(iTau).Mat{iLen}';
    end
    