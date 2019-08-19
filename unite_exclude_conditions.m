%
%%Step 5.1 : unite and/or exclude conditions from data (if requested)
%
%  inputs:
% MultiFileAchVecs - avalanche vectors from all sets
% SimilarityMat - similarity between avalanches distances matrix
% usedTauInfo - used taus info (optimal or others)
%
%  outputs:
% MultiFileAchVecs - avalanche vectors from all sets
% SimilarityMat - similarity between avalanches distances matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [MultiFileAchVecs, SimilarityMat] = unite_exclude_conditions(MultiFileAchVecs, SimilarityMat, usedTauInfo, saveFlg)

params_t = global_params();

if strcmp(params_t.taus_to_use, 'optimal_multi_files')
    usedTauInfo.tau_idxs = usedTauInfo.multi_files_tau_optimal_idx;
end

cond_save_str = [];
for iUni = 1:length(params_t.contitionsUnite)
    cond_save_str = [cond_save_str ' un_'];
    for iCond = 1:length(params_t.contitionsUnite{iUni})
        cond_save_str = [cond_save_str 'c' params_t.contitionsUnite{iUni}{iCond}(4) params_t.contitionsUnite{iUni}{iCond}(8)];
    end
end
if ~isempty(params_t.contitionsExclude)
    cond_save_str = [cond_save_str ' ex_'];
end
for iCond = 1:length(params_t.contitionsExclude)
    cond_save_str = [cond_save_str 'c' params_t.contitionsExclude{iCond}(4) params_t.contitionsExclude{iCond}(8)];
end
if isempty(cond_save_str)
    return;
end
    
if saveFlg
    fileInfo = MultiFileAchVecs{1}(usedTauInfo.multi_files_tau_optimal_idx).dataInfo.FileInfo;
    output_fp = [fileInfo.base_fp '3 similarity' cond_save_str '\' ];
    mkdir(output_fp);
end

%%%%%%%%%%%%%%%%%

for iTau = usedTauInfo.tau_idxs

    %exclude in MultiFileAchVecs
    if ~isempty(params_t.contitionsExclude)
        for iFile = 1:length(MultiFileAchVecs)
            if any(contains(params_t.contitionsExclude, MultiFileAchVecs{iFile}(iTau).file_id(10:17)))
                MultiFileAchVecs{iFile}(iTau).file_id = [];
            end
        end

        %exclude in SimilarityMat
        for iLen = 1:length(SimilarityMat(iTau).Id)
            if ~isempty(SimilarityMat(iTau).Id{iLen})
                ach_idx = contains(SimilarityMat(iTau).Id{iLen}, params_t.contitionsExclude);
                SimilarityMat(iTau).Mat{iLen} = SimilarityMat(iTau).Mat{iLen}(~ach_idx,~ach_idx);
                SimilarityMat(iTau).Id{iLen}(ach_idx) = [];
            end
        end
    end
    
    %%%%%%%
    
    for iUni = 1:length(params_t.contitionsUnite)
        
        %unite in MultiFileAchVecs
        for iFile = 1:length(MultiFileAchVecs)
            if ~isempty(MultiFileAchVecs{iFile}(iTau).file_id) && any(contains(params_t.contitionsUnite{iUni}(2:end), MultiFileAchVecs{iFile}(iTau).file_id(10:17)))
                
                for iEpoch = 1:length(MultiFileAchVecs{iFile}(iTau).epochs_vecs)
                    for iVec = 1:length(MultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch})
                        MultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch}(iVec).id(10:17) = params_t.contitionsUnite{iUni}{1};
                    end
                end
                
                for iLen = 1:length(MultiFileAchVecs{iFile}(iTau).Id)
                    for iAch = 1:length(MultiFileAchVecs{iFile}(iTau).Id{iLen})
                        MultiFileAchVecs{iFile}(iTau).Id{iLen}{iAch}(10:17) = params_t.contitionsUnite{iUni}{1};
                    end
                end                
                
                MultiFileAchVecs{iFile}(iTau).file_id(10:17) = params_t.contitionsUnite{iUni}{1};
                
                if strcmp(params_t.contitionsUnite{iUni}{1}(4),'1')
                    MultiFileAchVecs{iFile}(iTau).dataInfo.FileInfo.condition = '1rest';
                else
                    MultiFileAchVecs{iFile}(iTau).dataInfo.FileInfo.condition = '2imagine';
                end
                MultiFileAchVecs{iFile}(iTau).dataInfo.FileInfo.word_num = params_t.contitionsUnite{iUni}{1}(8);

            end
        end

        %unite in SimilarityMat
        for iLen = 1:length(SimilarityMat(iTau).Id)
            if ~isempty(SimilarityMat(iTau).Id{iLen})
                achIdx = find(contains(SimilarityMat(iTau).Id{iLen}, params_t.contitionsUnite{iUni}(2:end)));
                for iAch = achIdx
                    SimilarityMat(iTau).Id{iLen}{iAch}(10:17) = params_t.contitionsUnite{iUni}{1};
                end
            end
        end
            
    end
    
end

%%%%%%%%%%%%%%%%%
if saveFlg
    save([output_fp fileInfo.orig_fn '_similarity.mat'],'MultiFileAchVecs','SimilarityMat','usedTauInfo','-v7.3');
end
