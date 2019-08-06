%
%Step 4 : compare_avalanches - generate similarity matrix
%
%  inputs:
% MultiFileAchVecs - avalanche vectors from all sets
% usedTauInfo - used taus info (optimal or others)
% saveFlg
% plotFlg
%
%  outputs:
% SimilarityMat - similarity between avalanches distances matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SimilarityMat = compare_avalanches(MultiFileAchVecs, usedTauInfo, saveFlg, plotFlg)

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

for iTau = usedTauInfo.tau_idxs
    percent_waitbar = 0;
    f_waitbar = waitbar(percent_waitbar, ['calc similarity matrix \tau=' num2str(usedTauInfo.tau_vec(iTau)) '   ' num2str(100*percent_waitbar) '%'], 'Name',fileInfo.orig_fn );
     
    %compute similarity matrix   
    SimilarityMat(iTau).Mat = [];
    SimilarityMat(iTau).Id = [];
    SimilarityMat(iTau).tau = usedTauInfo.tau_vec(iTau);
    SimilarityMat(iTau).is_optimal_tau = iTau == usedTauInfo.multi_files_tau_optimal_idx;
    AllIds = [];
    for iFile = 1:length(MultiFileAchVecs)
        AllIds = [AllIds  MultiFileAchVecs{iFile}(iTau).Id{end}];
        max_avch_length_bins = length(MultiFileAchVecs{iFile}(iTau).Id)-1; %last one is all avalanche matrix
        if length(SimilarityMat(iTau).Id ) < max_avch_length_bins %alocate indexes for new avalanche length
            SimilarityMat(iTau).Id{max_avch_length_bins} = [];
        end
        for iLen = 1:max_avch_length_bins
            if ~isempty(MultiFileAchVecs{iFile}(iTau).Id{iLen})
                SimilarityMat(iTau).Id{iLen} = [SimilarityMat(iTau).Id{iLen}  MultiFileAchVecs{iFile}(iTau).Id{iLen}];
            end
        end
    end
    SimilarityMat(iTau).Id{end+1} = AllIds;
    SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, params_t.similarity_method, f_waitbar);

    if plotFlg %plot
        figure; 
        subplot(1,2,1); imagesc(SimilarityMat(iTau).Mat{end}); title(['all avalanches similarity,  \tau =' num2str(SimilarityMat(iTau).tau)]);
        subplot(1,2,2); imagesc(SimilarityMat(iTau).Mat{3}); title(['3 bins length avalanches similarity,  \tau =' num2str(SimilarityMat(iTau).tau)]);
    end

    close(f_waitbar);
end %for iTau

if saveFlg
    save([output_fp fileInfo.orig_fn '_similarity.mat'],'MultiFileAchVecs','SimilarityMat','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, similarity_method, f_waitbar)

%             fln = str2num(SimilarityMat(iTau).Id{iLen}{m}(21:23));
%             epc = str2num(SimilarityMat(iTau).Id{iLen}{m}(27:29));
%             ach = str2num(SimilarityMat(iTau).Id{iLen}{m}(33:36)); 

    nextLen = 1;
    
    if strcmp(similarity_method, 'jaccard') %faster computation
        for iLen = 1:length(SimilarityMat(iTau).Id)-1
            percent_waitbar = iLen/length(SimilarityMat(iTau).Id);
            waitbar(percent_waitbar,f_waitbar,['calc similarity matrix \tau=' num2str(SimilarityMat(iTau).tau) '   ' num2str(100*percent_waitbar) '%']);
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
            percent_waitbar = (iLen - 1 + m/length(SimilarityMat(iTau).Id{iLen}))/length(SimilarityMat(iTau).Id);
            waitbar(percent_waitbar,f_waitbar,['calc similarity matrix \tau=' num2str(SimilarityMat(iTau).tau) '   ' num2str(100*percent_waitbar) '%']);
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
    
%     plot_hist = plot_hist + v1;
%     figure;stem([plot_hist(max_inx,:)/max(plot_hist(max_inx,:))/0.5; v]');legend('cluster histogram','vector');title(['cluster: ' num2str(cluster_num(iLen)) '   similarity: ' num2str(cluster_sim(iLen))]);  
    