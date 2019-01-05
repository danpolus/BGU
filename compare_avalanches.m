%
%Step 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compare_avalanches()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude' -> set to appear in folder name as well

similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'  -> set to appear in folder name as well

%consider 60 or 64 electrodes vector
convert60to64channels_flg = 0;

taus_to_use = 'optimal_multi_files'; %'all' 'optimal_per_file' 'optimal_multi_files'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select avalanche files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

%calc optimal tau
AvalancheData = load([fp files{1}]);
for iTau = 1:length(AvalancheData.param.tau_vec)
    multi_files_epochs(iTau).av_size_vec = [];
    multi_files_epochs(iTau).sigma_vec = [];
end
for iFile = 1:length(files)
    AvalancheData = load([fp files{iFile}]);
    [~,tau_optimal_inxs(iFile)] = min(sqrt(([AvalancheData.all_epochs.sigma] - AvalancheData.param.optimal_sigma).^2 + (-[AvalancheData.all_epochs.alpha] - AvalancheData.param.optimal_alpha).^2));
    for iTau = 1:length(AvalancheData.param.tau_vec)
        multi_files_epochs(iTau).av_size_vec = [multi_files_epochs(iTau).av_size_vec  AvalancheData.all_epochs(iTau).av_size_vec];
        multi_files_epochs(iTau).sigma_vec = [multi_files_epochs(iTau).sigma_vec  AvalancheData.all_epochs(iTau).sigma_vec];
    end
end
for iTau = 1:length(AvalancheData.param.tau_vec)
    multi_files_epochs(iTau).alpha = estimateParamML(AvalancheData.param.ES.min,AvalancheData.param.ES.max,'zeta',-AvalancheData.param.optimal_alpha,multi_files_epochs(iTau).av_size_vec);
    multi_files_epochs(iTau).sigma = mean(multi_files_epochs(iTau).sigma_vec);
end
[~,multi_files_tau_optimal_idx] = min(sqrt(([multi_files_epochs.sigma] - AvalancheData.param.optimal_sigma).^2 + (-[multi_files_epochs.alpha] - AvalancheData.param.optimal_alpha).^2));
if strcmp(taus_to_use, 'all')
    tau_idxs = 1:length(AvalancheData.param.tau_vec);
else
    tau_idxs = unique([tau_optimal_inxs multi_files_tau_optimal_idx]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute similarity matrices for each file

multiple_files_fn = 'MultiFile';
MultiFileAchVecs = [];
for iFile = 1:length(files)
    
    AvalancheData = load([fp files{iFile}]);
    file_id = ['scn' AvalancheData.file_info.scenario(1) 'sbj' AvalancheData.file_info.subj_id 'cnd' ...
        AvalancheData.file_info.condition(1) 'wrd' AvalancheData.file_info.word_num 'fln' num2str(iFile,'%03d')];
    multiple_files_fn = [multiple_files_fn '_' file_id];
    
    AvalancheVectors = [];
    SimilarityMat = [];
    for iTau = tau_idxs
        
        if strcmp(taus_to_use, 'optimal_multi_files') && iTau ~= multi_files_tau_optimal_idx && iTau ~= tau_optimal_inxs(iFile)
            continue;
        end
        
        tic
        display(['calc similarity: ' files{iFile} '  tau: ' num2str(AvalancheData.param.tau_vec(iTau))]);
        
        %prepare avalanche vectors
        AvalancheVectors(iTau).epochs_vecs = [];
        AvalancheVectors(iTau).Id = [];
        AvalancheVectors(iTau).IdLen = [];
        AvalancheVectors(iTau).tau = AvalancheData.param.tau_vec(iTau);
        AvalancheVectors(iTau).is_optimal_tau = iTau == tau_optimal_inxs(iFile);
        if convert60to64channels_flg
            AvalancheVectors(iTau).nof_channels = 64;
        else
            AvalancheVectors(iTau).nof_channels = size(AvalancheData.all_epochs(1).av_raster_epochs,1);
        end
        for iEpoch = 1:size(AvalancheData.all_epochs(iTau).av_raster_epochs,3)
            AvalancheVectors(iTau).epochs_vecs{iEpoch} = raster2vectors(AvalancheData.all_epochs(iTau).av_raster_epochs(:,:,iEpoch), raster_input_type, convert60to64channels_flg);
            if ~isempty(AvalancheVectors(iTau).epochs_vecs{iEpoch})
                avch_length_bins = [AvalancheVectors(iTau).epochs_vecs{iEpoch}.length_bins];
                for iAvalanche=1:length(avch_length_bins)
                    AvalancheVectors(iTau).Id = [AvalancheVectors(iTau).Id {[file_id 'epc' num2str(iEpoch,'%03d') 'ach' num2str(iAvalanche,'%04d')]}];
                end
                if length(AvalancheVectors(iTau).IdLen) < max(avch_length_bins) %alocate indexes for new avalanche length
                    AvalancheVectors(iTau).IdLen{max(avch_length_bins)} = [];
                end
                for iLen = unique(avch_length_bins)
                    same_len_idx = find(avch_length_bins == iLen);
                    for iAvalanche = same_len_idx
                        AvalancheVectors(iTau).IdLen{iLen} = [AvalancheVectors(iTau).IdLen{iLen} {[file_id 'epc' num2str(iEpoch,'%03d') 'ach' num2str(iAvalanche,'%04d')]}];
                    end
                end
            end
        end
        
        MultiFileAchVecs{iFile} = AvalancheVectors;
        
        %compute similarity matrices
        SimilarityMat(iTau).MatFull = [];
        SimilarityMat(iTau).MatLen = [];
        SimilarityMat(iTau).tau = AvalancheVectors(iTau).tau;
        SimilarityMat(iTau).is_optimal_tau = AvalancheVectors(iTau).is_optimal_tau;
        SimilarityMat(iTau).Id = AvalancheVectors(iTau).Id;
        SimilarityMat(iTau).IdLen = AvalancheVectors(iTau).IdLen;
        SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, similarity_method);
                
        toc
    end %for iTau
    
    save([fp files{iFile}(1:end-4) '_similarity.mat'],'AvalancheVectors','SimilarityMat');
end % for iFile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute similarity matrices for all the files
SimilarityMat = [];
if strcmp(taus_to_use, 'optimal_multi_files')
    tau_idxs = multi_files_tau_optimal_idx;
end
for iTau = tau_idxs
    tic
    display(['calc similarity all files,  tau: ' num2str(AvalancheData.param.tau_vec(iTau))]);
    
    SimilarityMat(iTau).MatFull = [];
    SimilarityMat(iTau).MatLen = [];
    SimilarityMat(iTau).tau = AvalancheData.param.tau_vec(iTau);
    SimilarityMat(iTau).is_optimal_tau = iTau == multi_files_tau_optimal_idx;
    SimilarityMat(iTau).Id = [];
    SimilarityMat(iTau).IdLen = [];
    for iFile = 1:length(files)
        SimilarityMat(iTau).Id = [SimilarityMat(iTau).Id  MultiFileAchVecs{iFile}(iTau).Id];
        max_avch_length_bins = length(MultiFileAchVecs{iFile}(iTau).IdLen);
        for iLen = 1:max_avch_length_bins
            if length(SimilarityMat(iTau).IdLen ) < max_avch_length_bins %alocate indexes for new avalanche length
                SimilarityMat(iTau).IdLen{max_avch_length_bins} = [];
            end
            SimilarityMat(iTau).IdLen{iLen} = [SimilarityMat(iTau).IdLen{iLen}  MultiFileAchVecs{iFile}(iTau).IdLen{iLen}];
        end
    end
    SimilarityMat = calc_similarity_mat(SimilarityMat, MultiFileAchVecs, iTau, similarity_method);
    
toc
end %for iTau

save([fp multiple_files_fn '_similarity.mat'],'MultiFileAchVecs','SimilarityMat');

% figure;
% epoch_to_plot = 7;
% RasterPlot(AvalancheData.all_epochs(AvalancheData.tau_optimal_inx).av_raster_epochs(:,:,epoch_to_plot),AvalancheData.param.Fs,AvalancheData.param.tau_vec(AvalancheData.tau_optimal_inx));
% title(['Raster:  \tau = ' num2str(1000*AvalancheData.param.tau_vec(AvalancheData.tau_optimal_inx)) ' ms    epoch #' num2str(epoch_to_plot)]);



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
    