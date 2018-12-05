%
%Step 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude'

similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'

%consider 60 or 64 electrodes vector
convert60to64channels_flg = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select avalanche files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

for iFile = 1:length(files)
    AvalancheData = load([fp files{iFile}]);
    
    for iTau = 1:length(AvalancheData.param.tau_vec)
        tic
        
        %prepare avalanche vectors
        AvalancheVectors(iTau).epochs_vecs = [];
        AvalancheVectors(iTau).Id = [];
        AvalancheVectors(iTau).IdLen = [];
        AvalancheVectors(iTau).tau = AvalancheData.param.tau_vec(iTau);
        if convert60to64channels_flg
            AvalancheVectors(iTau).nof_channels = 64;
        else
            AvalancheVectors(iTau).nof_channels = size(AvalancheData.all_epochs(1).av_raster_epochs,1);
        end
        for iEpoch = 1:size(AvalancheData.all_epochs(iTau).av_raster_epochs,3)       
            AvalancheVectors(iTau).epochs_vecs{iEpoch} = raster2vectors(AvalancheData.all_epochs(iTau).av_raster_epochs(:,:,iEpoch), raster_input_type, convert60to64channels_flg);
            if ~isempty(AvalancheVectors(iTau).epochs_vecs{iEpoch})
                avch_length_bins = [AvalancheVectors(iTau).epochs_vecs{iEpoch}.length_bins];
                AvalancheVectors(iTau).Id = [AvalancheVectors(iTau).Id [iEpoch*ones(1,length(avch_length_bins)); 1:length(avch_length_bins)]];
                if length(AvalancheVectors(iTau).IdLen) < max(avch_length_bins) %alocate indexes for new avalanche length
                    AvalancheVectors(iTau).IdLen{max(avch_length_bins)} = [];
                end
                for iLen = unique(avch_length_bins)
                    same_len_idx = find(avch_length_bins == iLen);
                    AvalancheVectors(iTau).IdLen{iLen} = [AvalancheVectors(iTau).IdLen{iLen} [iEpoch*ones(1,length(same_len_idx)); same_len_idx]];
                end
            end           
        end
              
        %compute similarity matrices
        SimilarityMat(iTau).MatFull = [];
        SimilarityMat(iTau).MatLen = [];
        SimilarityMat(iTau).tau = AvalancheData.param.tau_vec(iTau);
        for m=1:size(AvalancheVectors(iTau).Id,2)
            for n=1:size(AvalancheVectors(iTau).Id,2)
                v1 = AvalancheVectors(iTau).epochs_vecs{AvalancheVectors(iTau).Id(1,m)}(AvalancheVectors(iTau).Id(2,m)).vec;
                v2 = AvalancheVectors(iTau).epochs_vecs{AvalancheVectors(iTau).Id(1,n)}(AvalancheVectors(iTau).Id(2,n)).vec;
                SimilarityMat(iTau).MatFull(m,n) = calc_sliding_similarity(v1, v2, AvalancheVectors(iTau).nof_channels, similarity_method);                
            end
        end
        for iLen = 1:length(AvalancheVectors(iTau).IdLen)
            ids = AvalancheVectors(iTau).IdLen{iLen};
            for m=1:size(AvalancheVectors(iTau).IdLen{iLen},2)
                for n=1:size(AvalancheVectors(iTau).IdLen{iLen},2)
                    v1 = AvalancheVectors(iTau).epochs_vecs{AvalancheVectors(iTau).IdLen{iLen}(1,m)}(AvalancheVectors(iTau).IdLen{iLen}(2,m)).vec;
                    v2 = AvalancheVectors(iTau).epochs_vecs{AvalancheVectors(iTau).IdLen{iLen}(1,n)}(AvalancheVectors(iTau).IdLen{iLen}(2,n)).vec;
                    SimilarityMat(iTau).MatLen{iLen}(m,n) = calc_vectors_similarity(v1, v2, similarity_method);
                end
            end
        end
          
        toc
    end %for iTau
    
    save([fp files{iFile}(1:end-4) '_similarity.mat'],'AvalancheVectors','SimilarityMat');
end % for iFile

% figure;
% epoch_to_plot = 7;
% RasterPlot(AvalancheData.all_epochs(AvalancheData.tau_optimal_inx).av_raster_epochs(:,:,epoch_to_plot),AvalancheData.param.Fs,AvalancheData.param.tau_vec(AvalancheData.tau_optimal_inx));
% title(['Raster:  \tau = ' num2str(1000*AvalancheData.param.tau_vec(AvalancheData.tau_optimal_inx)) ' ms    epoch #' num2str(epoch_to_plot)]);
