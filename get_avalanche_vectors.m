%
%Step 3
% prepare vectors and Ids
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_avalanche_vectors()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%consider 60 or 64 electrodes vector
convert60to64channels_flg = 0;

params_t = global_params();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select avalanche files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calc optimal tau

AvalancheFileData = load([fp files{1}],'all_epochs','param','file_info');
for iTau = 1:length(AvalancheFileData.param.tau_vec)
    multi_files_epochs(iTau).av_size_vec = [];
    multi_files_epochs(iTau).sigma_vec = [];
end
for iFile = 1:length(files)
    AvalancheFileData = load([fp files{iFile}]);
    [~,tau_optimal_inxs(iFile)] = min(sqrt(([AvalancheFileData.all_epochs.sigma] - AvalancheFileData.param.optimal_sigma).^2 + (-[AvalancheFileData.all_epochs.alpha] - AvalancheFileData.param.optimal_alpha).^2));
    for iTau = 1:length(AvalancheFileData.param.tau_vec)
        multi_files_epochs(iTau).av_size_vec = [multi_files_epochs(iTau).av_size_vec  AvalancheFileData.all_epochs(iTau).av_size_vec];
        multi_files_epochs(iTau).sigma_vec = [multi_files_epochs(iTau).sigma_vec  AvalancheFileData.all_epochs(iTau).sigma_vec];
    end
end
for iTau = 1:length(AvalancheFileData.param.tau_vec)
    multi_files_epochs(iTau).alpha = estimateParamML(AvalancheFileData.param.ES.min,AvalancheFileData.param.ES.max,'zeta',-AvalancheFileData.param.optimal_alpha,multi_files_epochs(iTau).av_size_vec);
    multi_files_epochs(iTau).sigma = mean(multi_files_epochs(iTau).sigma_vec);
end
[~,multi_files_tau_optimal_idx] = min(sqrt(([multi_files_epochs.sigma] - AvalancheFileData.param.optimal_sigma).^2 + (-[multi_files_epochs.alpha] - AvalancheFileData.param.optimal_alpha).^2));
if strcmp(params_t.taus_to_use, 'all')
    tau_idxs = 1:length(AvalancheFileData.param.tau_vec);
else
    tau_idxs = unique([tau_optimal_inxs multi_files_tau_optimal_idx]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract avalanche vectors
%compute similarity matrices for each file

multiple_files_fn = 'MultiFile';
MultiFileAchVecs = [];
for iFile = 1:length(files)
    
    AvalancheFileData = load([fp files{iFile}]);
    file_id = ['scn' AvalancheFileData.file_info.scenario(1) 'sbj' AvalancheFileData.file_info.subj_id 'cnd' ...
        AvalancheFileData.file_info.condition(1) 'wrd' AvalancheFileData.file_info.word_num 'fln' num2str(iFile,'%03d')];
    multiple_files_fn = [multiple_files_fn '_' file_id];
    
    AvalancheVectors = [];
    SimilarityMat = [];
    for iTau = tau_idxs
        
        if strcmp(params_t.taus_to_use, 'optimal_multi_files') && iTau ~= multi_files_tau_optimal_idx && iTau ~= tau_optimal_inxs(iFile)
            continue;
        end
        
%         tic
        display(['calc similarity: ' files{iFile} '  tau: ' num2str(AvalancheFileData.param.tau_vec(iTau))]);
        
        %prepare avalanche vectors
        AvalancheVectors(iTau).epochs_vecs = [];
        AvalancheVectors(iTau).Id = [];
        AvalancheVectors(iTau).IdLen = [];
        AvalancheVectors(iTau).tau = AvalancheFileData.param.tau_vec(iTau);
        AvalancheVectors(iTau).is_optimal_tau = iTau == tau_optimal_inxs(iFile);
        if convert60to64channels_flg
            AvalancheVectors(iTau).nof_channels = 64;
        else
            AvalancheVectors(iTau).nof_channels = size(AvalancheFileData.all_epochs(1).av_raster_epochs,1);
        end
        AvalancheVectors(iTau).file_id = file_id;
        for iEpoch = 1:size(AvalancheFileData.all_epochs(iTau).av_raster_epochs,3)
            AvalancheVectors(iTau).epochs_vecs{iEpoch} = raster2vectors(AvalancheFileData.all_epochs(iTau).av_raster_epochs(:,:,iEpoch), params_t.raster_input_type, convert60to64channels_flg);
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
        
%         %compute similarity matrices per file per tau (debug)
%         SimilarityMat(iTau).MatFull = [];
%         SimilarityMat(iTau).MatLen = [];
%         SimilarityMat(iTau).tau = AvalancheVectors(iTau).tau;
%         SimilarityMat(iTau).is_optimal_tau = AvalancheVectors(iTau).is_optimal_tau;
%         SimilarityMat(iTau).Id = AvalancheVectors(iTau).Id;
%         SimilarityMat(iTau).IdLen = AvalancheVectors(iTau).IdLen;
%         SimilarityMat = calc_similarity_mat(SimilarityMat, AvalancheVectors, iTau, params_t.similarity_method);
                
%         toc
    end %for iTau
    
    MultiFileAchVecs{iFile} = AvalancheVectors;
    
end % for iFile

tau.tau_vec = AvalancheFileData.param.tau_vec;
tau.tau_idxs = tau_idxs;
tau.multi_files_tau_optimal_idx = multi_files_tau_optimal_idx;
save([fp multiple_files_fn '.mat'],'MultiFileAchVecs','tau');


% %debug plot
% figure;
% epoch_to_plot = 7;
% RasterPlot(AvalancheFileData.all_epochs(multi_files_tau_optimal_idx).av_raster_epochs(:,:,epoch_to_plot),AvalancheFileData.param.Fs,AvalancheFileData.param.tau_vec(multi_files_tau_optimal_idx));
% title(['Raster:  \tau = ' num2str(1000*AvalancheFileData.param.tau_vec(multi_files_tau_optimal_idx)) ' ms    epoch #' num2str(epoch_to_plot)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% converts avalanche raster matrix to separate vectors
%
% INPUTS:
%   raster - avalanche raster matrix
%   raster_input_type - 'binary' 'events_cnt' 'amplitude'
%   convert60to64channels_flg
% OUTPUTS:
%   avalanche_vecs(i).vec
%   avalanche_vecs(i).first_bin_inx - avalanche onset bin index on the timeline
%   avalanche_vecs(i).length_bins - avalanche duration in bins
%
function avalanche_vecs = raster2vectors(raster, raster_input_type, convert60to64channels_flg)

avalanche_vecs = [];

if convert60to64channels_flg
    EOG_CHANNELS = [1 10 33 64];
    raster64 = zeros(64,size(raster,2));
    chan_inx = 1:64;
    chan_inx(EOG_CHANNELS) = [];
    raster64(chan_inx,:) = raster;
    raster = raster64;
end
if strcmp(raster_input_type,'binary')
    raster = double(raster ~= 0);
end

start_finish_mrk = diff([0 not(all(raster == 0,1)) 0]);
start_inx = find(start_finish_mrk == 1);
finish_inx = find(start_finish_mrk == -1) - 1;

for iAvalanche = 1:length(start_inx)
    avalanche_vecs(iAvalanche).vec = reshape(raster(:,start_inx(iAvalanche):finish_inx(iAvalanche)),1,[]);
    avalanche_vecs(iAvalanche).first_bin_inx = start_inx(iAvalanche);
    avalanche_vecs(iAvalanche).length_bins = finish_inx(iAvalanche) - start_inx(iAvalanche) + 1;
end
