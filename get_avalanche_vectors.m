%
%Step 3 : prepare vectors and Ids
%
%  inputs:
% AvalancheFileDataSets - avalanches extracted from each EEGlab set
% saveFlg
%
%  outputs:
% MultiFileAchVecs - avalanche vectors from all sets
% usedTauInfo - used taus info (optimal or others)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MultiFileAchVecs, usedTauInfo] = get_avalanche_vectors(AvalancheFileDataSets, saveFlg)

%consider 60 or 64 electrodes vector
convert60to64channels_flg = 0;

params_t = global_params();
dataInfo = AvalancheFileDataSets(1).dataInfo;

if saveFlg
    output_fp = [dataInfo.FileInfo.base_fp '2 avalanches\'];
    mkdir(output_fp);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calc optimal tau

for iTau = 1:length(dataInfo.tau_vec)
    multi_files_epochs(iTau).av_size_vec = [];
    multi_files_epochs(iTau).sigma_vec = [];
end
for iAvalancheDataSets = 1:length(AvalancheFileDataSets)
    all_epochs = AvalancheFileDataSets(iAvalancheDataSets).all_epochs;
    [~,tau_optimal_inxs(iAvalancheDataSets)] = min(sqrt(([all_epochs.sigma] - params_t.optimal_sigma).^2 + (-[all_epochs.alpha] - params_t.optimal_alpha).^2));
    for iTau = 1:length(dataInfo.tau_vec)
        multi_files_epochs(iTau).av_size_vec = [multi_files_epochs(iTau).av_size_vec  all_epochs(iTau).av_size_vec];
        multi_files_epochs(iTau).sigma_vec = [multi_files_epochs(iTau).sigma_vec  all_epochs(iTau).sigma_vec];
    end
end
for iTau = 1:length(dataInfo.tau_vec)
    multi_files_epochs(iTau).alpha = estimateParamML(dataInfo.ES.min,dataInfo.ES.max,'zeta',-params_t.optimal_alpha,multi_files_epochs(iTau).av_size_vec);
    multi_files_epochs(iTau).sigma = mean(multi_files_epochs(iTau).sigma_vec);
end
[~,multi_files_tau_optimal_idx] = min(sqrt(([multi_files_epochs.sigma] - params_t.optimal_sigma).^2 + (-[multi_files_epochs.alpha] - params_t.optimal_alpha).^2));
if strcmp(params_t.taus_to_use, 'all')
    tau_idxs = 1:length(dataInfo.tau_vec);
else
    tau_idxs = unique([tau_optimal_inxs multi_files_tau_optimal_idx]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract avalanche vectors

MultiFileAchVecs = [];
for iAvalancheDataSets = 1:length(AvalancheFileDataSets)
    
    all_epochs = AvalancheFileDataSets(iAvalancheDataSets).all_epochs;
    setDataInfo = AvalancheFileDataSets(iAvalancheDataSets).dataInfo;

    file_id = ['scn' setDataInfo.FileInfo.scenario(1) 'sbj' setDataInfo.FileInfo.subj_id 'cnd' ...
        setDataInfo.FileInfo.condition(1) 'wrd' setDataInfo.FileInfo.word_num 'fln' num2str(iAvalancheDataSets,'%03d')];
    
    AvalancheVectors = [];
    for iTau = tau_idxs
        
        if strcmp(params_t.taus_to_use, 'optimal_multi_files') && iTau ~= multi_files_tau_optimal_idx && iTau ~= tau_optimal_inxs(iAvalancheDataSets)
            continue;
        end
        
        %prepare avalanche vectors
        AvalancheVectors(iTau).epochs_vecs = [];
        AvalancheVectors(iTau).Id = [];
        AvalancheVectors(iTau).tau = dataInfo.tau_vec(iTau);
        AvalancheVectors(iTau).is_optimal_tau = iTau == tau_optimal_inxs(iAvalancheDataSets);
        if convert60to64channels_flg
            AvalancheVectors(iTau).nof_channels = 64;
        else
            AvalancheVectors(iTau).nof_channels = size(all_epochs(1).av_raster_epochs,1);
        end
        AvalancheVectors(iTau).file_id = file_id;
        AvalancheVectors(iTau).dataInfo = setDataInfo;
        
        AllIds = [];
        for iEpoch = 1:size(all_epochs(iTau).av_raster_epochs,3)
            AvalancheVectors(iTau).epochs_vecs{iEpoch} = raster2vectors(all_epochs(iTau).av_raster_epochs(:,:,iEpoch), params_t, convert60to64channels_flg);
            if ~isempty(AvalancheVectors(iTau).epochs_vecs{iEpoch})
                avch_length_bins = [AvalancheVectors(iTau).epochs_vecs{iEpoch}.length_bins];
                for iAvalanche=1:length(avch_length_bins)
                    AvalancheVectors(iTau).epochs_vecs{iEpoch}(iAvalanche).id = [file_id 'epc' num2str(iEpoch,'%03d') 'ach' num2str(iAvalanche,'%04d')];
                    AllIds = [AllIds {AvalancheVectors(iTau).epochs_vecs{iEpoch}(iAvalanche).id}];
                end  
                if length(AvalancheVectors(iTau).Id) < max(avch_length_bins) %alocate indexes for new avalanche length
                    AvalancheVectors(iTau).Id{max(avch_length_bins)} = [];
                end
                for iLen = unique(avch_length_bins)
                    same_len_idx = find(avch_length_bins == iLen);
                    for iAvalanche = same_len_idx
                        AvalancheVectors(iTau).Id{iLen} = [AvalancheVectors(iTau).Id{iLen} {AvalancheVectors(iTau).epochs_vecs{iEpoch}(iAvalanche).id}];
                    end
                end              
            end
        end
        AvalancheVectors(iTau).Id{end+1} = AllIds; %last one is the all avalanches ids
              
    end %for iTau
    
    MultiFileAchVecs{iAvalancheDataSets} = AvalancheVectors;
    
end % for iAvalancheDataSets

usedTauInfo.tau_vec = dataInfo.tau_vec;
usedTauInfo.tau_idxs = tau_idxs;
usedTauInfo.multi_files_tau_optimal_idx = multi_files_tau_optimal_idx;

if saveFlg
    save([output_fp dataInfo.FileInfo.orig_fn '_avalanches.mat'],'MultiFileAchVecs','usedTauInfo');
end


% %debug plot
% figure;
% epoch_to_plot = 7;
% RasterPlot(all_epochs(multi_files_tau_optimal_idx).av_raster_epochs(:,:,epoch_to_plot),dataInfo.fs,dataInfo.tau_vec(multi_files_tau_optimal_idx));
% title(['Raster:  \tau = ' num2str(1000*dataInfo.tau_vec(multi_files_tau_optimal_idx)) ' ms    epoch #' num2str(epoch_to_plot)]);


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
function avalanche_vecs = raster2vectors(raster, params_t, convert60to64channels_flg)

avalanche_vecs = [];

if convert60to64channels_flg
    raster64 = zeros(64,size(raster,2));
    chan_inx = 1:64;
    chan_inx(params_t.EOG_CHANNELS) = [];
    raster64(chan_inx,:) = raster;
    raster = raster64;
end
if strcmp(params_t.raster_input_type,'binary')
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
