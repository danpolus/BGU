%
function params_t = global_params()

params_t.EOG_CHANNELS = [1 10 33 64];

params_t.max_tau_sec = 0.04; %maximal tau is 40msec
params_t.optimal_alpha = -1.5;
params_t.optimal_sigma = 1;
params_t.std_TH = 3.0; %try 3.5?
params_t.zscore_mode = 'file'; % epoch / file / scenario -> play with it
params_t.zero_subj_scenario_bad_channels_flg = false; %combine bad channels from same subject -> play with it

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
params_t.raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude' -> set to appear in folder name as well
params_t.similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'  -> set to appear in folder name as well
params_t.taus_to_use = 'optimal_multi_files'; %'all' 'optimal_per_file' 'optimal_multi_files'

params_t.contitionsExclude = {'cnd1wrd1', 'cnd1wrd2', 'cnd1wrd3'}; 
% params_t.contitionsExclude = [];
% params_t.contitionsUnite{1} = {'cnd1wrd1', 'cnd1wrd2', 'cnd1wrd3'};
% params_t.contitionsUnite{2} = {'cnd2wrd1', 'cnd2wrd2', 'cnd2wrd3'};
params_t.contitionsUnite = [];

params_t.test_set_percent = 0.15; % consider equalizing trainig sets: 10% for 4 conditions, 40% for 6 conditions, or equalizing testing sets: 15% for 4 conditions, 10% for 6 conditions

params_t.max_nof_clusters_per_condition = 3; %10^6;%take max contrast
params_t.nof_clusters_optimization = 'limit'; % scaled limit
params_t.maximal_distance_threshold = 1 - 0;% 1 - 1/2; %for cases when there is same similarity all over the matrix

params_t.accTypes = {'TotalAccum', 'ThreshAccum', 'EpochAccum', 'SampLimitAccum'};
%minimal_similarity_threshold : to decide if matching cluster found and avoid anti-similarity. has to be >0
switch params_t.similarity_method
    case 'jaccard'
        params_t.minimal_similarity_threshold = 1/2; %1/2; %0; %1/3;%3 times the radius for normalized
    case 'correlation'
        params_t.minimal_similarity_threshold = 0;
    case 'levenshtein'
        params_t.minimal_similarity_threshold = 1/2;
end
params_t.condition_descision_threshold = 0:0.5:5;
params_t.condition_counter_limit = [1 2 4 8 16 32 64 128 512 1024];
params_t.check_prediction_salience = false;
