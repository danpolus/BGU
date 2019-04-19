
function params_t = global_params()

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
params_t.raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude' -> set to appear in folder name as well
params_t.similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'  -> set to appear in folder name as well
params_t.taus_to_use = 'optimal_multi_files'; %'all' 'optimal_per_file' 'optimal_multi_files'


params_t.test_set_percent = 0.1;

params_t.reasonable_nof_clusters = 100;
params_t.max_nof_clusters = Inf; % 15 Inf
params_t.minimal_contrast = 0.5;


params_t.compare_modes = {'Full', 'Len', 'ConcatLen'};%{'Len'};

%minimal_similarity_threshold : to decide if matching cluster found and avoid anti-similarity
switch params_t.similarity_method
    case 'jaccard'
        params_t.minimal_similarity_threshold = 0.5;
    case 'correlation'
        params_t.minimal_similarity_threshold = 0;
    case 'levenshtein'
        params_t.minimal_similarity_threshold = 0.5;
end
nof_thresh_to_test = 10;
% params_t.condition_descision_threshold = linspace(0,6,nof_thresh_to_test);  %average-based salience (to chose this method, search 'average-based salience' through project files)
params_t.condition_descision_threshold = linspace(0.5,1,nof_thresh_to_test);  %contrast-based salience (to chose this method, search 'contrast-based salience' through project files)
params_t.condition_counter_limit = linspace(1,30,nof_thresh_to_test);
