
function params_t = global_params()

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
params_t.raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude' -> set to appear in folder name as well

params_t.similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'  -> set to appear in folder name as well

params_t.test_set_percent = 0.1;

params_t.minimal_similarity_threshold = 0.5;
params_t.condition_descision_threshold = 0.5;
params_t.accum_sample_limit = 10; %max nof avalanches in accumulator for SampLimitAccum mode
