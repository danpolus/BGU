% find_bad_cahnnels - automaticaly finds outrageous channels and returns their indeces
%
% params:
%       channels_data
%       fs
% options:
%       'window_step_length_sec' - step size for analysis
%       'bad_time_percent' - what percent on the channel should be noisy (compared to other channels) so the channel would be marked as bad
%       'max_amp_score_thresh' - threshold (z-score based on median) for a noisy section (compared to other channels in this section)
%       'min_amp_score_thresh' - threshold (z-score based on median) for a silent section (compared to other channels in this section)
%
% return:
%       bad_channels - vector of bad_channels
%
function bad_channels = find_bad_cahnnels(channels_data, fs, varargin)

bad_channels = [];

%defaults
window_step_length = floor(60 * fs); % 1 minute window step
bad_time_percent = 0.5; % 50%
max_amp_score_thresh = 5;
min_amp_score_thresh = 0.2;


%parse input
if nargin <= 2
    warning('find_bad_cahnnels had no action');
    help find_bad_cahnnels;
    return;
end
if mod(nargin,2) == 1
    warning('Odd number of input arguments')
    help find_bad_cahnnels;
    return;
end

for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    %         if ~isstr(Param)
    %             error('Flag arguments must be strings')
    %         end
    %         Param = lower(Param);
    switch Param
        case 'window_step_length_sec'
            window_step_length = floor(Value * fs);
        case 'bad_time_percent'
            bad_time_percent = Value;
        case 'max_amp_score_thresh'
            max_amp_score_thresh = Value;   
        case 'min_amp_score_thresh'
            min_amp_score_thresh = Value;            
        otherwise
            warning(['Unknown input parameter ''' Param ''' ???'])
            help find_bad_cahnnels;
            return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channelStandardSmplitudes_m = [];
next_section_idx = 1;
for i = 1:window_step_length:size(channels_data,2)-2*window_step_length+1
    channelStandardSmplitudes_m = [channelStandardSmplitudes_m std(channels_data(:,i:i+window_step_length-1),[],2)];
    next_section_idx = next_section_idx + window_step_length;
end
channelStandardSmplitudes_m = [channelStandardSmplitudes_m std(channels_data(:,next_section_idx:end),[],2)];
medAmplitudes_m = repmat(median(channelStandardSmplitudes_m,1), size(channels_data,1), 1);

over_max_amp_time_percent = mean(channelStandardSmplitudes_m > max_amp_score_thresh*medAmplitudes_m, 2);
bad_channels = [bad_channels  find(over_max_amp_time_percent > bad_time_percent)'];

under_min_amp_time_percent = mean(channelStandardSmplitudes_m < min_amp_score_thresh*medAmplitudes_m, 2);
bad_channels = [bad_channels  find(under_min_amp_time_percent > bad_time_percent)']; 

bad_channels = unique(bad_channels);

