% find_bad_epoch_channels - automaticaly finds noisy channels in epochs
%
% params:
%       channels_data
% options:
%       'minimal_nof_bad_pnts' - minimal number of threshold crossed points in epoch in channel in order to reject it
%       'z_score_thresh' - threshold for a noisy point in a channel
%
% return:
%       bad_epoch_chan - matrix of bad channels in epochs
%
function bad_epoch_chan = find_bad_epoch_channels(channels_data, varargin)

[n_chan,n_pts,n_epoch] = size(channels_data);
bad_epoch_chan = zeros(n_chan,n_epoch);

%defaults
minimal_nof_bad_pnts = 1;
z_score_thresh = 5;

%parse input
if nargin <= 1
    warning('find_bad_epoch_channels had no action');
    help find_bad_epoch_channels;
    return;
end

for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};

    switch Param
        case 'minimal_nof_bad_pnts'
            minimal_nof_bad_pnts = Value;
        case 'z_score_thresh'
            z_score_thresh = Value;
        otherwise
            warning(['Unknown input parameter ''' Param ''' ???'])
            help find_bad_epoch_channels;
            return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc z score
channels_data_concat = reshape(channels_data,n_chan, []);
% epMean_m = repmat(mean(channels_data_concat,2),1,n_pts);
epMean_m = repmat(median(channels_data_concat,2),1,n_pts);
epStd_m = repmat(std(channels_data_concat,[],2),1,n_pts);

for ep=1:n_epoch
    epZscore_m = abs((channels_data(:,:,ep) - epMean_m)./epStd_m);
    high_zScore_pnts_v = sum(epZscore_m > z_score_thresh,2);
    bad_epoch_chan(high_zScore_pnts_v >= minimal_nof_bad_pnts,ep) = 1;
end
