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
