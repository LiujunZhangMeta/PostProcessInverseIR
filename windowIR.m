function [y,shifty] = windowIR(x,shiftx,fadeLen)
%WINDOWIR takes a time-domain signal, circularly shifts it, and applies a 
% hanning window to it. Useful for preparing an IR before zero-padding it.
% 
%USAGE
%     y = windowIR(x)
%     y = windowIR(x,shiftx)
%     y = windowIR(x,shiftx,fadeLen)
%
%INPUT PARAMETERS
%      x  : signal to be windowed
% fadeLen : length of the fade-in/fade-out [samples] - default value is 20%
%           of the signal length
%  shiftx : number of samples to circularly shift x before windowing - if 
%           empty, x will be centered around the highest peak
%
%OUTPUT PARAMETERS
%       y : windowed signal   
%  shifty : number of samples the signal was shifted
% 
%
%Author: Isaac Engel (August, 2018)

n = size(x,1);
nchannels = size(x,2);

if ~exist('fadeLen','var') || isempty(fadeLen)
    fadeLen = round(0.2*n);
end

if ~exist('shiftx','var') || isempty(shiftx)
    % First time-shift to center the IR around the highest peak
    [~,sigmax_ind]=max(mean(abs(x),2)); 
    shifty = -sigmax_ind+floor(n/2);
%     % Second time-shift to the right, to end the window at a zero-cross
%     y = circshift(x,shifty,1);
%     th = 0.25; % we will shift the window no more than 25% of its length
%     th_ind = length(y)-round(th*length(y)); 
%     [~,zerocross_ind]=min(mean(abs(y(th_ind:end,:)),2)); 
%     zerocross_ind = zerocross_ind + th_ind;
%     shifty = shifty + length(y) - zerocross_ind;
else
    shifty = shiftx;
end

y = circshift(x,shifty,1);  

if fadeLen > 0
    win_hann=hann(fadeLen*2-1);
    win_hann=[win_hann(1:fadeLen);
              ones(n-2*fadeLen,1);
              win_hann(fadeLen:end)];
    win_hann = repmat(win_hann,1,nchannels);

    y = y .* win_hann;
end

end

