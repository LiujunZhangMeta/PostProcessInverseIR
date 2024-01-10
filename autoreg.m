function [y, y_mp] = autoreg(x,config)
%AUTOREG calculates the regularized inverse of an impulse response,
%   according to [1] and [3]
% 
%USAGE
%          y = autoreg(x,config)
%  [y, y_mp] = autoreg(x,config)
%
%INPUT PARAMETERS
%        x : impulse response to be inverted (1 or more channels)
%   config : struct object containing the settings
%
% Config fields:
%       f1 : low cut-off frequency; if <= 0, no low cutoff will be used
%            (low-pass regularization)
%       f2 : high cut-off frequency; if >= srate/2, no high cutoff will be 
%            used (high-pass regularization)
%    srate : sample rate [Hz] (default = 48000)
%     taps : number of taps of the output inverse filter (default = 2048)
%    fType : type of filter to be used - 1: FIR (default), 2: Butterworth
%   order1 : order of the high-pass Butterworth filter used in the 
%            regularization (default = 5)
%   order2 : order of the low-pass Butterworth filter used in the 
%            regularization (default = 5)
%   config.option : 0: do not preserve low-frequencies
%            1: preserve low frequencies 
%            2: preserve high frequencies       
%            3: preserve both low AND high frequencies (default)
%  max_amp : maximum amplification allowed by the regularization [dB]
%            (default = 20)
%        k : smoothing window length used to obtain sigma
%      phi : minimum-phase factor (0<phi<1); if close to 1, it will reduce
%            the pre-ringing of the filter, but it will require a longer
%            filter [3] (default = 0)
% targetIR : IR whose magnitude we will try to approach (default = 1, which
%            produces flat EQ)
%
%OUTPUT PARAMETERS
%       y : regularized inverse impulse response   
%    y_mp : minimum phase equivalent of y
% 
%
%Author: Isaac Engel (August, 2018)
%
%REFERENCES
%    [1] Bolaños, Javier Gómez, Aki Mäkivirta, and Ville Pulkki. "Automatic
%        regularization parameter for headphone transfer function 
%        inversion." Journal of the Audio Engineering Society 64.10 (2016):
%        752-761.
%    [2] Schärer, Zora, and Alexander Lindau. "Evaluation of equalization
%        methods for binaural signals." Audio Engineering Society 
%        Convention 126. Audio Engineering Society, 2009.
%    [3] Bouchard, M., Norcross, S. G., & Soulodre, G. A. (2006, October). 
%        Inverse filtering design using a minimal-phase target function 
%        from regularization. In Audio Engineering Society Convention 121. 
%        Audio Engineering Society.

%% Check input arguments

if ~exist('config','var')
    warning('Usage: autoreg(x,config). Using default parameters...')
    config = struct();
end

if ~isfield(config,'srate')
    config.srate = 48000;
    warning('Sample rate not specified - defaulting to 48000...');
end

if ~isfield(config,'f1')
    config.f1 = 200; % recommended by [2]
end

if ~isfield(config,'f2')
    config.f2 = 16000; % recommended by [2]
end

if ~isfield(config,'taps')
    config.taps = 2048; % recommended by [2]
end

if mod(config.taps,2)
    config.taps = config.taps + 1; % force taps to be even
end

if ~isfield(config,'order1')
    config.order1 = 5;
end

if ~isfield(config,'order2')
    config.order2 = 5;
end

if ~isfield(config,'option')
    config.option = 3;
end

if ~isfield(config,'max_amp')
    config.max_amp = 20;
end

if ~isfield(config,'k')
    config.k = 0.5;
end

if ~isfield(config,'phi')
    config.phi = 0;
end

if ~isfield(config,'targetIR')
    config.targetIR = AKdirac(size(x,1));
end


%% Main variables
taps = config.taps;
srate = config.srate;
out_taps = taps; % number of taps of output filter
taps = max([2^11,2*taps,length(x)]); % number of taps to be used for the calculations
nchannels = size(x,2); % number of channels
win_len = length(x); % impulse response length in samples  
fvec = (srate*(0:taps/2)./taps)'; % vector of frequencies

%% Pre-processing

% Window the signal
if taps >= win_len
    x = [x;zeros(taps-win_len,nchannels)]; % zero pad
else
    fadesamples = round(taps*0.2);
    [~,~,ind] = trimSignal(x,taps); % keep the maximum amount of energy
    x = circshift(x,-(ind-fadesamples-1)); % shift to the right
    x = x(1:taps,:); % trim
    x = windowIR(x,0,fadesamples); % hannin
    % g window
end

% Same with the target IR
targetIR = config.targetIR;
if taps >= length(targetIR)
    targetIR = [targetIR;zeros(taps-size(targetIR,1),size(targetIR,2))];
else
    fadesamples = round(taps*0.2);
    [~,~,ind] = trimSignal(targetIR,taps);
    targetIR = circshift(targetIR,-(ind-fadesamples-1));
    targetIR = targetIR(1:taps,:);
    targetIR = windowIR(targetIR,0,fadesamples);
end

%% FFT

X = fft(x,taps); % impulse response in the frequency domain
H = abs(X(1:end/2+1,:)); % magnitude up to Nyquist freq.

%% Alpha parameter (equalization bandwidth)

% Low cutoff/high-pass filter
if config.f1 <= 0 || config.option == 1 || config.option == 3
    HP = ones(taps/2+1,1);
else
    [z,p,kk] = butter(config.order1,config.f1/(srate/2),'high');
    sos = zp2sos(z,p,kk); 
    [HP,~] = freqz(sos,(taps/2)+1);
    HP(end)=real(HP(end)); % ensure that it is symmetric    
end

% High cutoff/low-pass filter
if config.f2 >= srate/2 || config.option == 2 || config.option == 3
    LP = ones(taps/2+1,1);
else
    [z,p,kk] = butter(config.order2,config.f2/(srate/2),'low');
    sos = zp2sos(z,p,kk); 
    [LP,~] = freqz(sos,(taps/2)+1);
    LP(end)=real(LP(end)); % ensure that it is symmetric    
end

W = HP .* LP; % multiply both filters
W = repmat(W,1,nchannels); % repeat for every channel   

alpha = 1./abs(W.^2) - 1;
min_alpha = 1/(2*db2mag(config.max_amp))^2; % max_amp = 1/(2*sqrt(beta))
alpha = alpha + min_alpha; % limit the amplification

%% Sigma parameter (to avoid inversion of deep notches, following [1])

if config.k == 0
    Hs = H;
else
    % Obtain a smooth version of the magnitude curve
%     Hs = zeros(size(H));
%     for i=1:nchannels
%         Hs(:,i) = erbsmooth(H(:,i),srate,config.k); 
%     end
    Hs = AKfractOctSmooth(H,'amp',srate,2/config.k,false); % the factor 2 and the inversion was to try and make it behave similarly to the previous erbsmooth (which was removed because it was slow)
end

sigma = zeros(size(Hs));
ind = abs(Hs)>=abs(H);
sigma(ind) = abs(Hs(ind)) - abs(H(ind));

%% Beta parameter (total regularization)

beta = alpha + sigma.^2;

beta = [beta; flipud(beta(2:end-1,:))]; % reconstruct full spectrum

%% Apply regularization

% Factor in the target IR before the regularization
targetFR = fft(targetIR,taps);
X = X ./ targetFR; % if targetIR is a delta (= 1), this has no effect

% Calculate regularization factor
A = (abs(X).^2) ./ (abs(X).^2 + beta);

% Calculate minimum phase equivalent, following [3]
A(A==0)=eps;
phase = -imag(hilbert(log(A)));
A_minphase = A.*exp(1i*phase*config.phi);

Y = (1./X) .* A_minphase; % regularized inverse
y=ifft(Y,'symmetric'); % back to time domain

%% Experimental: magnitude preservation at the extremes

if config.option > 0

    Hy = abs(Y(1:end/2+1,:)); % magnitude up to Nyquist freq.

    % We want to cross over between two magnitude curves
    h1 = ones(size(Hy)); % 0dB
    h2 = Hy; % reg. inverse of X

    % Warp to logarithmic scale
    fvec_log = logspace(log10(1),log10(fvec(end)),length(h1))';
    h1_log = interp1(fvec,h1,fvec_log,'spline');
    h2_log = interp1(fvec,h2,fvec_log,'spline');

    % Window
    w = ones(size(Hy,1),1);

    % Below the low cutoff
    if config.option == 1 || config.option == 3
        fstart = config.f1 * 2^(-1); % start fade one octave before
        fend = config.f1 * 2^(0); % end fade at indicated frequency

        indstart = find(fvec_log<fstart,1,'last');
        indend = find(fvec_log<fend,1,'last');
        len = indend-indstart;

        if len > 0
            t = (0:len)'.*(pi/2)/len;
            slope = sin(t).^2; % ascending slope
            w(1:indstart,:) = 0;
            w(indstart:indend,:) = slope;
        end
    end

    % Above the high cutoff
    if config.option == 2 || config.option == 3
        fstart = config.f2 * 2^(0); % start fade at indicated frequency
        fend = config.f2 * 2^(1); % end fade one octave later

        indstart = find(fvec_log<fstart,1,'last');
        indend = find(fvec_log<fend,1,'last');
        len = indend-indstart;

        if len > 0
            t = (0:len)'.*(pi/2)/len;
            slope = 1 - sin(t).^2; % descending slope
            w(indstart:indend,:) = slope;
            w(indend:end,:) = 0;
        end
    end

    Hy_log = real( h1_log.^(1-w) .* h2_log.^w ); % real because sometimes we get imag values here even if we shouldn't
    Hy = interp1(fvec_log,Hy_log,fvec,'spline'); % back to linear
    Hy = [Hy; flipud(Hy(2:end-1,:))]; % reconstruct full spectrum
    Y = Hy .* exp(1i*angle(Y)); % keep the same phase as the original

    y=ifft(Y,'symmetric'); % back to time domain
    
end

%% Center, trim and window the signal

y = windowIR(y,[],0); % center IR (do not apply window yet)

% Trim (keep the central samples)
if rem(out_taps,2) == 0 % even taps
    ind = taps/2-out_taps/2+1 : taps/2+out_taps/2;
else % odd taps
    ind = taps/2-(out_taps-1)/2 : taps/2+(out_taps-1)/2;
end
y = y(ind,:);

y = y.*hann(out_taps); % apply a full hanning window
% y = windowIR(y,0,32); % apply a few samples of windowing

%% Get minimum phase equivalent

% Get minimum phase reconstruction from the inverse filter magnitude
halfmag = abs(Y(1:end/2+1,:));
y_mp = zeros(taps,nchannels);
for i=1:nchannels
    y_mp(:,i) = computeMinphaseReconstruction(halfmag(:,i),taps);
end

% Trim (keep the first samples)
y_mp = y_mp(1:out_taps,:);

win_hann = hann(2*out_taps-1);
y_mp = y_mp.*win_hann(out_taps:end); % apply hanning window

end
