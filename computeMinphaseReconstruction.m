%COMPUTEMINPHASERECONSTRUCTION Compute minphase reconstruction
%% INPUT
% input_mag : magnitude of left ear freq. domain signal 
%% OUTPUT
% output    : minphase reconstructed time-domain signal 
%% NOTES
% minimum-phase reconstruction:   Given the magnitude of freq. domain signal |X(f)| 
%                              1) minimum phase  phi = -Hilbert_transform(log(|X(f)|) (in MATLAB hilbert transform is imag(hilbert())
%                              2) X_R(f) = |X(f)| exp(i phi)
%                              3) x_R(t) = inverse_fourier_transform(X_R(f)) (in MATLAB use real(ifft())
function [ output ] = computeMinphaseReconstruction( input_mag, num_timesamples )
   input_mag(input_mag==0)=eps; % prevent zero values in magnitude
   isRow = size(input_mag,1) == 1;
    if( isRow )
        input_mag = transpose(input_mag);
    end
    if( mod(num_timesamples,2) == 0 ) % even
        complete_mag_response = [ input_mag(1:end);  conj(input_mag(end-1:-1:2)) ];
    else                        % odd 
        complete_mag_response = [ input_mag(1:end);  conj(input_mag(end:-1:2)) ];
    end
   phase =  -imag( hilbert( real(log(complete_mag_response)) ) );
   freq_response = complete_mag_response .* exp(1i*phase);
   freq_response(1) = real(freq_response(1));
   %output = convertFreqToTimeResponse( freq_response(1:round(end/2)+1), num_timesamples );
   output = ifft(freq_response,num_timesamples,1, 'symmetric');
   if( isRow )
        output = transpose(output);
    end
end

