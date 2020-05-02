function [out] = fft_real_matrix(x,inv)
%x is an input matrix containing the time series
%fft_real_vvt returns the fourier transform of the signal
%The length of the signal should be divisible by 2 or a power of 2.
%out is forward FFT if inv = 0; the imag part is the projection onto the
%-sin(2pi * f) axes. thus the -imag is the actual fourier coeficient.
%out is the inverse fft if inv = 1;
%out is the powerspectrum if inv = 2;
%T = N*dt;
%df = 1/T;
%frequency = (k-1) * df; with k = (1:Nyquist_frequency_ind);
%Nyquist_frequency = Nyquist_frequency_ind * df
% power is: pfx = fx .* conj(fx); %power is independent of dt. Power has unities of x squared.
%sum(pfx) = sum(x.*x) or sum(pfx) * dt = sum(x.*x) *dt is energy. dt = T/N
%pfx/df is power spectral density. This makes the amplitude of the spectrum
%independent of T or in turn of 0 padding.
%the first value of fft_real_vvt(x) / sqrt(N) is mean(x).


if inv == 0 || inv == 2 %do forward fourier transform
        si = size(x);
        col_vector=true;
    if si(1) < si(2) %then we have a row vector and chang it to a column
        x = x';
        col_vector=false;
    end
    N = size(x,1);
    Nyquist_frequency_ind = N/2 + 1;
    normalization_operator = zeros(Nyquist_frequency_ind,1) + sqrt(2/N);
    normalization_operator(1) = 1/sqrt(N);
    normalization_operator(end) = 1/sqrt(N);
    fx = fft(x);
    dum=repmat(normalization_operator,1,size(fx,2));
    out = fx(1:Nyquist_frequency_ind,:) .* dum; % this is real_fft of vincent
    if ~col_vector; out=out'; end
    if inv == 2
        out = out .* conj(out);%this is = pfx, power spectrum; x .* conj(x) = (abs(x))^2
    end
else %do inverse fourier transform
            si = size(x);
        col_vector=true;
    if si(1) < si(2) %then we have a row vector and chang it to a column
        x = x';
        col_vector=false;
    end

    Nyquist_frequency_ind = size(x,1);
    N = (Nyquist_frequency_ind -1) * 2;
    normalization_operator = zeros(Nyquist_frequency_ind,1) + sqrt(2/N);
    normalization_operator(1) = 1/sqrt(N);normalization_operator(end) = 1/sqrt(N);
    dum=repmat(normalization_operator,1,size(x,2));
    x = x ./ dum;
    x2 = flipud(x(2:Nyquist_frequency_ind -1,:));
    x = [x;conj(x2)];
    %a multiplication with a complec number e.g. when computing a time shift
    %causes the degenerated fourier coeficients to become complex. 
    %They have to be set to real.
    x(1,:) = real(x(1,:)); x(Nyquist_frequency_ind,:) = real(x(Nyquist_frequency_ind,:));
    out = ifft(x);
    if ~col_vector; out=out'; end

end

