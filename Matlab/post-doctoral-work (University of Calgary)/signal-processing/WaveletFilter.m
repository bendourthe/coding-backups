function [ sig_filtered, Filter_test ] = WaveletFilter( sig, dt, cf, mode, type )
%Written by vvt Oct 20th, 2006.%
% Computes the wavelet coefficients of:
% wvlt = (f/cf)^mode * exp( (-f/cf+1) * mode )
%The rows of the matrix sig will be filtered.
%type indicates if high or low pass the cases are:
%1 low pass retain 100% of what is below cf.
%2 low pass remove everithing above cf.
%3 high pass retain 100 % of what is above cf.
%4 high pass remove everithing that is below cf.
%5 band pass what is defined by the wavelet.

%vvt this has to be adapted for long signals by chopping the signals.

%Determine length of signal
[rows input_length] = size(sig);

%adjust signal to be as close to 0 at the ends. readjust at the end.
bias = (sig(:,1) + sig(:,end))/2;
sig = sig - repmat(bias,1,input_length);

%Determin next power of 2 and expand the matrix, and df.
input_length_expanded = 2^ceil(log(input_length)/log(2));
T = input_length_expanded * dt;
df = 1/T;
sig(rows,input_length_expanded) = 0;

% index of Nyquist frequency
IND_NYQ = input_length_expanded/2;
frequency_index = 1:IND_NYQ+1;
%remember (frequency_index -1)*df = frequency if index starts at 1.
index_of_cf = floor(cf/df + 1);
f_cf = (frequency_index - 1) * (df/cf);
f_cf = f_cf(2:end);

% find wavelet coefficients as a function of the frequency
wvlt = exp( ( -f_cf + 1 +log(f_cf)) * mode );
wvlt = [0 wvlt]; %because log(0) is not possible.
wvlt_neg = 1 - wvlt;

switch type
    case 1
        filter = [ones(1,index_of_cf) wvlt(index_of_cf+1:IND_NYQ+1)];
    case 2
        filter = [wvlt_neg(1:index_of_cf) zeros(1,IND_NYQ+1 - index_of_cf)];
    case 3
        filter = [wvlt(1:index_of_cf) ones(1,IND_NYQ+1 - index_of_cf)];
    case 4
        filter = [zeros(1,index_of_cf)  wvlt_neg(index_of_cf+1:IND_NYQ+1)];
    otherwise
        filter = wvlt;
end

Filter_test.frequency=(frequency_index-1)*df; 
Filter_test.value = filter; 

filter = [filter  fliplr(filter(2:IND_NYQ))];
fsig = fft(sig, [], 2);
Filter_test.PSD=fsig(frequency_index).*conj(fsig(frequency_index));
Filter_test.PSD(1)=0;
Filter_test.PSD=Filter_test.PSD/max(Filter_test.PSD);
fsig_filtered = fsig .* repmat(filter,rows,1);
sig_filtered = ifft(fsig_filtered, [], 2);
sig_filtered = sig_filtered(:,1:input_length);
sig_filtered = sig_filtered + repmat(bias,1,input_length);
%end of function.