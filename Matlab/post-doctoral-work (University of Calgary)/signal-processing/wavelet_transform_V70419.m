
function [power, cwt, pa] = ...
    wavelet_transform_V70419(sig,sampling_rate,new_sampling_rate,scale,control,nr_wavelets)
%Programmed by V. von Tscharner for Peter Brueggemann
%Does the wavelet transform of the electromyograms and creates the wavelets.
%sampling rate should be in the range 1000 to 10000 Hz.
%cwt is the complex wavelet transform of the signal;
%power is the Gauss filtered power ...
%   of the complex wavelet transformed signal cwt.
%power has the format time as rows and wavelet_power as columns.
%The sampling_rate automatically defines the number of wavelets; max 13.
%pa are the parameters for the computation. They are returned to the user.
%control=1 forces the wavelet creation. Has to be used if scale is changed.
%control=2 no gauss filter is applied to the power output.
%control=3, no effect, default.
%control=4 do not use power normalization and gauss filter.
%Use negative control values for Morlet wavelets. Not recommended.

%This function refers to itself. Change if name is changed.

persistent wavelets wavelets1 wavelets2 fwavelets wpa npa
%wpa are wavelet parameter transfered to the user
%npa are local parmeters used when the function is caled again.


if not(exist('control'))
    %display('no control submitted');
    control = 3;
else
    if sign(control) < 0
        control = -control;
        npa.flag.morlet = true;
    else
        npa.flag.morlet = false;
    end
end

if size(wavelets,1) == 0 % this is the first run of the function
control = 1; end

modulus = sampling_rate/new_sampling_rate;
if not(isempty(wavelets)) && not(wpa.modulus == modulus)
    wpa.modulus = modulus;    
    control = 1;
    disp('changed modulus in exist')
end

npa.flag.generate_wavelets = false;
switch control
    case 1
        if isempty(nr_wavelets)
            wpa.nr_of_wavelets = 18;
        else
            wpa.nr_of_wavelets = nr_wavelets;
        end
        npa.control = control;
        npa.flag.use_ck = true;
        npa.flag.use_gauss = true;
        npa.flag.generate_wavelets = true;
        npa.alpha = 150; % parameter of Morlet wavelet
        
    case 2 % exclude gauss filter
        npa.flag.use_gauss = false;
        
    case 3 % reinstall gauss filter
        npa.flag.use_gauss = true;
    
    case 4
        % if control = 4 do not use power normalization
        npa.control = control;
        npa.flag.use_ck = false;
        npa.flag.use_gauss = false;
        npa.flag.generate_wavelets = true;
end


if sampling_rate < 999
    disp('sampling_rate is too low');
    return
end

%change orientation of signal if col is used instead or row.
orientation=true; %signal is entered as row vector
[r,c] = size(sig);
if r > c
    sig = sig';
    orientation=false;
end

%Test input parameters
if isempty(wavelets)
    wpa.resolution.tr_ch = [];
    npa.tr_ch_present = 0;
    npa.crop = [];   
end


[tmax jmax] = size(wavelets);

% create wavelets for the first run or if control = 1
if npa.flag.generate_wavelets || tmax==0
    npa.split = 4; %split the wavelets in two groups one from wavelet 1:4
    % the other from 4:jmax. This is used in the wavelet transform and
    % for gauss filtering. 
    % On higher wavelets the duration can be reduced by croping the signal.
    % Computational load and time is saved.
    
    wpa.scale = scale;
    wpa.original_sampling_rate = sampling_rate;
    wpa.new_sampling_rate = new_sampling_rate;
    wpa.modulus = wpa.original_sampling_rate/wpa.new_sampling_rate;
    
    %increase length if wavelet 1 is not returning to 0.
    wpa.length_wavelets_log2 = ceil(log2(2^9/2000*sampling_rate*wpa.scale/0.3));
    wpa.width_factor_gauss = 3/8;
    
    %compute the theoretical maximum number of wavelets.
    j = 2;
    cfs1 = 1/wpa.scale*(1.45+j-2).^1.959;
    while cfs1 < sampling_rate/5
        j = j+1;
        cfs1 = 1/wpa.scale*(1.45+j-2).^1.959;
    end
    if j < wpa.nr_of_wavelets;
        wpa.nr_of_wavelets = j;
        disp(' maximal number of wavelets is '); disp(j);
    end
    %starting wavelet construction
    [wavelets, fwavelets, cfs] = Tscharner_wavelet(wpa,npa);
    [tmax jmax] = size(wavelets);
    wpa.cfs = cfs;
    wpa.resolution.tr_ch = [];
    npa.tr_ch_present = 0;
    wpa.resolution = resolution(wpa,wavelets,fwavelets,npa);
    wres = wpa.resolution.tr_ch * (wpa.original_sampling_rate/wpa.new_sampling_rate);
    npa.crop = floor((length(wavelets) - 4*wres(npa.split+1))/2);
    %Separate the wavelets for the split operation
    
    wavelets1 = wavelets(:,1:npa.split);
    wavelets2 = wavelets(1+npa.crop:end-npa.crop,npa.split+1:end);
    npa.tr_ch_present = 1;
    npa.flag.use_gauss = true;
end

% test if signal is long enough to do a wavelet transform
t1max = size(sig,2);
if t1max<2^wpa.length_wavelets_log2
%     disp('signal is too short');
end


le = tmax; % length of wavelets
%take every modulus'th point
t1bmax = ceil(t1max/wpa.modulus);
srb = wpa.new_sampling_rate; dtb = 1/srb;
t1b = [1:t1bmax]; %time points in the resampled array

%perform wavelet transform on sig
front_mirror = sig(le/2:-1:1);
end_mirror = sig(end:-1:end-le/2+1);
extended_sig = [front_mirror sig end_mirror];
%t2 is range of signal in extended_sig
t2 = [le/2+1:le/2+t1max];

%create a matrix where each col represent a signal to be transformed.
    %for very long signals an iterative process should be used.
sig_mat = zeros(t1bmax,le);
le_1 = le-1;
for n = t1b
    from = round((n-1)*wpa.modulus +1);
    sig_mat(n,:)=extended_sig(from:from+le_1);
end

%to use without split activate next line.
%npa.tr_ch_present = 0;

%perform wavelet transform. 
% Using wavelets1 and wavelets2 is faster.
% cwt contains real and complex values of wavelet transformed signal
if npa.tr_ch_present
    cwt = [sig_mat*wavelets1,sig_mat(:,1+npa.crop:end-npa.crop)*wavelets2];
else
    cwt = sig_mat*wavelets;
end
clear sig_mat

%prepare output of wavelet transform function
pa = wpa; % parameters
pa.wavelets = wavelets; %make wavelets available to user.
pa.fwavelets=fwavelets;
power = conj(cwt).*cwt; %here the power is computed

pa.resolution.gauss = []; %Clear gauss filter values, no output
%If control = 2 then do not use the gauss filter;
if npa.flag.use_gauss;
    leg = size(wpa.resolution.gauss,1);
    leg_1 = leg-1;
    split = npa.split;
    crop1 = wpa.resolution.g(1);
    gauss1 = wpa.resolution.gauss(1+crop1:end-crop1,1:split)';
    sum_op1 = ones(size(gauss1,2),1);
    crop2 = wpa.resolution.g(split+1);
    gauss2 = wpa.resolution.gauss(1+crop2:end-crop2,split+1:end)';
    sum_op2 = ones(size(gauss2,2),1);
    nmax = size(power,1);
    power_filtered = power;
    %mirror image pattern at both ends
    fr_mi = power(leg/2:-1:1,:);
    end_mi = power(end:-1:end-leg/2+1,:);
    power = [fr_mi;power;end_mi]';
    power1 = power(1:split,:);power2 = power(split+1:end,:);
    for n =1:nmax
        powerl1 = power1(:,n:n+leg_1);
        dum1 = (powerl1 .* gauss1) * sum_op1;
        powerl2 = power2(:,n+crop2:n+leg_1 - crop2);
        dum2 = (powerl2 .* gauss2) * sum_op2;
        power_filtered(n,:) = [dum1;dum2];
    end    
    power = power_filtered;
end

end % end of function wavelet transform        
        
%____________functions____________________________________________________

function [wave fwave cfs] = Tscharner_wavelet(wpa,npa)
%WAVELET creats an array of length 2 to power2 with intervals df around the
%center frequency cf with a width defined by the mode.
%wave are the wavelets in time domaine fwave in frequency domaine.
jmax_final = wpa.nr_of_wavelets;
scale = wpa.scale;
sampling_rate = wpa.original_sampling_rate;
power2 = wpa.length_wavelets_log2;
split = npa.split;

le = 2^power2;
dt = 1/sampling_rate;
tmax = le*dt;
df = sampling_rate/le;
fmax = le/2*df;
x = [1:le/2+1];
frequency = df*(x-1);

%compute center frequencies and mode for Tscharner_wavelets.
%we compute one wavelet more at the low and high end of the array and
%then delete them at the end. This way the two end wavelets are not
%distorted.
jmax = jmax_final+2;
j = [1:jmax];
cfs = 1/scale*(1.45+j-2).^1.959; %array of cf in Hz.
cfs_ch = round(cfs/df);


%make matrix with f/cfs(j) in each col indexed j and compute the wavelets
inv_cfs(j) = 1./cfs(j);
diag_inv_cfs = diag(inv_cfs);
f_mat = repmat(frequency',1,jmax);
f_mat(1,:) = eps;
frcf = f_mat*diag_inv_cfs;

% Morlet_wavelets
if npa.flag.morlet
    alpha = npa.alpha;
    fac_1 = (2*pi^2/alpha) * cfs;
    fac_2 = ((frcf - 1).^2);
    wavelet_exponent = - fac_2 * diag(fac_1);
else
% Cauchy_Wavelets
    modes = scale*cfs;
    diag_modes = diag(modes);
    log_frcf = log(frcf);
    wavelet_exponent = ((-frcf+1)+log_frcf)*diag_modes;
end

%limit wavelet_exponent otherwied wavelets may become NaN.
if any(any(wavelet_exponent < -300))
    dum = find(wavelet_exponent < -300);
    wavelet_exponent(dum) = -300; 
end
f_wavelets = exp(wavelet_exponent);

%normalize to plateau value.
sum_operator = ones(jmax,1);
plateau = f_wavelets * sum_operator;
mean_plateau = mean(plateau(cfs_ch(3): cfs_ch(8)));
f_wavelets = f_wavelets/mean_plateau;


if npa.flag.use_ck
% power normalization
% use ck correction to transform the Cauchy wavelets to power_wavelets
ck = (f_wavelets.*f_wavelets)*sum_operator;
ck = ck.^-0.5;
ck(find(ck>5)) = 5;
ck(1)=1;
ck_mat = repmat(ck,1,jmax)*1/sqrt(2);
else 
    ck_mat = 1;
end
fwave = f_wavelets.*ck_mat; 
% fwave are power normalized wavelets in frequency space.

% reduce wavelets to jmax_final, Eliminate first and last wavelet.
fwave = fwave(:,2:end-1);
cfs = cfs(2:end-1);

%With the addition of zeros the real part of wave is the real wavelet
%and the imag part the imag wavelet.
%disp(size(fwave));
wave = ifft([fwave;zeros(le/2,jmax_final)])*2;

%shift wavelet to center.
wave = [wave(le/2+1:le,:);wave(1:le/2,:)];
% To subtract mean activate next line. 
% This makes a small difference for Morlet wavelets.
wave = wave - repmat(mean(wave,1),size(wave,1),1);
end %of function Tscharner_wavelet.

function[res] = resolution(wpa,wavelets,fwavelets,npa)
%compute time resolution
wa = real(wavelets);
modulus = wpa.modulus;
jmax = wpa.nr_of_wavelets;
sr = wpa.new_sampling_rate;
sum_operator = ones(jmax,1);
e = exp(1+i*0);

for j = 1:jmax
    sig = wa(:,j)';
    length_wavelet = size(sig,2);
    sig = [sig*0 sig sig*0];
    %do not use gauss filter to compute time_resolution, set control = 2.
    
    w_sig = wavelet_transform_V70419(sig,...
        wpa.original_sampling_rate,wpa.new_sampling_rate,wpa.scale,2);
    total_intensity_sig = w_sig*sum_operator;
    limit = max(total_intensity_sig)/e;
    k = floor(length_wavelet/modulus);
    le2 = 2*length_wavelet;
    while (total_intensity_sig(k) < limit) & (k < le2)
        k = k + 1;
    end
    low_end = k;
    while (total_intensity_sig(k) >= limit) & (k < le2)
     k = k + 1;
    end
    high_end = k - 1;
    time_resolution(j) = high_end - low_end;
end
res.tr_ch = time_resolution;
res.tr_ms = time_resolution/sr*1000;

%compute frequency resolution.
% df = 1/T and T is length of wavelet times dt.
df = wpa.original_sampling_rate/size(wa,1);
for j = 1:jmax
    k = 1;
    fw = fwavelets(:,j);
    limit = max(fw)/e;
    while fw(k) < limit
        k = k+1;
    end
    low_freq = k;
    while fw(k) > limit;
        k = k+1;
    end
    high_freq = k-1;
    fres(j) = high_freq - low_freq;
end
res.fr_ch = fres;
res.fr_Hz = fres*df;

%define Gauss filters.
le_gauss = 2 * floor(3*res.tr_ms/1000 * wpa.new_sampling_rate/2);
res.g = (le_gauss(1) - le_gauss)/2;
x = le_gauss(1);
for j = 1:jmax
    %the (x/2 +1) is necessary to compensate for position.
    res.gauss(:,j) = dnorm(1:x,x/2 +1,res.tr_ch(j) * wpa.width_factor_gauss);
end

end %of function resolution

function [distribution] = dnorm(x,mean_pos,sigma)
%computs a gaussion or normal distribution around
%the mean_pos with a width sigma. 
%x reprsesents a vector of the values where the
%distribution is computed.

u = mean_pos;
gr1 = 1/(sqrt(2*pi)*sigma);
gr2 = -1/(2*sigma*sigma);
gr3 = (x - mean_pos).*(x - mean_pos)*gr2;
distribution = gr1*exp(gr3);
end %of function dnorm.


