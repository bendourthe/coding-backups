function [ GemFreq ] = get_mean_freq(magnitude_harmonics,freq_harmonics)
%GET_MEAN_FREQ Calculates the average frequency
%   Input Arguments:
%       1. magnutide_harmonics: this is a vector with the magnitude
%       (amplitude) of the harmonics
%       2. freq_harmonics: is a vector with the frequency of the harmonics.
%   Output Argument:
%       1. GemFreq: The average frequency 
%
%   --------Author: Maarten Afschrift (14/03/2014)-------

GemFreq = sum(magnitude_harmonics.*freq_harmonics)./sum(magnitude_harmonics);
end

