function [ Gemfreq] = average_freq_emg( data,n_steps )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


n_data=length(data);
step=round(n_data/n_steps);
vector_step=round(1:step:n_data);

Gemfreq=zeros(n_steps-1,1);

for i=1:length(vector_step)-1    
        [a b]=pwelch(data(vector_step(i):vector_step(i+1)),[],[],500,2000);
        % Determine meanfrequency.
        Gemfreq(i) = sum(a.*b)./sum(a);
end

end

