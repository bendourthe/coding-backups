%%

data=EMG_16(:,1);
n_data=length(data);
step=round(n_data/51);
vector_step=round(1:step:n_data);


for i=1:length(vector_step)-1    
        [a b]=pwelch(data(vector_step(i):vector_step(i+1)),[],[],500,1000);
        % Determine meanfrequency.
        Gemfreq(:,i) = sum(a.*b)./sum(a);
end
