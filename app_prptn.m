function [ app_prptn ] = app_prptn( n_tiss )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%   random appoptosis proportion of tissues. 
    [rnd_vec] = rand(1,n_tiss);
    app_prptn = rnd_vec/sum(rnd_vec); 

end

