function [ est_app ] = estimate_app_prop( ChIP_data, tiss_atlas )
%estimate_app_prop estimates tissue appoptosis proportion
%   This function receives a simulated ChIP-seq data and estimates the
%   appoptosis proportion of every cell. Let m = number of tissues in atlas, 
%   n = genome size (number of windows). 
%   The function computes [est_app = argmin(||X-Aw||)]
%   where X (1Xn) is the simulated ChIP-seq data (read count in every 
%   genomic window), A (mXn) - is the tissue atlax, and w (1Xm) - is the 
%   appoptosis proportion of every tissue. 

    % atlas = get_tiss_atlas();
    % Solve constrained linear least-squares problems (min(w)||X-Aw||^2).
    % comment: the function parameter of zeros is the lower bound of x  
    % representing the constraint that w>=0 in every coordination (since
    % appoptosis proportion of a cell can't be negative)
    %est_app = lsqlin(tiss_atlas, ChIP_data, [], [], [], [], zeros(1,size(tiss_atlas,2)));
    
    est_app = lsqlin(tiss_atlas, ChIP_data); % using the simpler version with no 
    % constraints do to increase in run time. TODO add constraint and solve
    % run time issue. 

    % normilize resault so that it sums up to 1. 
    est_app = est_app/sum(est_app);
    
end

