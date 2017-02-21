function [ tiss_atlas ] = get_tiss_atlas ( atlas_path )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%     % default path to atlas if no path was given 
%     if ~nargin
%         atlas_path = 'new_atlas_matrix';
%     end
    
    % rows of atlas = genomic window, coulmns of matrix = read count of
    % samples. Matrix[i][j] = coverage fraction of window at genomic site i,
    % in sample j. 
    atlas_mat = load(atlas_path);
    atlas_mat = atlas_mat.ta;
    % remove first 3 coulmns that are the chr, start cordination and end
    % cordination (remain matrix contains only readsXsamples).
    tiss_atlas = atlas_mat(:,4:size(atlas_mat,2));
    
    
% the next block of code was the function implementation for the matrix that
% was created from the bam files (using the function 'hist_analysis2.m). 
% the format of that matrix is described below. 

% % first row of window cordinations and first row of sample index. 
% % rows = samples, columns = windows. First row is window cordination,
%     % first column is sample index. For each (i>1, j>1): Atlas_mat[i,j] = 
%     % number of reads in the i'th sample in the j'th window.
%     atlas_mat = importdata(atlas_path);
%     % m = number of samples
%     n_smpls = size(atlas_mat,1) - 1;
%     % n = number of windows
%     n_win = size(atlas_mat, 2) - 1;
% 
%     % cut of matrix first row and column
%     tiss_atlas = atlas_mat(2:n_smpls + 1, 2:n_win + 1);

end

