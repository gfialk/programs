function [bld_smpl] = gnrt_bld_smpl(atls_mat, apptss_rate)
% Generate "cfDNA (modified histones)" by simulating cell appoptosis. 
%     Function generates "blood sample data" using histone modification atlas 
%     and cell appoptosis rate for every cell. The function receives path to 
%     atlas matrix (AM) - that contains histone modifications rate for tissues
%     throughout the genome and apopptosis rate vector (arv) and computes 
%     BS = arv * AM. 
%     params: 
%     atls_mat - path to atlas matrix (m+1 X n+1). Rows(AM) = tissue, Col(AM) = 
%     genomic window (more documentation can be found at 'mat_analysis.m').
%     apptss_rate - 1Xm vector representing appotosis rate of m tissues
%     (normilized).
%     returns:
%     BS = 1Xm vector representing blood sample (modified histones from m tissues 

    % flag - is the matrix a binary matrix
    is_bin = 0;

    % rows = samples, columns = windows. First row is window cordination,
    % first column is sample index. For each (i>1, j>1): Atlas_mat[i,j] = 
    % number of reads in the i'th sample in the j'th window.
    atlas_mat = importdata(atls_mat);

    % m = number of samples
    n_smpls = size(atlas_mat,1) - 1;
    % n = number of windows
    n_win = size(atlas_mat, 2) - 1;

    % cut of matrix first row and column
    atlas_mat = atlas_mat(2:n_smpls + 1, 2:n_win + 1);
    
    
    % generate random observation vector 1Xm (number of reads in every window).
    bld_smpl = apptss_rate * atlas_mat;

end