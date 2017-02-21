clear all;

% flag - is the matrix a binary matrix
is_bin = 1;

CORD_LINE = 1;
SAMPLE_INDEX = 1;
atlas_mat = importdata('bin_atlas_matrix');
% convert to matrix with lines representing samples and coulmns
% representing DNA cordination. 
%atlas_mat = atlas_mat';

%first line of matrix is cordination. first coulmn is sample index.
num_samples = size(atlas_mat,1) - 1;

% number of windows with peaks in sample - for binary matrix only
if is_bin
    n_pks_smpl = zeros(1, num_samples);
end

% iterate through samples (first line is cordinations)
for i=1:num_samples

    % plot samples (visualization for non-binary peaks)
    len_mat = length(atlas_mat);
    loc = atlas_mat(CORD_LINE, 2:len_mat);
    reads = atlas_mat(i+1, 2:len_mat);
    if ~is_bin
        subplot(num_samples,1,i);
        plot(loc,reads,'.')
        str_title = ['sample num ' num2str(i)]; 
        title(str_title);
    end
    
    if is_bin
        n_pks_smpl(i) = sum(reads);
    end
end

if is_bin
    % find windows with peaks in all windows - only for binary virsion
    h_keeping_gns = find(all(atlas_mat(2:num_samples+1,2:len_mat)));
    
    % alternativly can define house keeping genes as windows where alpha
    % percent of samples have peak.
    %h_keeping_gns = find(sum(atlas_mat(2:num_samples+1,2:len_mat)) == num_samples*alpha);
    length(h_keeping_gns)
    % uniqe windows (peaks in only one sample. notice the offset of 1) - only for binary virsion
    uniqe = find(sum(atlas_mat(2:num_samples+1,2:len_mat)) == 1); 
    length(uniqe)

end