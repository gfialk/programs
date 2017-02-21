function [ mod_atlas, tiss_list ] = win_file_to_atlas ( modification )
% function that generates atlas of a given modification
%   The function receives a specific modificatoin directory, goes through all
%   peak files of that modification and creates a atlas (tissues X genomic
%   window) where atlas[i][j] contains the fraction of the i'th genomic
%   window of the j'th tissue that was covered by a peak in the respective
%   peak file. The rows of the atlas matrix represent the different tissues
%   and the columns of the matrix represent the genomic windows.  
    
    WIN_SIZE = 10000;
    
    % number of genomic windows (number of lines in each win file)
    NUM_WIN = 309581;
    
    NUM_READ_FLD = 4;
    base_dir = '../win_files';
    mdf_dir = dir([base_dir '/' modification '/*.win']);
    
    num_tiss = size(mdf_dir,1);
    tiss_list = cell(1, num_tiss);
    tiss_inx = 0;
    
    % create matrix WIN_NUMBER X TISS_NUMBER. add 3 coulmns for chr,
    % window start cordination and window end cordination
    tiss_atlas = zeros(NUM_WIN, num_tiss + 3); 
    tic
    for win_file = mdf_dir'
        file_name = [base_dir '/' modification '/' win_file.name];
        tiss_inx = tiss_inx + 1
        tiss_list{tiss_inx} = file_name;
        w_file = importdata(file_name);
        
        % initialize atlas whith first three coulmns (chr, start_cord,
        % end_cord). coulmns 4-end are the read counts of the samples at
        % that cordination
        if tiss_inx == 1
            tiss_atlas(:, 1:4) = w_file;
        
        else
            tiss_atlas(:, tiss_inx + 3) = w_file(:,NUM_READ_FLD);
        end
        
    end
    toc
    
    mod_atlas = tiss_atlas;
    atlas_name = ['../mod_atlases/' modification '_atlas.mat']; 
    %dlmwrite(atlas_name, mod_atlas, '\t');
    save(atlas_name, 'mod_atlas');
    
              
end