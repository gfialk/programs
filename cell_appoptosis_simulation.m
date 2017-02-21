function [ bld_smpl ] = cell_appoptosis_simulation ( n_cells, app_prop, tiss_atlas )
%bld_chp_simulation generates simulated ChIP-seq data
%     Let m = number of tissues in atlas
%     n = genome size (number of windows) 
%     app_prop = appoptosis proportion of every tissue (vector 1Xm)
%     n_rds = number of estimated reads of the ChIP experament (simulating sequencing detph)
%     
%     This function simulates amount of modified nucleosomes in blood serum by computing
%     blood_smpl = atlas * app_prop, and then simulates the ChIP experament by sampling
%     from a poisson distribution of (normilized) blood_smpl. The function 
%     returns a vector (1Xn) which contains read count of modified nucleosome 
%     in every genomec window (representing a ChIP-seq experament data)
    
    n_cl_typs = size(tiss_atlas,2); % number of cells
    n_win = size(tiss_atlas,1); % number of genomic windows
    atls_vec = reshape(tiss_atlas,[],1); % turn atlas matrix to a column vector
    prob_vec = atls_vec;
    for c_typ=1:n_cl_typs
        cell_frc = app_prop(c_typ);
        prob_vec(c_typ*((c_typ-1)*n_win)+1:c_typ*n_win) = prob_vec(c_typ*((c_typ-1)*n_win)+1:c_typ*n_win)*cell_frc/n_win;
    end
    lambdas = (prob_vec/sum(prob_vec))*(n_win*n_cells);
    bld_smpl = poissrnd(lambdas);
    bld_smpl = reshape(bld_smpl,n_win,n_cl_typs);
    bld_smpl = sum(bld_smpl,2);
   % bld_smpl = bld_smpl/n_cells;
    
end