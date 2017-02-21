function [ ChIP_data ] = bld_chp_simulation ( n_rds, app_prop, tiss_atlas )
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
     
    % tiss_atlas = get_tiss_atlas(); 
    bld_smpl =  tiss_atlas * app_prop';
    lambdas = (bld_smpl/sum(bld_smpl))*n_rds;
    ChIP_data = poissrnd(lambdas);

end

