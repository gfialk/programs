function [ output_args ] = plot_results ( result_file )
%plot results of simulation
%   This function receives a results - a cell containg the results of the
%   simulation and plots the following:
%   1) scores (norm 1, norm 2, error fraction) as function of coverage
%   2) scores as function of initial cell fraction
%   3) cell estimation fraction as function of coverage (compared to initial fraction?)
%   4) initial/estimation ratio(difference?) per all cells
%   The input param results R, is a MxNxO cell array, where R[i][j][k] holds the
%   information of the simulation with the i'th cell fraction and the j'th
%   number of reads in the k'th repet. R[i][j][k] holds 'A' - a KxL cell array where
%   A[1][c] hold the initial appoptosis vector of the simulation where the
%   appoptosis rate of the c's cell was set. A[2][c] holds the estimated
%   appoptosis vector of that simulation.
%   Params:
%   n_rpt - double representing number of repets of the simulation
%   c_fract - vector of cell fractions gradient used in the simulatiom
%   cells - vector of indiceis of cells used in simulation
%   n_reads - vector of number of reads gradient used in the simulation

load(result_file);
MODIFICATION = modification;
N_FRC = length(C_FRACS);
N_CLS = length(CELLS);
N_RDS = length(N_READS);

figure();
for frc=1:N_FRC

    % store norm values 
    norms = zeros(3,N_RDS,N_CLS);
    
    for nr=1:N_RDS
        init_vecs = squeeze(mean(results(:,frc,nr,:,1,:),1));
        est_vecs = squeeze(mean(results(:,frc,nr,:,2,:),1));
        
        for cidx=1:N_CLS
            init_vec = init_vecs(cidx,:);
            est_vec = est_vecs(cidx,:);
            
            if ~any(init_vec)
                continue
            end
            
            norm1 = norm(init_vec - est_vec, 1);
            norm2 = norm(init_vec - est_vec);
            non_z = find(init_vec);
            err_frc = 1/length(non_z) * (sum(abs(init_vec(non_z)./est_vec(non_z)-1)));
            norms(1,nr,cidx) = norm1; norms(2,nr,cidx) = norm2; norms(3,nr,cidx) = err_frc; 
        end
    end
            
    subplot(1,3,1); hold on; errorbar(n_reads,mean(norms(1,:,:),3), std(norms(1,:,:),0,3)/2); hold off;
    subplot(1,3,2); hold on; errorbar(n_reads,mean(norms(2,:,:),3), std(norms(2,:,:),0,3)/2); hold off;
    subplot(1,3,3); hold on; errorbar(n_reads,mean(norms(3,:,:),3), std(norms(3,:,:),0,3)/2); hold off;
end

subplot(1,3,1);
xlabel('number of reads'); ylabel('score (norm 1)');
legend(strsplit(num2str(c_fracs)));
subplot(1,3,2);
xlabel('number of reads'); ylabel('score (norm 2)');
legend(strsplit(num2str(c_fracs)));
subplot(1,3,3);
xlabel('number of reads'); ylabel('score (error fraction)');
legend(strsplit(num2str(c_fracs)));

suptitle(['ChIP simulation with antibody for ' MODIFICATION]);


% initial fraction vs. estimation fraction as function of coverage 
    figure()
    for nr=1:12
        subplot(4,3,nr)
        xlabel('estimated fraction'); ylabel('initial fraction');
        refline([1 0]);
        title(sprintf('%d reads', N_READS(nr)));
        int_vecs = squeeze(mean(results(:,:,nr,:,1,:),1));
        est_vecs = squeeze(mean(results(:,:,nr,:,2,:),1));
        hold on;
        for cidx=1:98
            int_vec = squeeze(int_vecs(:,cidx,cidx))';
            est_vec = squeeze(est_vecs(:,cidx,cidx))';
            scatter(int_vec, est_vec)
        end
        hold off;
        suptitle('initial vs. estimated cell fraction');
    end
    
end

