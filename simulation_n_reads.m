function [results ] = simulation_n_reads(modification_atlas)
% Runs simulation of ChIP-seq experament and estimation model, and inspects
% the models preformance. 
%     This function simulates the effect of ChIP-seq depth on preformance
%     of the estimation. It runs the simulation with a set appoptosis
%     vector and a gradient of number of reads, and plots the score of the
%     simulation as a fuction of the number of reads. More infermation on
%     the simulation can be found in documentation of the 'run_simulation'
%     fucntion.  
        
%   indices of blood cells in the atlas. Taken from the roadmap project
%   (http://egg2.wustl.edu/roadmap/web_portal/meta.html)
    BLD_C = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 50 51 62];
    
    % set the fraction of one non-blood cell in the appoptosis vector. the
    % rest of non-blood cells are set to 0, and the fraction of all blood
    % cells is distributed randomly and sum up to 1-fraction.

    % gradient of cell fractions to simulate on
    C_FRACS = [0.001 0.005 0.01 0.04 0.07 0.1];
    N_READS = [100000 500000 1000000 2000000 3000000 4000000 5000000 6000000 ...
        7000000 8000000 9000000 10000000];
    tiss_atlas = get_tiss_atlas(['../mod_atlases/' modification_atlas]); 
    NUM_TISS = size(tiss_atlas,2); 
    CELLS = 1:NUM_TISS;
    NUM_RPT = 3; % number of repets of the simulation
        
    % results off the simulations. last dimention holds initial appoptosis
    % vector and estimated appoptosis vector
    results = zeros(NUM_RPT, length(C_FRACS), length(N_READS), length(CELLS), 2, length(CELLS));
    
    tic
    for rp=1:NUM_RPT
        for cf=1:length(C_FRACS)
            CELL_FRAC = C_FRACS(cf);

            for c_idx=1:NUM_TISS 
                % the cell is a blood cell 
                if any(BLD_C == c_idx)
                    continue
                end
                
                % initialize random vector representing the proportion of
                % appoptosis of the blood cells
                bld_vec = randi([0, 1000], 1, length(BLD_C));
                app_vec = zeros(1,NUM_TISS);
                app_vec(c_idx) = CELL_FRAC;
                % normilize blood cell 
                bld_vec = (bld_vec * (1 - CELL_FRAC))/sum(bld_vec);

                % set values of blood cells in the appoptosis vector
                for k=1:length(BLD_C)
                    app_vec(BLD_C(k)) = bld_vec(k);
                end

                % at this point, the app_vec is a normilized vector that
                % contains the fraction value for the current cell, random
                % values for blood cells, and zeros for all other cells.
                
                for nr=1:length(N_READS)
                    [rp, cf, c_idx, nr]
                    n_reads = N_READS(nr);
                    ChIP_data = bld_chp_simulation(n_reads, app_vec, tiss_atlas);
                    est_app_vec = estimate_app_prop(ChIP_data, tiss_atlas);
                    % save initial and estimated appoptosis vector of simulation 
                    results(rp,cf,nr,c_idx,1,:) = app_vec;
                    results(rp,cf,nr,c_idx,2,:) = est_app_vec'; % I transposed est_app_vec. make sure doesnt cause problems
                    
% 
%                     figure(est_vs_int_f);
%                     subplot(ceil(length(N_READS)/2), 2, nr);
%                     non_bld_c = app_vec ~=0; % scatter plot only non zero values (to avoid clutter);
%                     scatter(CELLS(non_bld_c), app_vec(non_bld_c), 'filled');
%                     hold on;
%                     scatter(CELLS, app_vec - est_app_vec');
%                     refline([0 0]);
%                     hold off;
%                     ylim([-max(est_app_vec)/50 max(est_app_vec)/50]); 
%                     title(['number of reads ' num2str(n_reads)], 'FontSize', 10);
% 
%                     figure(score_f);
%                     scatter(n_reads, est_app_vec(c_idx));
                end
            end
        end
    end
    toc
%     figure(score_f);
%     min_frc = N_READS(1);
%     max_frc = N_READS(length(N_READS));
%     xlim([-min_frc max_frc+min_frc]);
%     ylim([-CELL_FRAC*6 CELL_FRAC*6]);
%     hline = refline([0 CELL_FRAC]);
%     hline.Color = 'r';
%     legend(hline, 'initial fraction');
%     xlabel('number of sequencing reads', 'FontSize', 14, 'FontWeight', 'bold');
%     ylabel('estimated cell fraction', 'FontSize', 14, 'FontWeight', 'bold');
%     title(['Estimated cells fraction using ' MODIFICATOIN], 'FontSize', 17, 'FontWeight', 'bold'); 
%     hold off;
%     
%     figure(est_vs_int_f);
% %     text(0.5, 0, 'initial and estimated fraction of cells', 'FontSize', 14', ...
% %         'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' );
%     xlabel('cell index');
%     ylabel('fraction');
%     legend('int frac', 'int frac - est frac');
%     hold off;
%     
%     norm1 = norm1/length(NUM_TISS);
%     norm2 = norm2/length(NUM_TISS);
%     error_frac = error_frac/length(NUM_TISS);
%     est_frctns/length(NUM_TISS);
%     figure();
%     plot(N_READS, est_frctns/length(NUM_TISS), '.');
%     hold on;
%     hline = refline([0 CELL_FRAC]);
%     hline.Color = 'r';
%     %plot(NUM_READS, std(est_frctns)); % need to save all values of runs
%     %for this
%     hold off;
%     title(['average estimated cell fraction. initial fraction = ' num2str(CELL_FRAC)]);
%     xlabel('number of sequencing reads', 'FontSize', 14, 'FontWeight', 'bold');
%     ylabel('estimated cell fraction', 'FontSize', 14, 'FontWeight', 'bold');
%     
%     figure()
%     subplot(1,3,1);
%     plot (N_READS, norm1, '-o');
%     title('norm 1', 'FontSize', 14, 'FontWeight', 'bold');
%     xlabel('number of reads', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('score', 'FontSize', 12, 'FontWeight', 'bold');
% 
%     subplot(1,3,2);
%     plot(N_READS, norm2, '-o');
%     title('norm 2', 'FontSize', 14, 'FontWeight', 'bold');
%     xlabel('number of reads', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('score', 'FontSize', 12, 'FontWeight', 'bold');
% 
%     subplot(1,3,3);
%     plot(N_READS, error_frac, '-o');
%     title('error fraction', 'FontSize', 14, 'FontWeight', 'bold');
%     xlabel('number of reads', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('score', 'FontSize', 12, 'FontWeight', 'bold');
%     suptitle(['Score as function of # reads (average over cells). Initial cell fraction = ' num2str(CELL_FRAC)]); % 'FontSize', 17, 'FontWeight', 'bold');
res_file = [modification '_results'];
%save(res_file);
end
   