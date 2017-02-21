function [est_app_prop] = run_simulation( n_reads )
% Runs simulation of ChIP-seq experament and estimation model, and inspects
% the models preformance. 
%     This function creates different vectors representing appoptosis
%     proportion of the tissues, then simulates the ChIP-seq exeprament by 
%     running bld_chp_simulation.m with that vector and estimate_app_prop.m
%     to get the models estimate of the appoptosis proportion of the
%     tissues. It then compares the original appoptosis proportion and the
%     estimated appoptosis proportion to receive a feedback of the model
%     preformance. The simulation recieves the number of reads per
%     experament. 
        
%   indices of blood cells in the atlas. Taken from the roadmap project
%   (http://egg2.wustl.edu/roadmap/web_portal/meta.html)
    BLD_C = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 50 51 62];
    INT_FRAC = 0.1;
    MIN_FRAC = 0.0001;
    
    tiss_atlas = get_tiss_atlas('../mod_atlases/H3K27ac_atlas');  % make sure that matrix is passed by 
            %reference when calling bld_chp_simulation and estimate_app_prop 
    NUM_TISS = size(tiss_atlas,2); 
    
    %   go through all non-blood cells and for each cell, set a gradient of
    %   decreasing fractions of that cell. plot the score of the simulation
    %   as a function of the cell fraction. 
    for cell=1:1 % change to 1:NUM_TISS
        % the cell is a blood cell 
        if any(BLD_C == cell)
            continue
        end
        
        % initialize random vector representing the proportion of
        % appoptosis of the blood cells
        bld_vec = randi([0, 1000], 1, length(BLD_C));
        
        % loop from initial fraction to minimum fraction in steps the size
        % of the ratio between the two.
        fractions=INT_FRAC:-(INT_FRAC-MIN_FRAC)/(INT_FRAC/MIN_FRAC)*100:MIN_FRAC;
        norm1 = zeros(1,length(fractions));
        norm2 = zeros(1,length(fractions));
        norm_nir = zeros(1,length(fractions)); % change name
        est_frctns = zeros(1,length(fractions)); % vector of the franctions
        % that were estimated. Reflects how well was the model in finding 
        % minor fraction among majority of blood cells. 
        for i=1:2  % change back to 1:length(fractions)
            tic
            frctn = fractions(i);
            % initialize appoptosis vector 
            app_vec = zeros(1,NUM_TISS);
            app_vec(i) = frctn;
            % normilize blood cell 
            bld_vec = (bld_vec * (1 - frctn))/sum(bld_vec);
            
            % set values of blood cells in the appoptosis vector
            for k=1:length(BLD_C)
                app_vec(BLD_C(k)) = bld_vec(k);
            end
            
            % at this point, the app_vec is a normilized vector that
            % contains the fraction value for the current cell, random
            % values for blood cells, and zero for all other cells.
            
            ChIP_data = bld_chp_simulation(n_reads, app_vec, tiss_atlas);
            est_app_prop = estimate_app_prop(ChIP_data, tiss_atlas);
            norm1(i) = norm(app_vec - est_app_prop', 1);
            norm2(i) = norm(app_vec - est_app_prop');
            non_z = find(est_app_prop'); % non zero elements of estimation vector
            norm_nir(i) = sum(abs((app_vec(non_z)./est_app_prop(non_z)') -1).^2);
            est_frctns(i) = est_app_prop(cell); % saving the value of estimated
            % fraction of the non blood-cell.
            toc
                       
        end
        
        %norm1 = log(norm1)
        %norm2 = log(norm2)
        %norm_nir = log(norm_nir)
        norm1;
        norm2;
        norm_nir
        abs(fractions./est_frctns-ones(1,length(fractions))).^2
        clf
        figure()
        title(['score as function of fraction: ' 'cell ' num2str(i) ' ' num2str(n_reads) ' reads']);
 
        subplot(1,3,1);
        plot (fractions, norm1);
        title('norm1');
        xlabel('appoptosis fraction');
        ylabel('score');
        
        subplot(1,3,2);
        plot(fractions, norm2);
        title('norm2');
        xlabel('appoptosis fraction');
        ylabel('score');
        
        subplot(1,3,3);
        plot(fractions, norm_nir);
        title('norm nir');
        
        % title(['score as function of fraction: ' 'cell ' num2str(i) ' ' num2str(n_reads) ' reads']);
        xlabel('appoptosis fraction');
        ylabel('score');
        %legend('norm1', 'norm2', 'norm nir');
   
    end
end
   

        
        
        
        
        
        
        
        
        
%         
%     N_RUNS = 9;
%     runs=zeros(N_RUNS,2);
%     % Clear current figure window
%     clf
%     figure();
%     color_vec = hsv(N_RUNS);
%     
%     hold on;
%     for i=1:N_RUNS
%         tic
%         i
%         
%         % Create a 1-by-NUM_TISS array of random integer values drawn from a
%         % discrete uniform distribution. randi([imin imax]).
%         rnd_app = randi([0, 1000], 1, NUM_TISS); % TODO check about the values of imin, imax
%         rnd_app = rnd_app/sum(rnd_app);
%                 
%         % Create a 1-by-NUM_TISS array of random integer values drawn from a
%         % exponential distribution with parameter mu. 
%         mu = 10;
%         ex_app = exprnd(mu, 1, NUM_TISS);
%         ex_app = ex_app/sum(ex_app);
%         
%         
%         ChIP_data = bld_chp_simulation(n_reads, ex_app, tiss_atlas);
%         est_app_prop = estimate_app_prop(ChIP_data, tiss_atlas);
%         
%         plot(abs(ex_app-est_app_prop'), 'Color', color_vec(i,:));
%         hline = refline([0 min(ex_app)]);
%         hline = refline([0 max(ex_app)]);
%         
% %   calculate the norm2 and norm inf of the appotosis and estimated
% %   appotosis rate. 
%         nrm = norm(ex_app-est_app_prop'); 
%         nrm_inf = norm(ex_app-est_app_prop', Inf);
%         runs(i,1) = nrm;
%         runs(i,2) = nrm_inf;
%         toc
%     end
%     hold off;