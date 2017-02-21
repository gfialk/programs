close all;
figure()
for nr=1:12
    subplot(4,3,nr)
    xlabel('estimated fraction'); ylabel('initial fraction');
    xlim([0 0.12]); ylim([0 0.12]);
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