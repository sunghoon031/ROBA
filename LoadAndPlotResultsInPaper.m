clc; close all;

% Add all subfolders.
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Synthetic dataset results

% result_mat is [11 x nRuns] matrix, consisting of
% - row 1: edge percent
% - row 2: mean rotation error of edges 
% - row 3: median rotation error of edges
% - row 4: mn1 of RA
% - row 5: mn1 of RA + ROBA
% - row 6: md1 of RA
% - row 7: md1 of RA + ROBA
% - row 8: mn2 of RA
% - row 9: mn2 of RA + ROBA
% - row 10: md2 of RA
% - row 11: md2 of RA + ROBA

figure;
for config = 1:12
    fname = sprintf('synthetic_results_%d.mat', config);
    load(['results\', fname]);

    med_ep = median(result_mat(1,:));
    med_mnee = median(result_mat(2,:));
    med_mdee = median(result_mat(3,:));
    med_RA = median(result_mat(4,:));
    med_ROBA = median(result_mat(5,:));

    error_reduced = result_mat(5,:) - result_mat(4,:);
    error_reduced = sum(error_reduced < 0);

    disp(['config ', num2str(config), ...
        ', md edge: ', num2str(med_ep), ...
        '%, md mn edge error: ', num2str(med_mnee), ...
        'deg, md md edge error: ', num2str(med_mdee), 'deg'])

    disp(['config ', num2str(config), ...
        ', md RA: ', num2str(med_RA), ...
        'deg, md ROBA: ', num2str(med_ROBA), ...
          'deg, error reduced: ', num2str(error_reduced), '%'])



    subplot(2,6,config)


    boxplot([result_mat(4,:)', result_mat(5,:)'], 'Widths', 0.8)
%         boxplot([result_mat(6,:)', result_mat(7,:)'], 'Widths', 0.8)
%         boxplot([result_mat(8,:)', result_mat(9,:)'], 'Widths', 0.8)
%         boxplot([result_mat(10,:)', result_mat(11,:)'], 'Widths', 0.8)
    grid on
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    set(gca,'XTickLabel',{'RA', 'RA + ROBA'});

    switch config
        case 1
            title('Baseline')
        case 2
            title('More points')
        case 3
            title('Fewer views')
        case 4
            title('More views')
        case 5
            title('Closer points')
        case 6
            title('Farther points')
        case 7
            title('Less noise')
        case 8
            title('More noise')
        case 9
            title('Planar scene')
        case 10
            title('Pure rotations')
        case 11
            title('Pure + Planar')
        case 12
            title('Mixed rotations')
    end
    
end

load synthetic_results_baseline_evolution.mat;
figure;
subplot(1,2,1)
boxplot(baseline_cost_ROBA,  'Widths', 1, 'Whisker', 0, 'symbol','')
ylim([0 6])
ylabel('Total cost')
xlabel('Number of iterations')
set(gca, 'YGrid', 'on', 'XGrid', 'off')

subplot(1,2,2)
boxplot(baseline_mn1_ROBA, 'Widths', 1, 'Whisker', 0, 'symbol','')
ylim([0 3])
ylabel('Mean error (deg)')
xlabel('Number of iterations')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
suptitle('Baseline')    

%% Real data results
load real_results_main.mat;

for i = 1:15

    figure
    switch i
        case 1
            suptitle('Alamo')
        case 2
            suptitle('Ellis Island')
        case 3
            suptitle('Gendarmenmarkt')
        case 4
            suptitle('Madrid Metropolis')
        case 5
            suptitle('Montreal Notre Dame')
        case 6
            suptitle('Notre Dame')
        case 7
            suptitle('NYC Library')
        case 8
            suptitle('Piazza del Popolo')
        case 9
            suptitle('Piccadilly')
        case 10
            suptitle('Roman Forum')
        case 11
            suptitle('Tower of London')
        case 12
            suptitle('Trafalgar')
        case 13
            suptitle('Union Square')
        case 14
            suptitle('Vienna Cathedral')
        case 15
            suptitle('Yorkminster')
    end
    
    
    subplot(1,5,1)
    plot(0:100, ROBA_cost(i,:))
    axis square
    grid on
    set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
    xticks(0:10:100)
    xlabel('Iterations')
    title('Total cost')

    subplot(1,5,2)
    plot(0:100, ROBA_mean_L1(i,:))
    axis square
    grid on
    set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
    xticks(0:10:100)
    xlabel('Iterations')
    ylabel('Error (deg)')
    title('mn1')
    

    subplot(1,5,3)
    plot(0:100, ROBA_med_L1(i,:))
    axis square
    grid on
    set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
    xticks(0:10:100)
    xlabel('Iterations')
    ylabel('Error (deg)')
    title('md1')

    subplot(1,5,4)
    plot(0:100, ROBA_mean_L2(i,:))
    axis square
    grid on
    set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
    xticks(0:10:100)
    xlabel('Iterations')
    ylabel('Error (deg)')
    title('mn2')

    subplot(1,5,5)
    plot(0:100, ROBA_med_L2(i,:))
    axis square
    grid on
    set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
    xticks(0:10:100)
    xlabel('Iterations')
    ylabel('Error (deg)')
    title('md2')

    
end


%% Relative errors of the aggregated results from the real-world datasets. 
ROBA_mean_L1_normalized = ROBA_mean_L1;
ROBA_med_L1_normalized = ROBA_med_L1;
ROBA_mean_L2_normalized = ROBA_mean_L2;
ROBA_med_L2_normalized = ROBA_med_L2;

for i = 1:15
    ROBA_mean_L1_normalized(i,:) = ROBA_mean_L1_normalized(i,:)/ROBA_mean_L1_normalized(i,1)*100;
    ROBA_med_L1_normalized(i,:) = ROBA_med_L1_normalized(i,:)/ROBA_med_L1_normalized(i,1)*100;
    ROBA_mean_L2_normalized(i,:) = ROBA_mean_L2_normalized(i,:)/ROBA_mean_L2_normalized(i,1)*100;
    ROBA_med_L2_normalized(i,:) = ROBA_med_L2_normalized(i,:)/ROBA_med_L2_normalized(i,1)*100;
end

figure;
subplot(1,4,1)
boxplot(ROBA_mean_L1_normalized(:,11:10:101), 'whisker', 1)
ylim([0 100])
xticklabels({'10','20','30', '40', '50', '60', '70', '80', '90', '100'})
xlabel('Number of iterations')
ylabel('Relative error (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
axis square
title('mn1')

subplot(1,4,2)
boxplot(ROBA_med_L1_normalized(:,11:10:101), 'whisker', 0.5)
ylim([0 100])
xticklabels({'10','20','30', '40', '50', '60', '70', '80', '90', '100'})
xlabel('Number of iterations')
ylabel('Relative error (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
axis square
title('md1')

subplot(1,4,3)
boxplot(ROBA_mean_L2_normalized(:,11:10:101), 'whisker', 1)
ylim([0 100])
xticklabels({'10','20','30', '40', '50', '60', '70', '80', '90', '100'})
xlabel('Number of iterations')
ylabel('Relative error (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
axis square
title('mn2')

subplot(1,4,4)
boxplot(ROBA_med_L2_normalized(:,11:10:101), 'whisker', 1)
ylim([0 100])
xticklabels({'10','20','30', '40', '50', '60', '70', '80', '90', '100'})
xlabel('Number of iterations')
ylabel('Relative error (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
axis square
title('md2')




