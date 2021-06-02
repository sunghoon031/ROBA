clear all; close all; clc;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Choose options:
plotResult = true;
datasets = 2;%1:15; 
% 1: Alamo, 2: Ellis_Island, 3: Gendarmenmarkt
% 4: Madrid_Metropolis, 5: Montreal_Notre_Dame, 6: Notre_Dame
% 7: NYC_Library, 8: Piazza_del_Popolo, 9: Piccadilly
% 10: Roman_Forum, 11: Tower_of_London, 12: Trafalgar
% 13: Union_Square, 14: Vienna_Cathedral, 15: Yorkminster

square_rooting = true;
switch_alpha = true;
approximate_gradient = true;

%%
tStart = tic; 

dev_angle = 10^(-4);
nIterations = 100;

RA_mean_L1 = nan(15,1);
RA_med_L1 = nan(15,1);
RA_mean_L2 = nan(15,1);
RA_med_L2 = nan(15,1);

ROBA_cost = nan(15,nIterations+1);
ROBA_mean_L1 = nan(15,nIterations+1);
ROBA_med_L1 = nan(15,nIterations+1);
ROBA_mean_L2 = nan(15,nIterations+1);
ROBA_med_L2 = nan(15,nIterations+1);

times_RA = nan(15,1);
times_ROBA_init = nan(15,1);
times_ROBA_per_iteration = nan(15,1);

for dataset = datasets
    
    switch dataset
        case 1
            load Alamo.mat;
            disp(['========== ', num2str(dataset), '. Alamo =========='])
        case 2
            load Ellis_Island.mat;
            disp(['========== ', num2str(dataset), '. Ellis_Island =========='])
        case 3
            load Gendarmenmarkt.mat;
            disp(['========== ', num2str(dataset), '. Gendarmenmarkt =========='])
        case 4
            load Madrid_Metropolis.mat;
            disp(['========== ', num2str(dataset), '. Madrid_Metropolis =========='])
        case 5
            load Montreal_Notre_Dame.mat;
            disp(['========== ', num2str(dataset), '. Montreal_Notre_Dame =========='])
        case 6
            load Notre_Dame.mat;
            disp(['========== ', num2str(dataset), '. Notre_Dame =========='])
        case 7
            load NYC_Library.mat;
            disp(['========== ', num2str(dataset), '. NYC_Library =========='])
        case 8
            load Piazza_del_Popolo.mat;
            disp(['========== ', num2str(dataset), '. Piazza_del_Popolo =========='])
        case 9
            load Piccadilly.mat;
            disp(['========== ', num2str(dataset), '. Piccadilly =========='])
        case 10
            load Roman_Forum.mat;
            disp(['========== ', num2str(dataset), '. Roman_Forum =========='])
        case 11
            load Tower_of_London.mat;
            disp(['========== ', num2str(dataset), '. Tower_of_London =========='])
        case 12
            load Trafalgar.mat;
            disp(['========== ', num2str(dataset), '. Trafalgar =========='])
        case 13
            load Union_Square.mat;
            disp(['========== ', num2str(dataset), '. Union_Square =========='])
        case 14
            load Vienna_Cathedral.mat;
            disp(['========== ', num2str(dataset), '. Vienna_Cathedral =========='])
        case 15
            load Yorkminster.mat;
            disp(['========== ', num2str(dataset), '. Yorkminster =========='])
    end




    %% Rotation averaging
    tic
    R_avg_mat = AverageSO3Graph(RR,edge_IDs_reverse, 'Method', 'L0.5');
    times_RA(dataset) = toc;
    
    R_avg = cell(1,nViews);
    for i = 1:nViews
        R_avg{i} = R_avg_mat(:,:,i);
    end
    
    [~,indvidual_errors_averaging, mean_error_L1, median_error_L1] = AlignRotationL1(views.R, R_avg);
    [~, mean_error_L2, median_error_L2] = AlignRotationL2(views.R, R_avg);
    disp(['Rot Avg error (deg): mean_L1 ', num2str(mean_error_L1), ...
        ', med_L1 = ', num2str(median_error_L1), ...
        ', mean_L2 = ', num2str(mean_error_L2), ...
        ', med_L2 = ', num2str(median_error_L2), ...
        ', Took ', num2str(times_RA(dataset)), 's.'])
    RA_mean_L1(dataset) = mean_error_L1;
    RA_med_L1 = median_error_L1;
    RA_mean_L2 = mean_error_L2;
    RA_med_L2 = median_error_L2;
    
    %% Initialize ROBA
    tic;
    precomputed_mat = cell(1, nEdges);
    total_cost_gt = 0;
    total_cost = 0;
    for i = 1:size(edge_IDs, 2) 
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        [~, point_idx1, point_idx2] = intersect(views.points_ID{j}, views.points_ID{k});
        rays_j = views.rays_SIFT{j}(:,point_idx1);
        rays_k = views.rays_SIFT{k}(:,point_idx2);
        precomputed_mat{i} = GetNeighborData(rays_j, rays_k);
        total_cost = total_cost + GetPairwiseCost(R_avg{j}, R_avg{k}, precomputed_mat{i}, square_rooting);
        total_cost_gt = total_cost_gt + GetPairwiseCost(views.R{j}, views.R{k}, precomputed_mat{i}, square_rooting);
    end
    
    R_est = R_avg;
    r = nan(3*nViews,1);
    for i = 1:nViews
        r(3*i-2:3*i) = LogMap(R_est{i});
    end
    
    times_ROBA_init(dataset) = toc;
    
    disp(['Initial cost = ', num2str(total_cost), ', GT cost = ', num2str(total_cost_gt)])
    disp(['Our initialization took ', num2str(times_ROBA_init(dataset)), 's (init)'])

    %% ROBA

    min_total_cost = total_cost;
    R_best = R_est;
    best_it = 1;
    [~,rot_individual_errors_init, mean_error_L1, median_error_L1] = AlignRotationL1(views.R, R_est);
    [~, mean_error_L2, median_error_L2] = AlignRotationL2(views.R, R_est);
    
    ROBA_cost(dataset, 1) = total_cost;
    ROBA_mean_L1(dataset, 1) = mean_error_L1;
    ROBA_med_L1(dataset, 1) = median_error_L1;
    ROBA_mean_L2(dataset, 1) = mean_error_L2;
    ROBA_med_L2(dataset, 1) = median_error_L2;

    %ADAM parameters:
    alpha = 0.01; beta1 = 0.9; beta2 =0.999; epsilon = 10^(-8);   
    m_prev = 0;
    v_prev = 0;
    r_prev = zeros(1, 3*nViews);
    
    costUpCount = 0;

    tic;
    for it = 1:nIterations  
            
        total_cost_prev = total_cost;
        [total_cost, g] = ComputeCostGradient(precomputed_mat, edge_IDs, R_est, dev_angle, approximate_gradient, square_rooting);

        % ADAM:
        m = beta1*m_prev + (1-beta1)*g;
        v = beta2*v_prev + (1-beta2)*(g.^2);
        m_prev = m;
        v_prev = v;
        m_hat = m/(1-beta1^it);
        v_hat = v/(1-beta2^it);
        r = r - alpha*m_hat./(sqrt(v_hat)+epsilon);

        for i = 1:nViews
            R_est{i} = ExpMap(r(3*(i-1)+1:3*(i-1)+3));
        end

        if (switch_alpha)
            if (total_cost > total_cost_prev)
                costUpCount = costUpCount + 1;
                if (costUpCount == 5)
                    alpha = 0.001;
                    costUpCount = 0;
                end
            else
                costUpCount = 0;
            end
        end

        total_cost_prev = total_cost;

        if (total_cost < min_total_cost)
            min_total_cost = total_cost;
            R_best = R_est;
            best_it = it+1;
        end
        
        ROBA_cost(dataset, it+1) = total_cost;
        [~,individual_errors, mean_error_L1, median_error_L1] = AlignRotationL1(views.R, R_est);
        [~, mean_error_L2, median_error_L2] = AlignRotationL2(views.R, R_est);
        
        ROBA_mean_L1(dataset, it+1) = mean_error_L1;
        ROBA_med_L1(dataset, it+1) = median_error_L1;
        ROBA_mean_L2(dataset, it+1) = mean_error_L2;
        ROBA_med_L2(dataset, it+1) = median_error_L2;
        
        if (mod(it,10)==0)
            disp(['ROBA (it=', num2str(it), ') error (deg) = ', num2str(mean_error_L1), ...
                ', med_L1 = ', num2str(median_error_L1), ...
                ', mean_L2 = ', num2str(mean_error_L2), ...
                ', med_L2 = ', num2str(median_error_L2)])
        end
        
    end
    time_optimization = toc;

    times_ROBA_per_iteration(dataset) = time_optimization/nIterations;
    disp(['Our optimization took ', num2str(time_optimization), 's.'])
    
    if (plotResult)
        figure;
        subplot(2,2,1)
        plot(ROBA_cost(dataset,:))
        ylim([0 inf])
        xlim([0 nIterations])
        title('Cost')
        xlabel('Iteration')
        subplot(2,2,2)
        plot( ROBA_mean_L1(dataset, :))
        ylim([0 inf])
        xlim([0 nIterations])
        title('Average Rot Error')
        ylabel('deg')
        xlabel('Iteration')
        subplot(2,2,[3 4])
        b = bar(1:nViews, [rot_individual_errors_init', individual_errors']);
        xlim([0.5 nViews+0.5])
        title('Individual Rot Errors')
        ylabel('deg')
        xlabel('Camera ID')
        legend({'Initial', 'Final'})
        b(1).FaceColor = [1 0 0];
        b(2).FaceColor = [0 1 0];

        switch dataset
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
    end

end
    
    
    
tEnd = toc(tStart);
disp([newline, 'Total experiment time = ', num2str(tEnd/3600), 'hrs.'])
%save('results\real_results_main.mat')



%% Function definitions

function [total_cost, g] = ComputeCostGradient(precomputed_mat, edge_IDs, R_est, dev_angle, approximate_gradient, square_rooting)


    total_cost = 0;
    nEdges = size(edge_IDs,2);
    nViews = length(R_est);
    
    R_est_x = nan(3,3,nViews);
    R_est_y = nan(3,3,nViews);
    R_est_z = nan(3,3,nViews);
    
    for i = 1:nViews
        rotvec = LogMap(R_est{i});
        R_est_x(:,:,i) = ExpMap(rotvec+[dev_angle;0;0]);
        R_est_y(:,:,i) = ExpMap(rotvec+[0;dev_angle;0]);
        R_est_z(:,:,i) = ExpMap(rotvec+[0;0;dev_angle]);
    end

    g = zeros(3*nViews, 1);
    for i = 1:nEdges
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        
        cost = GetPairwiseCost(R_est{j},R_est{k}, precomputed_mat{i}, square_rooting);
        total_cost = total_cost + cost;
        
        delta_cost_jx = GetPairwiseCost(R_est_x(:,:,j), R_est{k}, precomputed_mat{i}, square_rooting);
        delta_cost_jx = delta_cost_jx - cost;
        delta_cost_jy = GetPairwiseCost(R_est_y(:,:,j), R_est{k}, precomputed_mat{i}, square_rooting);
        delta_cost_jy = delta_cost_jy - cost;
        delta_cost_jz = GetPairwiseCost(R_est_z(:,:,j), R_est{k}, precomputed_mat{i}, square_rooting);
        delta_cost_jz = delta_cost_jz - cost;
                
        g(3*j-2) = g(3*j-2) + delta_cost_jx;
        g(3*j-1) = g(3*j-1) + delta_cost_jy;
        g(3*j)   = g(3*j)   + delta_cost_jz;

        if (approximate_gradient)
            g(3*k-2) = g(3*k-2) - delta_cost_jx;
            g(3*k-1) = g(3*k-1) - delta_cost_jy;
            g(3*k)   = g(3*k)   - delta_cost_jz;
        else
            delta_cost_kx = GetPairwiseCost(R_est{j}, R_est_x(:,:,k), precomputed_mat{i}, square_rooting);
            delta_cost_kx = delta_cost_kx - cost;
            delta_cost_ky = GetPairwiseCost(R_est{j}, R_est_y(:,:,k), precomputed_mat{i}, square_rooting);
            delta_cost_ky = delta_cost_ky - cost;
            delta_cost_kz = GetPairwiseCost(R_est{j}, R_est_z(:,:,k), precomputed_mat{i}, square_rooting);
            delta_cost_kz = delta_cost_kz - cost;
            
            g(3*k-2) = g(3*k-2) + delta_cost_kx;
            g(3*k-1) = g(3*k-1) + delta_cost_ky;
            g(3*k)   = g(3*k)   + delta_cost_kz;
        end

    end   
    g = g/dev_angle;
end


function precomputed_mat = GetNeighborData(rays_ref, rays_neighbor)
    xxF = zeros(3,3);
    xyF = zeros(3,3);
    xzF = zeros(3,3);
    yyF = zeros(3,3);
    yzF = zeros(3,3);
    zzF = zeros(3,3);
    for i =1:size(rays_ref,2)
        ray1 = rays_ref(:,i);
        ray2 = rays_neighbor(:,i);
        F = ray2*ray2';
        fx = ray1(1); fy = ray1(2); fz = ray1(3);

        xxF = xxF + fx*fx*F;
        xyF = xyF + fx*fy*F;
        xzF = xzF + fx*fz*F;
        yyF = yyF + fy*fy*F;
        yzF = yzF + fy*fz*F;
        zzF = zzF + fz*fz*F;
    end

    precomputed_mat.xxF = xxF;
    precomputed_mat.xyF = xyF;
    precomputed_mat.xzF = xzF;
    precomputed_mat.yyF = yyF;
    precomputed_mat.yzF = yzF;
    precomputed_mat.zzF = zzF;
    precomputed_mat.nPoints = size(rays_ref,2);
    
end

function cost = GetPairwiseCost(R_ref, R_neighbor, precomputed_mat, square_rooting)
    xxF = precomputed_mat.xxF;
    xyF = precomputed_mat.xyF;
    xzF = precomputed_mat.xzF;
    yyF = precomputed_mat.yyF;
    yzF = precomputed_mat.yzF;
    zzF = precomputed_mat.zzF;  

    cost = computeEdgeCost(R_ref, R_neighbor, xxF, xyF, xzF, yyF, yzF, zzF);
    
%     R12_est = R_ref*R_neighbor';
%     r1 = R12_est(1,:); r2 = R12_est(2,:); r3 = R12_est(3,:); 
%     m11 = r3*yyF*r3'-2*r3*yzF*r2'+r2*zzF*r2';
%     m22 = r1*zzF*r1'-2*r1*xzF*r3'+r3*xxF*r3';
%     m33 = r2*xxF*r2'-2*r1*xyF*r2'+r1*yyF*r1';
%     m12 = r1*yzF*r3'-r1*zzF*r2'-r3*xyF*r3'+r3*xzF*r2';
%     m21 = m12;
%     m13 = r2*xyF*r3'-r2*xzF*r2'-r1*yyF*r3'+r1*yzF*r2';
%     m31 = m13;
%     m23 = r1*xzF*r2'-r1*yzF*r1'-r3*xxF*r2'+r3*xyF*r1';
%     m32 = m23;
%     M = [m11 m12 m13; m21 m22 m23; m31 m32 m33];
%     [~,D] = eig(M);
%     cost = sqrt(abs(D(1,1)))

    if (~square_rooting)
        cost = cost^2;
    end

end


function [R_geo1, errors, mean_error_L1, median_error_L1] = AlignRotationL1(R_true, R_est)

    
    nViews = size(R_true,2);
    
    errors = zeros(1,nViews);
    R_transform = cell(1, nViews);
    
    for i = 1:nViews
        R_transform{i} = R_est{i}'*R_true{i};
    end

    % L1 averaging
  
    vectors_total = zeros(9,nViews);
    for i = 1:nViews
        vectors_total(:,i)= R_transform{i}(:);
    end
    med_vectors_total = median(vectors_total,2);
    [U,~,V] = svd(reshape(med_vectors_total, [3 3]));
    R_med = U*V.';
    if (det(R_med) < 0)
        V(:,3) = -V(:,3);
        R_med = U*V.';
    end

    R_geo1 = R_med;
    for j = 1:10
        step_num = 0;
        step_den = 0;
        for i = 1:nViews
            v =  LogMap(R_transform{i}*R_geo1');
            v_norm = norm(v);
            step_num = step_num + v/v_norm;
            step_den = step_den + 1/v_norm;
        end
        delta = step_num/step_den;
        delta_angle = norm(delta);

        delta_axis = delta/delta_angle;
        so3_delta = SkewSymmetricMatrix(delta_axis);
        R_delta = eye(3)+so3_delta*sin(delta_angle)+so3_delta^2*(1-cos(delta_angle));
        R_geo1 = R_delta*R_geo1;
        if (delta_angle < 0.001)
            break;
        end
    end
    
    

    for i = 1:nViews
        error = abs(acosd((trace(R_true{i}*(R_est{i}*R_geo1)')-1)/2));
        errors(i) = error;
    end
    mean_error_L1 = mean(errors);
    median_error_L1 = median(errors);
end

function [R_geo2, mean_error, median_error] = AlignRotationL2(R_true, R_est)
    
    nViews = size(R_true,2);
    errors = zeros(1,nViews);
    R_transform = cell(1, nViews);
    for i = 1:nViews
         R_transform{i} = R_est{i}'*R_true{i};
    end

    R_sum = zeros(3,3);
    for i = 1:nViews
        R_sum = R_sum + R_transform{i};
    end
    R_geo2 = ProjectOntoSO3(R_sum);
  
    for j = 1:10
        v = zeros(3,1);
        for i = 1:nViews
            v =  v + LogMap(R_transform{i}*R_geo2');
        end
        v = v/nViews;
        
        delta_angle = norm(v);
        
        R_delta = ExpMap(v);
        R_geo2 = R_delta*R_geo2;
        if (delta_angle < 0.001)
            break;
        end
    end
    
    for i = 1:nViews
        error = abs(acosd((trace(R_true{i}*(R_est{i}*R_geo2)')-1)/2));
        errors(i) = error;
    end
    
    mean_error = mean(errors);
    median_error = median(errors);
end


function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
end

function out = ExpMap(in)
    angle = norm(in);
    if (angle == 0)
        out = eye(3);
        return;
    end
    axis = in/angle;
    so3 = SkewSymmetricMatrix(axis);
    R = eye(3)+so3*sin(angle)+so3^2*(1-cos(angle));
    out = R;
end

function out = LogMap(in)
    if (in(1,1) == 1 && in(2,2) == 1 && in(3,3) == 1)
        out = [0;0;0];
        return;
    end
    
    cos_theta = min(1, max(-1, (trace(in)-1)/2));
    sin_theta = sqrt(1-cos_theta^2);
    
    if (sin_theta == 0)
        out = [0;0;0];
        return;
    end
    
    theta = acos(cos_theta);
    ln_R = theta/(2*sin_theta)*(in-in');
    out = [ln_R(3,2);ln_R(1,3);ln_R(2,1)];
end

function R = ProjectOntoSO3(M)   
    [U,~,V] = svd(M);
    R = U*V.';
    if (det(R) < 0)
        V(:,3) = -V(:,3);
        R = U*V.';
    end
end
