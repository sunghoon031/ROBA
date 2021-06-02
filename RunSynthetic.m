clear all; close all; clc;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Choose configuration:

draw3D = false;
config = 1;
% 1: Baseline
% 2: More points
% 3: Fewer views
% 4: More views
% 5: Closer points
% 6: Farther points
% 7: Less noise
% 8: More noise
% 9: Planar scene
% 10: Pure rotations
% 11: Pure + Planar
% 12: Mixed rotations
% 13: without square rooting
% 14: without switching alpha
% 15: without gradient approximation

%% Set configuration parameters:

max_angle_deg = 20; 
max_angle_est = 20;
focal_length = 525;
img_width = 640;
img_height = 480;
K = [focal_length 0 0; 0 focal_length 0; 0 0 1];
    
n_views = 100;
n_samples = 100;
n_min_cov = 50;
img_noise = 1;
d_min = 2;
d_max = 5;
pure_rotation = false;
mixed_rotation = false;
square_rooting = true;
switch_alpha = true;
approximate_gradient = true;
    
switch config
    case 2 % more points
        n_min_cov = 100;
    case 3 % fewer cameras
        n_views = 30;
    case 4 % more cameras
        n_views = 300;
    case 5 % closer points
        d_max = 3;
    case 6 % farther points
        d_max = 10;
    case 7 % less noisy points
        img_noise = 0.5;
    case 8 % more noisy points
        img_noise = 2;
    case 9 % planar scene
        d_min = 5;
    case 10 % pure rotation
        pure_rotation = true;
    case 11 % pure rotation + planar scene
        pure_rotation = true;
        d_min = 5;
    case 12 % mixed rotation (5 cams at 20 positions)
        mixed_rotation = true;
    case 13 % without square rooting
        square_rooting = false;
    case 14 % without switching alpha
        switch_alpha =false;
    case 15 % without gradient approximation
        approximate_gradient = false;
end

%% Simulate camera poses:

% Evenly distribute the cameras on a circle.
% Randomly perturb their orientations (within some bound).

theta = 2*pi/n_views;
if (mixed_rotation)
    theta = 2*pi/20;
end
    
radius = 1/norm([cos(2*theta)-cos(theta);sin(2*theta)-sin(theta)]);
scene_dim = 2*(radius+d_max*tand(60));
    
pos = cell(1, n_views);
for i = 1:n_views
    if (pure_rotation)
        pos{i} = [0;0;0];
    elseif (mixed_rotation)
        j = ceil(i/5);
        pos{i} = [radius*cos(j*theta);radius*sin(j*theta);0];
    else
        pos{i} = [radius*cos(i*theta);radius*sin(i*theta);0];
    end
end

R = cell(1, n_views);
t = cell(1,n_views);
for i = 1:n_views
    R{i} = RandomRotation(rand(1)*max_angle_deg);
    t{i} = -R{i}*pos{i};
end

%% Generate points:

points = [];
b_view_ready = zeros(1,n_views);
n_covisible_points_in_view = zeros(1,n_views);
while(sum(b_view_ready)<n_views)
    p = [scene_dim*(rand(1)-0.5); scene_dim*(rand(1)-0.5); d_min+rand(1)*(d_max-d_min)];

    % Add the point p if it is visible in two neighboring views that do not
    % have enough covisible points yet. 
    
    good_point = false;
    for i = 1:n_views
        if (b_view_ready(i))
            continue;
        end
        pc = R{i}*p+t{i};
        uv = K*[pc(1)/pc(3); pc(2)/pc(3);1];
        if (abs(uv(1)) > img_width/2 || abs(uv(2)) > img_height/2 || pc(3) < 0)
            continue;
        end

        % Make sure two neighboring views share enough points.
        if (i == n_views)
            j = 1;
        else
            j = i+1;
        end

        pc = R{j}*p+t{j};
        uv = K*[pc(1)/pc(3); pc(2)/pc(3);1];
        if (abs(uv(1)) > img_width/2 || abs(uv(2)) > img_height/2 || pc(3) < 0)
            continue;
        end

        good_point = true;

        n_covisible_points_in_view(i) = n_covisible_points_in_view(i) + 1;
        if (n_covisible_points_in_view(i) >= n_min_cov)
            b_view_ready(i) = 1;
            %disp(['Good cams = ', num2str(sum(b_good))])
        end
    end

    if (good_point)
        points = [points, p];
    end
end

n_points = size(points,2);

%% Add noise to image and obtain the rays:

points_in_each_view = cell(1, n_views);
rays_in_each_view = cell(1, n_views);
for i = 1:n_points
    p = points(:,i);
    for j = 1:n_views
        pc = R{j}*p+t{j};
        uv = K*[pc(1)/pc(3); pc(2)/pc(3);1];
        if (abs(uv(1)) > img_width/2 || abs(uv(2)) > img_height/2 || pc(3) < 0)
            continue;
        end
        points_in_each_view{j} = [points_in_each_view{j}, i];

        noise_mag = abs(normrnd(0, img_noise));
        while(true)
            noise_dir = rand(2,1)-0.5;
            noise_dir = noise_dir/norm(noise_dir);
            uv_pert = uv(1:2) + noise_mag*noise_dir;
            if (abs(uv_pert(1)) < img_width/2 && abs(uv_pert(2)) < img_height/2)
                break;
            end
        end
        ray = K\[uv_pert;1];
        ray = ray/norm(ray);
        rays_in_each_view{j} = [rays_in_each_view{j}, ray];
    end
end

%% Get edges:
edges = [];
edge_idx = 0;
for i = 1:n_views
    for j = 1:n_views
        if (i<=j)
            continue;
        end

        intersection = intersect(points_in_each_view{i},points_in_each_view{j});

        if (length(intersection) >= 10)
            edges = [edges, [i;j]];
            edge_idx = edge_idx + 1;
            points_for_each_edge{edge_idx} = intersection;
        end
    end
end

n_edges = size(edges,2);
n_inliers = [];
rot_errors = [];
edges2 = [];
c = 0;

%% Remove weak edges and get relative rotations:
for i = 1:n_edges
    %disp([num2str(i/n_edges*100), '%'])

    j = edges(1,i);
    k = edges(2,i);

    [~, loc_j] = ismember(points_for_each_edge{i},points_in_each_view{j});
    [~, loc_k] = ismember(points_for_each_edge{i},points_in_each_view{k});

    rays_j = rays_in_each_view{j}(:,loc_j);
    rays_k = rays_in_each_view{k}(:,loc_k);

    thr = 0.02;
    [R_jk, t_out, inlier_idx] = RotationFromGT(R{j}*R{k}',  -R{j}*R{k}'*t{k}+t{j}, rays_j, rays_k, n_samples, thr, max_angle_est);

    if (length(inlier_idx) >= 10)
        c = c + 1;
        n_inliers = [n_inliers, length(inlier_idx)];
        rot_errors = [rot_errors, DegBtnRotations(R_jk, R{j}*R{k}')];
        edges2 = [edges2, edges(:,i)];
        precomputed_mat{c} = GetNeighborData(rays_j(:, inlier_idx), rays_k(:, inlier_idx));
        RR(:,:,c) = R_jk;
    end
end
mean_rot_error = mean(rot_errors);
median_rot_error = median(rot_errors);

edges = edges2;
n_edges = c;
edge_percent = n_edges/(n_views*(n_views-1)/2)*100;

edges_reverse = edges;
edges_reverse(1,:) = edges(2,:);
edges_reverse(2,:) = edges(1,:);

disp('Simulation ready!')

%% Rotation averaging:
R_avg_mat = AverageSO3Graph(RR,edges_reverse, 'Method', 'L0.5');
R_avg = cell(1,n_views);
for i = 1:n_views
    R_avg{i} = R_avg_mat(:,:,i);
end
[~,~, mn1_RA, md1_RA] = AlignRotationL1(R, R_avg);
[~, mn2_RA, md2_RA] = AlignRotationL2(R, R_avg);

disp('RA complete!')

%% ROBA:
dev_angle = 10^(-4);
n_iterations = 100;

R_est = R_avg;
r = nan(3*n_views,1);
for i = 1:n_views
    r(3*i-2:3*i) = LogMap(R_est{i});
end

cost_history = nan(1,n_iterations+1);
total_cost = 0;
for i = 1:size(edges, 2) 
    j = edges(1,i);
    k = edges(2,i);
    total_cost = total_cost + GetPairwiseCost(R_est{j}, R_est{k}, precomputed_mat{i}, square_rooting);
end
cost_history(1) = total_cost;

%ADAM parameters:
alpha = 0.01; beta1 = 0.9; beta2 =0.999; epsilon = 10^(-8);   
m_prev = 0;
v_prev = 0;
r_prev = zeros(1, 3*n_views);

costUpCount = 0;

for it = 1:n_iterations  

    total_cost_prev = total_cost;
    [total_cost, g] = ComputeCostGradient(precomputed_mat, edges, R_est, dev_angle, approximate_gradient, square_rooting);
    cost_history(it+1) = total_cost;

    % ADAM:
    m = beta1*m_prev + (1-beta1)*g;
    v = beta2*v_prev + (1-beta2)*(g.^2);
    m_prev = m;
    v_prev = v;
    m_hat = m/(1-beta1^it);
    v_hat = v/(1-beta2^it);
    r = r - alpha*m_hat./(sqrt(v_hat)+epsilon);

    for i = 1:n_views
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

end

disp('ROBA complete!')

[~,~, mn1_ROBA, md1_ROBA] = AlignRotationL1(R, R_est);
[~, mn2_ROBA, md2_ROBA] = AlignRotationL2(R, R_est);



%% Evaluate results

disp(['config: ', num2str(config), ...
    ', edge: ', num2str(edge_percent), ...
    '%, mn edge err: ', num2str(mean_rot_error), ...
    'deg, md edge err: ', num2str(median_rot_error), ...
    'deg, RA: ', num2str(mn1_RA), ...
    'deg, ROBA: ', ...
    num2str(mn1_ROBA), 'deg'])

figure;
plot(cost_history)
ylabel('Total Cost')
xlabel('Number of iterations')

if (~draw3D)
    return;
end

figure;
hold on
axis equal
for i = 1:n_points
    scatter3(points(1,i), points(2,i),points(3,i), 'k.')
end
for i = 1:n_views
    scatter3(pos{i}(1),pos{i}(2),pos{i}(3), 'bo')
    
    ep = R{i}'*([0;0;5] - t{i});
    sp = pos{i};
    plot3([sp(1),ep(1)],[sp(2),ep(2)],[sp(3),ep(3)], 'r-')
end
for i = 1:n_edges
    j = edges(1,i);
    k = edges(2,i);
    
    ep = pos{j};
    sp = pos{k};
    plot3([sp(1),ep(1)],[sp(2),ep(2)],[sp(3),ep(3)], 'g-')
end




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

function out = RandomRotation(angle_deg)
    axis = rand(3,1)-0.5;
    axis = axis/norm(axis);
    angle = angle_deg/180*pi;
    rotvec = angle*axis;
    out = ExpMap(rotvec);
end

function out = RandomRotation2(unit_axis, angle_deg)
    angle = angle_deg/180*pi;
    rotvec = angle*unit_axis;
    out = ExpMap(rotvec);
end

function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
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

function theta = DegBtnRotations(R1, R2)
    if (isnan(R1(1,1)) || isnan(R2(1,1)))
        theta = nan;
        return;
    end
    cos_theta = (trace(R1*R2')-1)/2;
    cos_theta = max(-1, min(1, cos_theta));
    theta = acos(cos_theta);
    theta = theta/pi*180;
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

function [R_out, t_out, inlier_idx] = RotationFromGT(R_gt, t_gt, rays1, rays2, nSamples, thr, max_angle_est)
    % Obtain pose samples around the ground-truth relative pose.
    % Then, choose the one that yileds the most inliers based on L1-optimal
    % angular reprojection errors.
    
    inlier_idx = [];
    R_out = nan(3,3);
    t_out = nan(3,1);
    
    if (norm(t_gt)==0) % pure rotation
        t_gt = rand(3,1)-0.5;
        t_gt = t_gt/norm(t_gt);
    end
    
    nPoints = size(rays1,2);

    maxInliers = 0;
    for it = 1:nSamples
        R_temp = R_gt*RandomRotation(rand(1)*max_angle_est);
        axis = cross(t_gt, rand(3,1)-0.5);
        axis = axis/norm(axis);
        t_temp = RandomRotation2(axis, rand(1)*max_angle_est)*t_gt;
        E_temp = SkewSymmetricMatrix(t_temp)*R_temp;
        
        % Count consensus
        n_inliers = 0;
        inliers = [];
        for i = 1:nPoints
            f1 = rays1(:,i);
            f2 = rays2(:,i);
            nee = abs(f1'*E_temp*f2);
            Rf2 = R_temp*f2;
            phi0 = AngleRadBtn(Rf2, t_temp);
            phi1 = AngleRadBtn(f1, t_temp);
            sin_max = min(1, sin(max(phi0,phi1)));
            sin_L1 = nee/sin_max;
            L1 = abs(asin(sin_L1));
            
            if (L1 < thr)
                n_inliers = n_inliers + 1;
                inliers = [inliers, i];
            end
        end

        if (n_inliers > maxInliers)
            maxInliers = n_inliers;
            R_out = R_temp;
            t_out = t_temp;
            inlier_idx = inliers;
        end
        
    end
end

function a = AngleRadBtn(b,c)
    if (b(1)==c(1) && b(2) == c(2) && b(3) == c(3))
        a = 0;
    elseif (b(1)==-c(1) && b(2) == -c(2) && b(3) == -c(3))
        a = 0;
    else
        a = atan2(norm(cross(b,c)),dot(b,c));
        a = min(a, pi-a);
    end
end

function R = ProjectOntoSO3(M)   
    [U,~,V] = svd(M);
    R = U*V.';
    if (det(R) < 0)
        V(:,3) = -V(:,3);
        R = U*V.';
    end
end