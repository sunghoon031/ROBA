clear all; close all; clc;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Load dataset

dataset = 2; 
% 1: Alamo, 2: Ellis_Island, 3: Gendarmenmarkt
% 4: Madrid_Metropolis, 5: Montreal_Notre_Dame, 6: Notre_Dame
% 7: NYC_Library, 8: Piazza_del_Popolo, 9: Piccadilly
% 10: Roman_Forum, 11: Tower_of_London, 12: Trafalgar
% 13: Union_Square, 14: Vienna_Cathedral, 15: Yorkminster



nMinCovisiblePoints = 10;
draw2verify = true;

tic;
[views, points, tracks, edge_info, edge_IDs, edge_IDs_reverse, RR, edge_errors, nViews, nPoints, nEdges, nEdgesForEachView, nPointsForEachEdge, nPointsForEachView] =...
    LoadDataset(dataset, draw2verify, nMinCovisiblePoints);
toc

disp(['Dataset loaded. nViews = ', num2str(nViews), ', nEdges = ', num2str(nEdges)])


%%
function [views, points, tracks, edge_info,  edge_IDs, edge_IDs_reverse, RR, edge_errors, nViews, nPoints, nEdges, nEdgesForEachView, nPointsForEachEdge, nPointsForEachView] ...
    = LoadDataset(dataset, draw2verify, nMinCovisiblePoints)
    
    switch dataset
        case 1
            filename1 = 'Alamo.out';
            filename2 = 'Alamo_EGs.txt';
            filename3 = 'Alamo_coords.txt';
            filename4 = 'Alamo_tracks.txt';
        case 2
            filename1 = 'Ellis_Island.out';
            filename2 = 'Ellis_Island_EGs.txt';
            filename3 = 'Ellis_Island_coords.txt';
            filename4 = 'Ellis_Island_tracks.txt';
        case 3
            filename1 = 'Gendarmenmarkt.out';
            filename2 = 'Gendarmenmarkt_EGs.txt';
            filename3 = 'Gendarmenmarkt_coords.txt';
            filename4 = 'Gendarmenmarkt_tracks.txt';
        case 4
            filename1 = 'Madrid_Metropolis.out';
            filename2 = 'Madrid_Metropolis_EGs.txt';
            filename3 = 'Madrid_Metropolis_coords.txt';
            filename4 = 'Madrid_Metropolis_tracks.txt';
        case 5
            filename1 = 'Montreal_Notre_Dame.out';
            filename2 = 'Montreal_Notre_Dame_EGs.txt';
            filename3 = 'Montreal_Notre_Dame_coords.txt';
            filename4 = 'Montreal_Notre_Dame_tracks.txt';
        case 6
            filename1 = 'Notre_Dame.out';
            filename2 = 'Notre_Dame_EGs.txt';
            filename3 = 'Notre_Dame_coords.txt';
            filename4 = 'Notre_Dame_tracks.txt';
        case 7
            filename1 = 'NYC_Library.out';
            filename2 = 'NYC_Library_EGs.txt';
            filename3 = 'NYC_Library_coords.txt';
            filename4 = 'NYC_Library_tracks.txt';
        case 8
            filename1 = 'Piazza_del_Popolo.out';
            filename2 = 'Piazza_del_Popolo_EGs.txt';
            filename3 = 'Piazza_del_Popolo_coords.txt';
            filename4 = 'Piazza_del_Popolo_tracks.txt';
        case 9
            filename1 = 'Piccadilly.out';
            filename2 = 'Piccadilly_EGs.txt';
            filename3 = 'Piccadilly_coords.txt';
            filename4 = 'Piccadilly_tracks.txt';
        case 10
            filename1 = 'Roman_Forum.out';
            filename2 = 'Roman_Forum_EGs.txt';
            filename3 = 'Roman_Forum_coords.txt';
            filename4 = 'Roman_Forum_tracks.txt';
        case 11
            filename1 = 'Tower_of_London.out';
            filename2 = 'Tower_of_London_EGs.txt';
            filename3 = 'Tower_of_London_coords.txt';
            filename4 = 'Tower_of_London_tracks.txt';
        case 12
            filename1 = 'Trafalgar.out';
            filename2 = 'Trafalgar_EGs.txt';
            filename3 = 'Trafalgar_coords.txt';
            filename4 = 'Trafalgar_tracks.txt';
        case 13
            filename1 = 'Union_Square.out';
            filename2 = 'Union_Square_EGs.txt';
            filename3 = 'Union_Square_coords.txt';
            filename4 = 'Union_Square_tracks.txt';
        case 14
            filename1 = 'Vienna_Cathedral.out';
            filename2 = 'Vienna_Cathedral_EGs.txt';
            filename3 = 'Vienna_Cathedral_coords.txt';
            filename4 = 'Vienna_Cathedral_tracks.txt';
        case 15
            filename1 = 'Yorkminster.out';
            filename2 = 'Yorkminster_EGs.txt';
            filename3 = 'Yorkminster_coords.txt';
            filename4 = 'Yorkminster_tracks.txt';
    end

    data_ = fileread(filename1);
    data_ = cellstr(data_);
    data_ = cellfun(@(newline) strsplit(newline, '\n'), data_, 'UniformOutput', false);
    data_ = [data_{:}];
    data_ = data_';        
    [data, ~] = sscanf(data_{2}, '%f');
    nViews = data(1);
    nPoints = data(2);

    % Find the number of views with groundtruth
    nViewsGt = 0;
    for i = 3:(nViews+1)*5-3
        if (mod(i-3,5)>0)
            continue;
        end
        [data, ~] = sscanf(data_{i}, '%f');
        if(data(1)==0.0 && data(2)==0.0 && data(3)==0.0)
            continue;
        end
        nViewsGt = nViewsGt + 1;
    end

    views.ID = nan(1, nViewsGt);
    views.ID_old = nan(1, nViewsGt);
    views.width = nan(1, nViewsGt);
    views.height = nan(1, nViewsGt);
    views.R = cell(1, nViewsGt);
    views.t = cell(1, nViewsGt);
    views.K = cell(1, nViewsGt);
    views.radialDistortion = cell(1, nViewsGt);
    views.pos = cell(1, nViewsGt);
    views.points_ID = cell(1, nViewsGt);
    views.points_uv = cell(1, nViewsGt);
    views.SIFT_ID = cell(1, nViewsGt);
    views.SIFT_uv = cell(1, nViewsGt);
    views.connected = zeros(1, nViewsGt);

    %% Parse camera data
    cam_idx = 0;
    cam_idx_old = 0;
    row_count = -1;
    for i = 3:(nViews+1)*5-3
        [data, ~] = sscanf(data_{i}, '%f');
        row_count = row_count + 1;

        if (mod(i-2,5) == 0)
            cam_idx_old = cam_idx_old +1;
        end

        if(data(1)==0.0 && data(2)==0.0 && data(3)==0.0)
            continue;
        end


        switch mod(row_count,5)
            case 0
                cam_idx = cam_idx + 1;
                views.ID(cam_idx) = cam_idx;
                views.ID_old(cam_idx) = cam_idx_old;
                
                views.K{cam_idx} = [data(1), 0 0; 0 data(1) 0; 0 0 1]; 
                views.radialDistortion{cam_idx} = [data(2), data(3)];

            case 1
                views.R{cam_idx} = eye(3);
                views.R{cam_idx}(1,:) = data';
            case 2
                views.R{cam_idx}(2,:) = data';
            case 3
                views.R{cam_idx}(3,:) = data';
            case 4
                views.t{cam_idx} = data;
        end
    end

    i_prev = i;

    % Change the camera reference such that z-axis is forward!
    c_all = nan(3,nViewsGt);
    for i = 1:nViewsGt
        views.R{i} = diag([-1 1 -1])*views.R{i};
        views.t{i} = diag([-1 1 -1])*views.t{i};
        views.pos{i} = -views.R{i}'*views.t{i};
        c_all(:,i) = views.pos{i};
    end
    
    
    %% Parse point data
    points.ID = nan(1, nPoints);
    points.pos = cell(1, nPoints);
    points.views = cell(1, nPoints);
    
    p_all = nan(3,nPoints);
    color_all = nan(3, nPoints);

    point_idx = 0;
    row_count = 0;
    for i = i_prev+1:length(data_)
        [data, ~] = sscanf(data_{i}, '%f');
        row_count = row_count + 1;
        views_temp = [];
        switch mod(row_count,3)
            case 1
                point_idx = point_idx + 1;
                point = [data(1); data(2); data(3)];
                p_all(:,point_idx) = point;
            case 2
                color_all(:,point_idx) = data/256;
            case 0
                data = data(2:end);
                for j = 1:length(data)
                    switch mod(j, 4)
                        case 1
                            cam_idx = data(j)+1;
                            views_temp = [views_temp, cam_idx];
                            
                            point_c = views.R{cam_idx}*point+views.t{cam_idx};
                            uv = views.K{cam_idx}*point_c/(point_c(3));
                            uv = uv(1:2);
                            
                            views.points_ID{cam_idx} = [views.points_ID{cam_idx}, point_idx];
                            views.points_uv{cam_idx} = [views.points_uv{cam_idx}, uv];     
                    end
                end
        end

        points.ID(point_idx) = point_idx;
        points.pos{point_idx} = point;
        points.views{point_idx} = views_temp;
    end
    
    %% Get image size and SIFT
    data_ = fileread(filename3);
    data_ = cellstr(data_);
    data_ = cellfun(@(newline) strsplit(newline, '\n'), data_, 'UniformOutput', false);
    data_ = [data_{:}];
    data_ = data_';      

    header_idx = 1;
    uvs_SIFT = [];
    for i = 1:length(data_)
        if (i == header_idx)           
            C = textscan(data_{i},'%s %s %s %s %s %s', 'Delimiter',',');
            ID_old = sscanf(string(C{1}), '#index = %f');
            nKeys =  sscanf(string(C{3}), 'keys = %f');
            width = 2*sscanf(string(C{4}), 'px = %f');
            height = 2*sscanf(string(C{5}), 'py = %f');

            cam_idx = views.ID(find(views.ID_old == ID_old));
            if (~isempty(cam_idx))
                views.width(cam_idx) = width;
                views.height(cam_idx) = height;
            end

            header_idx = header_idx + nKeys +1;
        elseif (~isempty(cam_idx))
            [data, ~] = sscanf(data_{i}, '%f');
            views.SIFT_ID{cam_idx} = [views.SIFT_ID{cam_idx}, data(1)];
            uv_SIFT = [data(2)-width/2, data(3)-height/2];
            uvs_SIFT = [uvs_SIFT; uv_SIFT];     
            
            if (i == header_idx-1)
                disp(['Undistortion: ', num2str(i/length(data_)*100), '%, number of SIFT points in current img = ', num2str(length(uvs_SIFT))])
                
                cameraParams = cameraParameters('IntrinsicMatrix',views.K{cam_idx},'RadialDistortion',views.radialDistortion{cam_idx}); 

                % Partition the points for undistortion to prevent 
                % excessive memory use. We divide them into groups of 1000.
                k_prev = 0;
                k = min(1000, size(uvs_SIFT,1));
                for j = 1:size(uvs_SIFT,1)
                    if (j==k)
                        %disp(['k = ', num2str(k)])
                        views.SIFT_uv{cam_idx} = [views.SIFT_uv{cam_idx} , (-undistortPoints(uvs_SIFT(k_prev+1:k,:),cameraParams))'];
                        k_prev = k;
                        k = min(k + 1000, size(uvs_SIFT,1));
                    end
                end
                uvs_SIFT = [];
            end
        end
    end


        

    %% Get edges
    data_ = fileread(filename2);
    data_ = cellstr(data_);
    data_ = cellfun(@(newline) strsplit(newline, '\n'), data_, 'UniformOutput', false);
    data_ = [data_{:}];
    data_ = data_';      

    edge_IDs = [];
    edge_rs = [];
    edge_info = cell(1, nViewsGt);
    connectedViewIDs = [];
    for i = 1:length(data_)
        [data, ~] = sscanf(data_{i}, '%f');
        j = data(1);
        k = data(2);
        r = data(3:11);
        R = reshape(r, [3,3])';
        R = diag([-1 1 -1])*R*diag([-1 1 -1]);
        r = R(:);
        
        jj = find(views.ID_old==j);
        kk = find(views.ID_old==k);
        if (isempty(jj) || isempty(kk))
            continue;
        end

        % Skip the edges with insufficient covisible points.
        % They exist probabily because Bundler didn't reconstruct some of the points.
        intersection = intersect(views.points_ID{jj}, views.points_ID{kk});
        if (length(intersection) < nMinCovisiblePoints)
            continue;
        end
        
        if (jj > kk)
            jj_ = jj;
            jj = kk;
            kk = jj_;
        end

        edge_IDs = [edge_IDs, [jj;kk]];
        edge_rs = [edge_rs, r];
        edge_info{jj} = [edge_info{jj}, kk];
        edge_info{kk} = [edge_info{kk}, jj];
        
        if (~ismember(jj, connectedViewIDs))
            connectedViewIDs = [connectedViewIDs, jj];
        end

        if (~ismember(kk, connectedViewIDs))
            connectedViewIDs = [connectedViewIDs, kk];
        end
    end
    
    %% Remove views that are not connected to the main spanning tree
    nTotalConntected = length(connectedViewIDs);
    for i = 1:length(connectedViewIDs)
        id = connectedViewIDs(i); % Root node.
        connections = id;
        
        % Incrementally add nodes stemming from the root node:
        % (1) Add all the nighbors of the root node to 'not_checked_yet'.
        % (2) Add the first node in 'not_checked_yet' 
        %     (i.e., 'first_not_checked') to 'connections'.
        % (3) Add all the neighbors of 'first_not_checked' in
        %     'not_checked_yet', unless it already exists in 'connections'
        % (4) Remove 'first_not_checked' from 'not_checked_yet'.
        % (5) Repeat Step 2-4 until 'not_checked_yet' is empty.
        % (6) If more than half of the total nodes are in 'connections', it
        %     means that we found the main spanning tree, so terminate.
        %     Otherwise, choose the next node as the root node, and repeat.
        
        not_checked_yet = edge_info{id}; % Step 1
        while (~isempty(not_checked_yet))
            first_not_checked = not_checked_yet(1);
            connections = [connections, first_not_checked]; % Step 2
            for j = edge_info{first_not_checked}
                if (~ismember(j, connections) && ~ismember(j, not_checked_yet))
                    not_checked_yet = [not_checked_yet, j]; % Step 3
                end
            end
            not_checked_yet(1) = []; % Step 4
        end
        if (length(connections) > nTotalConntected*0.5)
            break;
        end
    end
    
    
    connectedViewIDs  = connections;
    
    edge_IDs_ = [];
    edge_rs_ = [];
    for i = 1:length(edge_IDs)
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        r = edge_rs(:,i);
        
        if (ismember(j, connectedViewIDs) && ismember(k, connectedViewIDs))
            edge_IDs_ = [edge_IDs_, [j;k]];
            edge_rs_ = [edge_rs_, r];
        end
    end
    edge_IDs = edge_IDs_;
    edge_rs = edge_rs_;
    
    for i = 1:nViewsGt
        if (~ismember(i, connectedViewIDs))
            edge_info{i} = [];
            continue;
        end
        
        edge_info{i} = intersect(edge_info{i}, connectedViewIDs);
    end
    
    
    
    %% Sort
    connectedViewIDs = sort(connectedViewIDs);
    edge_IDs_ = [];
    edge_rs_ = [];
    for i = 1:nViewsGt
        if (isempty(edge_info{i}))
            continue;
        end
        edge_info{i} = sort(edge_info{i});
        
        for j = 1:size(edge_IDs,2)
            if (edge_IDs(1,j)==i)
                edge_IDs_ = [edge_IDs_, edge_IDs(:,j)];
                edge_rs_ = [edge_rs_, edge_rs(:,j)];
            end
        end
    end
    edge_IDs = edge_IDs_;
    edge_rs = edge_rs_;
    
    %% Reindex views
    
    oldViewGT_to_newViewGT = nan(1, nViewsGt);
    
    for i = 1:nViewsGt
        new_id = find(connectedViewIDs==i);
        if (~isempty(new_id))
            oldViewGT_to_newViewGT(i) = new_id;
        end
    end
    
    j = 0;
    for i = 1:nViewsGt
        if (isnan(oldViewGT_to_newViewGT(i)))
            continue;
        end
        
        j = j + 1;
        views_.ID(j) = j;
        views_.ID_old(j) = views.ID_old(i);
        views_.width(j) = views.width(i);
        views_.height(j) = views.height(i);
        views_.R{j} = views.R{i};
        views_.t{j} = views.t{i};
        views_.K{j} = views.K{i};
        views_.radialDistortion{j} = views.radialDistortion{i};
        views_.pos{j} = views.pos{i};
        views_.points_ID{j} = views.points_ID{i};
        views_.points_uv{j} = views.points_uv{i};
        views_.SIFT_ID{j} = views.SIFT_ID{i};
        views_.SIFT_uv{j} = views.SIFT_uv{i};
        
        
        edge_info_{j} = oldViewGT_to_newViewGT(edge_info{i});
    end

    edge_IDs_ = nan(size(edge_IDs));
    edge_rs_ = nan(size(edge_rs));
    j = 0;
    for i = 1:size(edge_IDs,2)
        edge_ID1 = oldViewGT_to_newViewGT(edge_IDs(1,i));
        edge_ID2 = oldViewGT_to_newViewGT(edge_IDs(2,i));
        if (isnan(edge_ID1) || isnan(edge_ID2))
            continue;
        end
        j = j + 1;
        edge_IDs_(:,j) = [edge_ID1;edge_ID2];
        edge_rs_(:,j) = edge_rs(:,i);
    end
   
    views = views_;
    edge_info = edge_info_;
    edge_IDs = edge_IDs_;
    edge_rs = edge_rs_;
    
    nEdgesForEachView = zeros(1, length(edge_info));
    for i = 1:length(edge_info)
        nEdgesForEachView(i) = length(edge_info{i});
    end
    
    for i = 1:nPoints
        points.views{i} = oldViewGT_to_newViewGT(points.views{i});
        points.views{i} = points.views{i}(~isnan(points.views{i}));
    end
    
    old_IDplus1_to_new_ID = nan(1,cam_idx_old);
    for i = 0:cam_idx_old-1
        for j = 1:length(views.ID_old)
            if (views.ID_old(j)==i)
                old_IDplus1_to_new_ID(i+1) = views.ID(j);
                break;
            end
        end
    end
    
    
    %% Reindex points
    oldPoint_to_newPoint = nan(1, nPoints);
    
    j = 0;
    for i = 1:nPoints 
        if (length(points.views{i})<2)
            continue;
        end
        
        j = j + 1;
        oldPoint_to_newPoint(i) = j;

        points_.ID(j) = j;
        points_.ID_old(j) = i;
        points_.pos{j} = points.pos{i};
        points_.views{j} = points.views{i};
    end
    points = points_;
    
    for i = 1:length(views.ID)
        views.points_ID{i} = oldPoint_to_newPoint(views.points_ID{i});
        views.points_uv{i} = views.points_uv{i}(:,~isnan(views.points_ID{i}));
        views.points_ID{i} = views.points_ID{i}(~isnan(views.points_ID{i}));
    end
   
    
    %% Make the first camera location the origin
    % I found that this is necessary in order for the Matlab's
    % bundleadjustment function to work... I think it's because there is a
    % problem when the points coordinates are too large.
    
    first_view_pos = views.pos{1};
    for i = 1:length(views.ID)
        views.pos{i} = views.pos{i} - first_view_pos;
        views.t{i} = -views.R{i}*views.pos{i};
    end
    
    for i = 1:length(points.ID)
        points.pos{i} = points.pos{i} - first_view_pos;
    end
    
    
    %% Get Tracks
    data_ = fileread(filename4);
    data_ = cellstr(data_);
    data_ = cellfun(@(newline) strsplit(newline, '\n'), data_, 'UniformOutput', false);
    data_ = [data_{:}];
    data_ = data_';      

    track_idx = 0;
    for i = 1:length(data_)
        [data, ~] = sscanf(data_{i}, '%f');

        if (i==1)
            tracks.views = cell(1, data(1));
            tracks.SIFT_IDs = cell(1, data(1));
            continue;
        end

        track_idx = track_idx + 1;

        for j = 2:length(data)
            if (mod(j,2)==0)
                view_ID = old_IDplus1_to_new_ID(data(j)+1);
                SIFT_ID = data(j+1);
                if (~isnan(view_ID))
                    tracks.views{track_idx} = [tracks.views{track_idx}, view_ID];
                    tracks.SIFT_IDs{track_idx} = [tracks.SIFT_IDs{track_idx}, SIFT_ID];
                end
            end
        end
    end

    j = 0;
    for i = 1:length(tracks.views)
        if (length(tracks.views{i})>=2)
            j = j + 1;
            tracks_.views{j} = tracks.views{i};
            tracks_.SIFT_IDs{j} = tracks.SIFT_IDs{i};
        end
    end
    tracks = tracks_;
    
    
     %% Associate points and tracks

    points.track_ID = nan(1,length(points.ID));
    points.mean_reproj_error = nan(1,length(points.ID));
    points.SIFT_IDs = cell(1, length(points.ID));

    tracks_length = nan(1, length(tracks.views));
    tracks_view_max = nan(1, length(tracks.views));
    tracks_view_min = nan(1, length(tracks.views));
    for j = 1:length(tracks.views)
        track_view_IDs = tracks.views{j};
        tracks_length(j) = length(track_view_IDs);
        tracks_view_max(j) = max(track_view_IDs);
        tracks_view_min(j) = min(track_view_IDs);
    end

    for i = 1:length(points.ID)
        view_IDs = points.views{i};
        nViews = length(view_IDs);

        uvs_reproj = nan(2,nViews);
        for k = 1:nViews
            x_w = points.pos{i};
            view_IDs = points.views{i};
            view_ID = view_IDs(k);
            x_c = views.R{view_ID}*x_w+views.t{view_ID};
            uv_reproj = views.K{view_ID}*(x_c/x_c(3));
            uv_reproj = uv_reproj(1:2);
            uvs_reproj(:,k) = uv_reproj;
        end
        points.uvs_reproj{i} = uvs_reproj;
    end

    tStart = tic;
    points_to_check_next = 1:length(points.ID);
    for thr = [5 inf]
        points_to_check = points_to_check_next;
        points_to_check_next = [];

        c = 0;
        for i = points_to_check
            c = c+1;
            view_IDs = points.views{i};
            nViews = length(view_IDs);
            uvs_reproj = points.uvs_reproj{i};

            minCost = inf;
            for j = 1:length(tracks.views)  
                if (tracks_length(j) < nViews)
                    continue;
                end
                if (tracks_view_max(j) < view_IDs(end))
                    continue;
                end
                if (tracks_view_min(j) > view_IDs(1))
                    continue;
                end

                track_view_IDs = tracks.views{j};
                track_SIFT_IDs = tracks.SIFT_IDs{j};
                [isMember, foundLocation] = ismember(view_IDs, track_view_IDs);
                isSubset = all(isMember);
                if (~isSubset)
                    continue;
                end

                cost = 0;
                SIFT_IDs = nan(1,nViews);
                for k = 1:nViews
                    view_ID = view_IDs(k);
                    SIFT_ID = track_SIFT_IDs(foundLocation(k));
                    uv_SIFT = views.SIFT_uv{view_ID}(:,SIFT_ID+1);

                    uv_diff = uv_SIFT-uvs_reproj(:,k);
                    uv_diff_sq = uv_diff.^2;
                    if (uv_diff_sq(1) > thr^2 || uv_diff_sq(2) > thr^2)
                        cost = inf;
                        break;
                    end

                    cost = cost + sqrt(uv_diff_sq(1)+uv_diff_sq(2));
                    SIFT_IDs(k) = SIFT_ID;

                    if (cost > minCost)
                        break;
                    end

                end

                if (cost < minCost)
                    minCost = cost;
                    points.track_ID(i) = j;
                    points.mean_reproj_error(i) = minCost/length(view_IDs);
                    points.SIFT_IDs{i} = SIFT_IDs;

                    if (points.mean_reproj_error(i) < 3)
                        break;
                    end
                end
            end    

            if (isnan(points.mean_reproj_error(i)))
                points_to_check_next = [points_to_check_next, i];
                %disp(['nan! points to check without limit = ', num2str(length(points_to_check_without_limit))])
            end

            
            if (rand(1) < 0.1)
                tNow = toc(tStart);
                estimated_time_remaining = tNow/c*(length(points_to_check)-c)/3600;
                disp([num2str(estimated_time_remaining), ' hrs to go, thr = ', num2str(thr), ', ', num2str(c/length(points_to_check)*100), '%, track ID = ', num2str(points.track_ID(i)), ', mean reproj error = ', num2str(points.mean_reproj_error(i)), ' pix'])
            end
        end
    end

    
     
    %% Obtain image masurements in each view

    for i = 1:length(views.ID)
        views.points_SIFT_uv{i} = nan(2, length(views.points_ID{i}));
        views.rays_SIFT{i} = nan(3, length(views.points_ID{i}));

        for j = 1:length(views.points_ID{i})
            point_ID = views.points_ID{i}(j);
            foundLocation = points.views{point_ID}==i;
            SIFT_ID = points.SIFT_IDs{point_ID}(foundLocation);
            uv_SIFT = views.SIFT_uv{i}(:,SIFT_ID+1);
            views.points_SIFT_uv{i}(:,j) = uv_SIFT;
            f = views.K{i}\[uv_SIFT; 1];
            f = f/norm(f);
            views.rays_SIFT{i}(:,j) = f;
        end
    end

    %% Parse the relative rotations.
    nEdges = size(edge_IDs,2);
    edge_errors = nan(1, nEdges);
    nPointsForEachEdge = nan(1, nEdges);
    nPointsForEachView = zeros(1, length(views.ID));

    RR = nan(3,3, nEdges);
    for i = 1:nEdges
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        r = edge_rs(:,i);
        
        intersection = intersect(views.points_ID{j}, views.points_ID{k});
        if (length(intersection) < nMinCovisiblePoints)
            error(['Something is wrong in the dataset! There is an edge with less than ', num2str(nMinCovisiblePoints), ', covisible points!'])
        end

        R_jk = reshape(r, [3,3]);
        
        theta = DegBtnRotations(R_jk, views.R{j}*views.R{k}');
        %disp(['edgeID = ', num2str(i), ', rot error = ', num2str(theta), ' deg.'])
        edge_errors(i) = theta;
        nPointsForEachEdge(i) = length(intersection);
        RR(:,:,i) = R_jk;

        nPointsForEachView(j) = nPointsForEachView(j) + length(intersection);
        nPointsForEachView(k) = nPointsForEachView(k) + length(intersection);
    end
    
    %% Additional output
    nViews = length(views.ID);
    nPoints = length(points.ID);
    nEdges = size(edge_IDs,2);

    edge_IDs_reverse = edge_IDs;
    edge_IDs_reverse(1,:) = edge_IDs(2,:);
    edge_IDs_reverse(2,:) = edge_IDs(1,:);
    
    
    %% Draw
    if (draw2verify)
        % Plot the reconstruction.
        figure;
        hold on
        scatter3(c_all(1,:),c_all(2,:),c_all(3,:))
        scatter3(p_all(1,:),p_all(2,:),p_all(3,:), 5, color_all', 'filled')
        axis equal
        set(gca,'Visible','off')

        % Plot 4 edges with one of their points.
        figure;
        for i = 1:4
            subplot(2,2,i)
            while (true)
                edgeID_to_plot = edge_IDs(:, randi(size(edge_IDs,2)));
                viewID1_to_plot = edgeID_to_plot(1);
                viewID2_to_plot = edgeID_to_plot(2);

                [pointIDs, point_idx1, point_idx2] = intersect(views.points_ID{viewID1_to_plot}, views.points_ID{viewID2_to_plot});
                if (~isempty(pointIDs))
                    point_idx = randi(length(pointIDs));
                    pointID_to_plot = pointIDs(point_idx);
                    point_idx_in_view1 = point_idx1(point_idx);
                    point_idx_in_view2 = point_idx2(point_idx);
                    break;
                end
            end

            hold on
            c1 = views.pos{viewID1_to_plot};
            c2 = views.pos{viewID2_to_plot};
            scatter3(c1(1),c1(2),c1(3), 'x')
            scatter3(c2(1),c2(2),c2(3), 'x')
            p = points.pos{:,pointID_to_plot};
            scatter3(p(1),p(2),p(3), 'o')

            uv1 = views.points_uv{viewID1_to_plot}(:,point_idx_in_view1);
            uv2 = views.points_uv{viewID2_to_plot}(:,point_idx_in_view2);
            f1 = views.R{viewID1_to_plot}'*(views.K{viewID1_to_plot}\[uv1;1]);
            f2 = views.R{viewID2_to_plot}'*(views.K{viewID2_to_plot}\[uv2;1]);
            endpoint1 = views.pos{viewID1_to_plot} + norm(c1-p)*f1/norm(f1);
            endpoint2 = views.pos{viewID2_to_plot} + norm(c2-p)*f2/norm(f2);

            line([c1(1), endpoint1(1)], [c1(2), endpoint1(2)], [c1(3), endpoint1(3)])
            line([c2(1), endpoint2(1)], [c2(2), endpoint2(2)], [c2(3), endpoint2(3)])
            title(['View IDs = ', num2str(viewID1_to_plot), ', ', num2str(viewID2_to_plot), ', point ID = ', num2str(pointID_to_plot)]);
        end
        
    end
    
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

function out = RandomRotation(angle_deg)
    axis = rand(3,1)-0.5;
    axis = axis/norm(axis);
    angle = angle_deg/180*pi;
    rotvec = angle*axis;
    out = ExpMap(rotvec);
end

function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
end

function out = ExpMap(in)
    angle = norm(in);
    axis = in/angle;
    so3 = SkewSymmetricMatrix(axis);
    R = eye(3)+so3*sin(angle)+so3^2*(1-cos(angle));
    out = R;
end
