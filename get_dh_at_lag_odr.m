%%
% Loops through query points (nodes).
% For each node, run a fixed radius neighbour search.
% Find points X within f*L of the node:
%  - f: detrending radius factor
%  - L: baseline
%
% Calculates orthogonal height of points X w.r.t best fit plane using
% orthogonal distance regression.
% differences between points separated by a lag distance `lag_radius`. Point
% heights are determined as the height above the best fit plane calculated
% using orthogonal distance regression (see paper below) using points within
% a local radius `detrend_radius`.
%
% Ryan M. Pollyea, Jerry P. Fairley; Estimating surface roughness of
%   terrestrial laser scan data using orthogonal distance regression.
%   Geology 2011;; 39 (7): 623–626. doi: https://doi.org/10.1130/G32078.1
%

function output = get_dh_at_lag_odr(X_search, X_query, baseline, detrend_factor)

% identify the detrending radius to use around query points: 
detrending_radius = detrend_factor * baseline;

% transform array X_search into a pointCloud object: 
pcd_search = pointCloud(X_search);

% identify the number of query points (nodes):
n_query = size(X_query, 1);

% preallocate an array to store [x, y, z, dh, nn]:
output = zeros([n_query, 5]);

% loop through all query points: 
for jj = 1:n_query

    % identify this iteration's query point: 
    node = X_query(jj, :);

    % find points in the search pointCloud within 1 detrending radius of the node:     
    [X_idx, X_dists] = findNeighborsInRadius(pcd_search, node, detrending_radius);


    % extract the neighbours in radius: (x, y, z)
    Xn = X_search(X_idx,:);

    % calculate ODR heights of all points used to detrend:
    Xn_heights = calculate_odr_heights_xyz([Xn; node]);


    % identify points within 1L +/- 5% of the node
    pad_percent = 0.05;
    target_neighbour_idx = X_dists < (baseline + pad_percent*baseline) & ...
        X_dists > (baseline - pad_percent* baseline);

    % if any points within 1L +/- 5% of the node: 
    if any(target_neighbour_idx)
        h1 = Xn_heights(end);  % height of node
        h2s = Xn_heights(target_neighbour_idx);  % height of first neighbour
        dh = abs(h1 - h2s(1));  % height difference (absolute)
        output(jj, :) = [node, dh, 1];  % save [x, y, z, dh, nn]
        % nn is number of points used to calculate dh, in this case always
        % 1 or 0. 

    else
           continue
    end

end

end



