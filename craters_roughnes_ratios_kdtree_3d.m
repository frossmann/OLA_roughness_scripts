
%
% Script to calculate interior-exterior roughness ratio of Bennu's impact
% craters using OLA ROI point clouds. Point clouds used here should have
% their .Intensity fields filled with measurements of delta_h(L).
%
% This script takes the RMS deviation of delta_h(L) values in the crater
% interior (interior_factor * radius) and compares it with the RMS
% deviation of delta_h(L) values in the crater exterior (ext_factor1 *
% radius to ext_factor2 * radius).



clear
close all
clc

% The source folder should contain crater ROI pointclouds with Intensity
% field equal to measurements of delta_h(L). 
pcd_source = '';

% Destination folder for results: 
output_folder = '';

%% Preliminaries:

% Iterate over baselines: 
baselines = [20:10:100, 200];

% set a target size to downsample large pointclouds:
target_size = 2e5;

% Set extent for interior and exterior regions
% interior region is interior_factor * radius
interior_factor = 0.6;
% exterior region is between [exterior_factors(1)*radius, exterior_factors(2)*radius]
exterior_factors = [2, 3];

for bl = 1:numel(baselines)
    tic
    baseline = baselines(bl);
    fprintf('Running for %i cm baseline...\n', baseline)
    
    % Point to my directory:
    dir_name = strcat(pcd_source, 'crater*_baseline',...
        num2str(baseline), 'cm*.pcd');
    d = dir(dir_name);
    
    
    % Load the Bierhaus 2020 crater catalog:
    T = readtable('../data/2020_06_bt_final_ma.xlsx');
    
    
    % Initialize an output variable:
    ratio_info = zeros([length(d), 7]);  % For saving a list of ratios
    
    % Loop through each file in the directory (roughness map) and find it's
    % inside - outside median roughness ratio with MADM variance.
    
    for id = 1:length(d)
        
        name = d(id).name;
        crater_name = extractBefore(name, '_');      % Name
        crater_id = str2double(crater_name(7:end));   % Index
        
        fprintf('\tProcessing %s...\n',d(id).name);
        filename = fullfile(d(id).folder,name);
        pc = pcread(filename);
        
        if pc.Count > target_size
            target_percent = target_size / pc.Count;
            pc = pcdownsample(pc, 'random', target_percent);
            fprintf('\t\tPointcloud downsampled to %i points (%4.2f %%)\n', pc.Count, target_percent*100);
        end
        
        
        % Extract info from the table:
        radius = T(crater_id,:).diameter_km_/2;
        center = [T(crater_id,:).x_km_,...
            T(crater_id,:).y_km_, T(crater_id,:).z_km_];
        
        
        [lon, lat, ~] = cart2sphd(pc.Location);
        [loc, lac, ~] = cart2sphd(center);
        k = dsearchn([lon, lat],[loc, lac]);
        nearest_center = pc.Location(k,:);
        
        % find points within radius of center for inside:
        [inside,~] = findNeighborsInRadius(pc,nearest_center,0.6*radius);
        inner_points = pc.Intensity(inside).*1e5;
        
        % find points within radius of center for outside:
        [in_outside,~] = findNeighborsInRadius(pc, nearest_center, 2*radius);  % these are indices!
        [in_limit,~] = findNeighborsInRadius(pc, nearest_center, 3*radius);
        outside = setdiff(in_limit, in_outside); % Set difference between 3R and 2R
        
        outer_points = pc.Intensity(outside).*1e5;
        
        % Calculate RMS Deviation of interior and exterior height deltas, and their ratiop
        rmsd_i = rms(inner_points, 'omitnan');
        rmsd_o= rms(outer_points, 'omitnan');
        rmsd_i_over_o = rmsd_i / rmsd_o;
        
        n_inside = numel(inner_points);
        n_outside= numel(outer_points);
        
        % Update the output array:
        ratio_info(id,:) = [crater_id,...    % 1: crater id
            radius*2000,...                  % 2: diameter in metres
            rmsd_i,...                       % 3: interior rmsd
            n_inside,...                     % 4: interior count
            rmsd_o,...                       % 5. exterior rmsd
            n_outside,...                    % 6. exterior count
            rmsd_i_over_o];                  % 7. interior/exterior ratio
        
        fprintf('\t\t%4.2f %% done.\n',100*id/numel(d));
        
    end
    % ========================================================================
    
    fprintf('Finished processing for %scm baseline.\n',num2str(baseline));
    toc
    
    
    savename = fullfile(output_folder,...
        sprintf('RRie_%scm_%s-%s-%sR_C.csv',...
        num2str(baseline),...
        num2str(interior_factor),...
        num2str(exterior_factors(1)),...
        num2str(exterior_factors(2))));
    writematrix(ratio_info, savename);
end
