% Script to create ROI pointclouds for a selection of Bennu's craters. 
%

% Set source of pointcloud(s) that the ROI(s) will be cropped from.
%
% This could be a folder of pointclouds which only have spatial information
% (XYZ) or pointclouds that include Intensity values for delta_h(L) or the
% local RMSD of delta_h(L:  
%
pcd_source = './';
% glob = {'*_ola_odr2L_height_deltas_baseline', 'cm*.pcd'};  % usage: strcat(glob(1), baseline, glob(2))


%  Set destination for output: 
pcd_destination = './';

% Load the crater catalog: 
crater_T_fn = '../data/2020_06_bt_final_ma.xlsx';
crater_T = readtable(crater_T_fn, 'VariableNamingRule','preserve');
%%
% *************** ROI EXTRACTION FOR ROUGHNESS RATIOS ****************
% =========================================================================
% % Load the information for craters that were used in Bierhaus et al.,
% % (2023)
% bierhaus2023_crater_T_fn = '../data/bierhaus2023_craters_gt5m_lt80m_olav21_roi_flags.xlsx;
% bierhaus2023_crater_T = readtable(bierhaus2023_crater_T_fn, 'VariableNamingRule','preserve');
% 
% %  Find information for craters which passed QC: 
% idx = find(bierhaus2023_crater_T.flags == 0);
% crater_ids = bierhaus2023_crater_T.id(idx);
% =========================================================================



% ********** ROI EXTRACTION FOR ROUGHNESS PROFILES *************
% =========================================================================
% Find craters with a diameter > 80 m
min_diameter = 80 / 1000;  % km
idx = find(crater_T.('diameter [km]') > min_diameter);
crater_ids = idx;
% =========================================================================  

% Isolate a new table
candidates = crater_T(crater_ids,:);
diameters = candidates.('diameter [km]');
% Plot histograms of diameter, longitude, latitude:
figure;
subplot(311);
histogram(1000*diameters,'normalization','pdf');
xlabel('m')
title('diameter');
subplot(312);
histogram(candidates.('longitude [deg]'),'binwidth', 15,'normalization','pdf');
xlabel('deg');
ylabel('PDF');
title('longitude [deg]');
subplot(313)
histogram(candidates.('latitude [deg]'),'binwidth', 10,'normalization','pdf');
xlabel('deg');
title('latitude [deg]');

% Plot as global scatter
figure;
f1 = scatter(candidates.('longitude [deg]'), candidates.('latitude [deg]'),...
    10, 1000*diameters,...
    'filled');
% f1.MarkerFaceAlpha = 0.5;
cbar = colorbar;
grid on
title(cbar,'diameter [m]');
title('Global distribution (not to scale)');
xlabel('longitude');
ylabel('latitude');
axis equal
xlim([0 360])
ylim([-90, 90])

[F, X] = ecdf(diameters);
figure;
plot(X*1000, F,'linewidth',3);
grid on
xlabel('Diameter [m]');
ylabel('ECDF');
title('ECDF of diameters');

%% Loop through baselines

baselines = [20:10:100, 200:100:1000, 2000];

for jj = 1:numel(baselines)
    baseline = baselines(jj);
    
   
    pc_dir_name = strcat(pcd_source, glob{1}, num2str(baseline), glob{2});
    pc_dir = dir(pc_dir_name);
    
    filename = fullfile(pc_dir.folder, pc_dir.name);
    
    [folder, file, ext] = fileparts(filename);
    fprintf('Loading %s...\n',strcat(file, ext));
    tic
    pc = pcread(filename);
    toc
    
    
    [lon, lat, ~] = cart2sphd(pc.Location);

    
    total_time = 0;
    loop_time = tic;
    for ii = 1:numel(crater_ids)
        loop_start = tic;
        ID = crater_ids(ii);
        fprintf('    crater %i (%i of %i).', ID, ii, numel(crater_ids));
        
        savename = fullfile(pcd_destination, strcat('crater',num2str(ID),'_', file, ext));
        
        if isfile(savename)
            fprintf('    File already exists).\n');
            continue
        end
        
        
        fprintf('\n');
        crater = candidates(ii, :);
        diameter = crater.('diameter [km]');
        fprintf('        diameter: %8.8f m\n', diameter * 1000)
        
%          if 1e5*diameter < baseline / 0.2
%              fprintf('(    Skipping for size constraint).\n')
%            continue 
%         end
        
        radius = diameter/2;
        center = [crater.('x [km]'), crater.('y [km]'), crater.('z [km]')];

        % Find the nearest point (in lat/lon) to the center coordinate of
        % the crater. This will be the centre point of the ROI:
        [loc, lac, ~] = cart2sphd(center);
        k = dsearchn([lon, lat],[loc, lac]);
        nearest_center = pc.Location(k,:);

        % Find points inside of 3 crater radii.
        fprintf('        searching tree...\n')
        [indices, ~] = findNeighborsInRadius(pc, nearest_center, 3*radius);

        
        fprintf('        selecting...\n')
        pc_out = select(pc, indices);
        
        % Save the circular ROI:
        fprintf('        writing cloud...\n')
        pcwrite(pc_out, savename, encoding='binary');
        loop_end = toc(loop_start);
        total_time = total_time + loop_end;
        
        %         fprintf('%4.2f %% done.\n',100*ii/numel(craterIDs));
    end  % end of for crater in craters
    fprintf('    Finished processing for %i cm baseline.\n',baseline);
    fprintf('Elapsed time in loop is %4.2fs\n',toc(loop_time));
end  % end of for baseline in baselines
