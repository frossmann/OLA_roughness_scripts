%
% Script to calculate the RMS deviation in local circular patches about
% Bennu's surface globally.


% Path to folder containing pointClouds with .Intensity field set to 
% delta_h(L) 
source_folder = '';

d = dir(fullfile(source_folder, '*.pcd'));

destination_folder = source_folder; 

min_radius = 2 / 1000;


%%
for ii = 1:numel(d)
    
    filename = fullfile(d(ii).folder,d(ii).name);
    
    [root, name, ext] = fileparts(filename);
    fprintf('\nLoading %s... ', name);
    
    baseline = str2double(extractBefore(extractAfter(name, 'baseline'), 'cm')) / 100 / 1000;

    radius = max(baseline, min_radius);
    % radius = min_radius;

    voxel_size_coarse = radius;
    
    prefix = strcat('rms_deviation_local', num2str(radius*1000), 'm_voxel', num2str(voxel_size_coarse*1000),'m_');
    pcd_savename = fullfile(destination_folder, strcat(prefix,name, ext));
    mat_savename = fullfile(destination_folder, strcat(prefix,name, '.mat'));
    
%     if isfile(pcd_savename) & isfile(mat_savename);
%         fprintf('File already exists.\n');
%         continue 
%     end
    
    pc_full = pcread(filename);
    fprintf('Done. \n')
  
    
    fprintf('Filtering... ');
    keep = (vecnorm(pc_full.Location,2,2) > 0.2) & ~isnan(pc_full.Intensity);
    pc_in = select(pc_full, keep);
    fprintf('Done.\n');
    %% Downsample
    
    fprintf('Downsampling... ');
    pc_down = pcdownsample(pc_in, gridAverage=voxel_size_coarse);
    fprintf('Done.\n');
    
    %% Accumulate points and calculate RMSD
    fprintf('Searching neighbours...\n');
    tic;
    [pcd_nodes, n_neighbours] = radial_search_RMSD(pc_in, pc_down, radius);
    toc;
    fprintf('Done.\n ');
    

    xyz = pcd_nodes.Location;
    rmsd = pcd_nodes.Intensity;

    %% Show the point cloud: 
    pcd_show = copy(pcd_nodes);
    pcd_show.Intensity = pcd_show.Intensity .* 1e5;
    
    figure;
    r = min(vecnorm(pcd_show.Location,2,2));
    [x,y,z] = sphere;
    hs1 = surf(x*r,y*r,z*r);
    set(hs1, 'FaceColor', [0, 0, 0], linestyle='none')
    axis equal
    hold on;
    pcshow2(pcd_show);
    % clim([prctile(pcd_show.Intensity, 1), prctile(pcd_show.Intensity, 99)]);
    colorbar(color='w');
    title(sprintf('%s\n%s', prefix, name), color='w', interpreter='none')
    
    %% Save the point cloud: 
    fprintf('Writing .pcd... ');
    pcwrite(pcd_nodes, pcd_savename, encoding='binary');
    fprintf('Done.\n ');

    fprintf('Writing .mat... ');
    save(mat_savename, 'xyz', 'rmsd', 'n_neighbours');

    fprintf('Done.\n ');
    
end
