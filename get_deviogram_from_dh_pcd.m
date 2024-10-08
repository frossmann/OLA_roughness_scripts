%
% Calculates the deviogram of RMS deviation (nu(L)) as a function of
% delta_h(L). 
%

clear

% path to pointClouds with delta_h(L) as Intensity values. 
path_to_clouds = '';
d = dir(fullfile(path_to_clouds, '*.pcd'));

% output results into the data's parent folder 
output_path = fullfile(d(1).folder);

% empty arrays to hold L and RMSD
baselines = zeros(numel(d), 1);
rmsds = zeros(numel(d), 1);

% Loop through files and calculate RMSD(pc.Intensity) for each: 
for ii = 1:numel(d)

    filename = fullfile(d(ii).folder, d(ii).name);
    
    [~, name, ext] = fileparts(filename);
    fprintf('%s\n', name);
    
    baseline = str2double(extractBefore(extractAfter(name, 'baseline'), 'cm'));

    pc = pcread(filename);
    
    % save baseline and RMSD
    baselines(ii) = baseline;
    rmsds(ii) = rms(pc.Intensity);
    
end


% sort the entries by ascending baseline: 
[sorted_baselines, sorted_idx] = sort(baselines);
sorted_rmsds = rmsds(sorted_idx);

L_m = sorted_baselines ./ 100;  % baselines in m
rmsd_m = sorted_rmsds * 1000; % rmsd in m

%% Show the deviogram:

figure;
plot(L_m, rmsd_m)
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
grid on;
xlim([L_m(1) * 0.1, L_m(end) * 10])
ylim([rmsd_m(1) * 0.1, rmsd_m(end) * 10])
xlabel('L, m')
ylabel('RMSD, m')


%% Save the deviogram data: 
% save('ola_l2a_v21_odr2L_dh_L0.20m-20m_rmsd.mat', 'rmsd_m', "L_m")