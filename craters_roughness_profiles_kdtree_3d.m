% 
% Script to calculate radial roughness profiles of Bennu's impact
% craters using OLA ROI point clouds. Point clouds used here should have
% their .Intensity fields filled with measurements of RMSD(L).


clear
clc


SAVE = true;

% Path to folder containing crater ROI point clouds with Intensity values
% set to the local RMS deviation.
source_folder = '';

% Where to save the results:
destination_folder = '';

% Bootstrap confidence intervals in parallel
optset = statset('UseParallel',true);
n_boot = 1000;  % bootstrap 1000 iterations
%%

% Load the craters spreadsheet from Bierhaus 2023 (subsurface),
T = readtable('../data/2020_06_bt_final_ma.xlsx','VariableNamingRule','preserve');

% Pick a roughness baseline (or baselines):
baselines = [20:10:100, 200:100:1000, 2000];


% Set up the framework for roughness profiles:
dr_percent = 0.1;  % each radial bin is 10% of the crater radius
annuli_edge_percents = [0, 0.2:dr_percent:3];  % initialize bin edges
% find bin centers
annuli_center_percents = annuli_edge_percents(1:end-1) + diff(annuli_edge_percents) / 2;
n_bins = numel(annuli_center_percents);  % number of bins


fprintf('starting...\n')
for bb = 1:numel(baselines)  % Loop over baselines
    
    tic
    
    baseline = baselines(bb);
    fprintf('batching for %icm...\n', baseline)

    dir_name = strcat(...
        source_folder,...
        '*crater*baseline',num2str(baseline),...
        'cm*.pcd');
    d = dir(dir_name);
    
    % Initialize output array:
    %     profiles = zeros([numel(d) - size(outliers,1), 3*n_bins + 1]);
    profiles = zeros([numel(d), 3*n_bins + 1]);
    
    
    % Loop over craters and make roughness profiles:
    for ij = 1:length(d)
        
        crater_name = extractBefore(d(ij).name,'_');
        crater_id = str2double(crater_name(7:end));
        fprintf('    crater %i:  %i of %i\n', crater_id, ij, length(d));
        
        % Load the pointcloud
        pc = pcread(fullfile(d(ij).folder, d(ij).name));
        
        % Extract info from the table:
        radius = T(crater_id,:).('diameter [km]')/2;  % radius in km
        %             loc = T(crater_id,:).('longitude [deg]');  % crater center longitude
        %             lac = T(crater_id,:).('latitude [deg]');  % crater center latitude
        % center = mean(pc.Location);
        center = [T(crater_id, :).('x [km]'),...
            T(crater_id, :).('y [km]'),...
            T(crater_id, :).('z [km]')];
        
        
        % Find closest point of pointcloud in spherical coordinates to
        % the listed center in the table:
        [loc, lac, ~] = cart2sphd(center);
        [lon, lat, ~] = cart2sphd(pc.Location);
        k = dsearchn([lon, lat],[loc, lac]);
        nearest_center = pc.Location(k,:);
        
        % Set up the roughness profiles:
        edge_radii = radius .* annuli_edge_percents;
        
        % Initialize some variables to hold median roughness as a function of
        % radial distance from crater centre, and it's variability.
        med_roughness = zeros(1, n_bins);
        med_roughness_ci_lower = zeros(1, n_bins);
        med_roughness_ci_upper = zeros(1, n_bins);
        
        % Loop through concentric bins from the inside-outwards:
        parfor ii = 1:n_bins
            
            in_inside = findNeighborsInRadius(pc,nearest_center,edge_radii(ii));  % in inner ring
            in_outside = findNeighborsInRadius(pc,nearest_center,edge_radii(ii+1)); % in outer ring
            in_between = setdiff(in_outside, in_inside);  % in between
            
            sample =  pc.Intensity(in_between) .* 1e5;
            
            %             fprintf('N = %i\n',numel(inBetween));
            %             pcshow(theSample);
            %             pause
            
            % NaN value for annuli containing fewer than 5 points:
            if numel(sample) < 5
                med_roughness(ii) = NaN;
                median_roughness_ci_lower(ii) = NaN;
            else
                % save the median:
                med_roughness(ii) = median(sample, 'omitnan');
                % approximate the uncertainty in the median:
                confidence_interval = bootci(n_boot, {@median, sample}, 'Options', optset);
                med_roughness_ci_lower(ii) = abs(confidence_interval(1) - med_roughness(ii));
                med_roughness_ci_upper(ii) = abs(confidence_interval(2) - med_roughness(ii));
                %                 med_roughness_ci_lower(ii) = NaN;
                %                 med_roughness_ci_upper(ii) = NaN;
               
            end
        end  % end of radial histogram
        
        
        % Pack everything into one output array:
        profiles(ij,1) =  crater_id;
        profiles(ij,2:n_bins + 1) = med_roughness;
        profiles(ij, (n_bins + 2):(2*n_bins + 1)) = med_roughness_ci_lower;
        profiles(ij, (2*n_bins + 2):(3*n_bins + 1)) = med_roughness_ci_upper;
        
%         fprintf('\t%2.1f %% done.\n',100*ij/numel(d))
    end % end of for crater in craters
    
    if SAVE
        % Save the output variable:
        the_date = string(datetime(now,'ConvertFrom','datenum','Format','yyyy-MM-dd'));
        savename = strcat(destination_folder, the_date,...
            '_roughness_profiles_baseline',...
            num2str(baseline),'cm.mat');
        save(savename,'profiles')
    end
    
    % Update:
    fprintf('Finished processing for %s cm baseline\n',num2str(baseline));
    toc
end  % end of for baseline in baselines.
