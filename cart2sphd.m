function [LON, LAT, R] = cart2sphd(XYZ)
% Runs MATLAB's cart2sph() but outputs in degrees instead of radians to
% save some keystrokes. 
% 
% Input:
% - XYZ: an [Nx3] matrix in cartesian coordinates. 

    [LON, LAT, R] = cart2sph(XYZ(:,1), XYZ(:,2), XYZ(:,3));
    LON = rad2deg(LON);
    LON(LON < 0) = LON(LON<0) + 360;
    LAT = rad2deg(LAT);
end