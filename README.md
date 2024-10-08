# metre_scale_roughness_of_asteroid_bennu_from_OLA

## Files: 
- `cart2sphd`: function, cartesian to spherical coordinates output in units of degrees.
- `craters_roughness_ratios_kdtree_3d`: script to calculate interior/exterior roughness ratio ($RR_{i/e} = \nu(L)_{interior} / \nu(L)_{exterior}$) for crater ROIs.
- `craters_roughness_profiles_kdtree_3d`: script to calculate radial roughness profiles ($\nu(L)$ vs $R_C$) for crater ROIs.
- `crop_crater_roi_pcd`: script to create crater ROIs from pre-existing point clouds.
- `get_deviogram_from_dh_pcd` script to compute deviogram ($\nu(L)$ vs $L$) using point clouds.
- `get_dh_at_lag_odr`: function, calculates differences in (detrended) height between points separated by a lag distance, using orthogonal distance regression to define height and some detrending radius ( `= 2L` for paper).
- `make_global_rmsd_pcd`: script to generate point cloud of RMS deviation roughness from point clouds with $\Delta h(L)$ as `Intensity` field.

## Data
View only link to data repository containg $\Delta h(L)$ pointclouds for $L$ = 0.2 - 20.0 m: https://zenodo.org/records/13900222?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImFiMzEwOWQyLTNmNDUtNGIxMS1iNDYxLTJhZDU5OTdlNTEyZiIsImRhdGEiOnt9LCJyYW5kb20iOiIwMTIyNjQ4ZjNiZGI4MjZjMjg0NmM4OGFjZmU4OTYxOCJ9.qXW8d94elfJvC-vQzfiLrYtPVmUsNTRIXOFaZ-pmXwwytW3vjVJ81blwFtGh1WSUnelSzqnuz4UIej5b58KI0Q
