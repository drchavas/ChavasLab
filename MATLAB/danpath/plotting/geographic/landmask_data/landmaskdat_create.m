%%landmaskdat_create.m

%% 5 Mar 2014, Dan Chavas

coastal_res = 10;   %output resolution [pts/deg]
make_plot = 0;
[isLand,lon_grid,lat_grid] = land_mask_global(coastal_res,make_plot);
save landmaskdat isLand lon_grid lat_grid