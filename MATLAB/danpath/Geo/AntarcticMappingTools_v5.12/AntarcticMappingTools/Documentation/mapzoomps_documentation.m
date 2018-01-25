%% |mapzoomps| documentation
% |mapzoomps| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017). Click <List_of_functions.html here>
% for a complete list of functions in AMT. 
% 
% The |mapzoomps| function zooms a south polar stereographic map to a specified location and extent.
% This is an adaptation of |mapzoom|, but does not require Matlab's Mapping Toolbox. 
% Syntax for |mapzoomps| is similar to |mapzoom|, but differs slightly for some options.
% 
%% Syntax
% 
%  mapzoomps
%  mapzoomps(lat,lon) 
%  mapzoomps(x,y)
%  mapzoomps('SCAR location') 
%  mapzoomps(...,'size',mapsizekm)
%  mapzoomps(...,InsetLocation)
%  mapzoomps(...,'insetsize',sizefraction)
%  mapzoomps(...,'frame','off')
%  mapzoomps(...,'km')
%  h = mapzoomps(...) 
% 
%% Description 
% 
% |mapzoomps(lat,lon)| centers a 500 km wide map about the georeferenced 
% location given by |lat, lon|. 
% 
% |mapzoomps(x,y)| centers a 500 km by 500 km map about the polar stereographic 
% eastings and northings |x| and |y|.  Polar stereographic coordinates are
% automatically determined by the |islatlon| function. 
% 
% |mapzoomps('SCAR location')| uses <scarloc_documentation.html |scarloc|> to find the coordinates corresponding
% to the string |'SCAR location'|. 
% 
% |mapzoomps(...,'size',mapsizekm)| specifies size of the map in kilometers given 
% |mapsizekm|, which can be a scalar to create a square map or a two-element array
% for a rectangular map in the form |[mapwidthkm mapheightkm]|, where |mapwidthkm| and
% |mapheightkm| are the dimensions of the map in kilometers. 
%
% |mapzoomps(...,InsetLocation)| creates an inset map at the location |InsetLocation|,
% which can be 
% 
% * |'southeast'| or |'se'| (lower right)
% * |'northwest'| or |'nw'|  (upper left) 
% * |'northeast'| or |'ne'|  (upper right)
% * |'southwest'| or |'sw'|  (lower left) 
%
% |mapzoomps(...,'insetsize',sizefraction)| specifies size of the inset as a
% fraction of the width of the current map. Default |sizefraction| is |0.25|. 
%
% |mapzoomps(...,'frame','off')| removes frame from the inset. 
%
% |mapzoomps(...,'km')| is for plots in polar stereographic kilometers rather than the default meters.
% 
% |h = mapzoomps(...)| returns a handle |h| of inset map axes. 
% 
%% Example 1 
% Initialize a 500 km by 500 km map centered on Pitman Fracture Zone, plot <http://www.mathworks.com/matlabcentral/fileexchange/52333-ibcso-international-bathymetric-chart-of-the-southern-ocean-toolbox/content/IBCSO/documentation/html/ibcso_documentation.html |ibcso|> 
% bathymetry, and place a graphical reference scale with <scalebarps_documentation.html |scalebarps|>.  When calling |mapzoomps|
% specify |'ne'| to place an inset map in the upper right-hand corner. 

mapzoomps('pitman fracture zone','ne')
ibcso('image','xy')
scalebarps
axis image % tightens up null space

%% Example 2 
% Plot Bedmap2 bed elevation, apply relief shading with <shadem_documentation.html |shadem|>, and center a 
% 1600 km wide by 800 km tall map on the Lambert Glacier/Amery Ice Shelf system at (71°S,69°N).  Place 
% frameless inset map in the upper left hand corner: 

figure
bedmap2('bed','xy') 
shadem([225 80],3)
mapzoomps(-71,69,'size',[1600 800],'nw','frame','off')
scalebarps
scarlabel({'Lambert Glacier','Cape Child','Landing Bluff'},...
   'fontangle','italic','fontsize',8)

%% Fixing a misaligned inset map 
% Sometimes the inset gets misaligned. Here's an example: 

figure
mapzoomps('pine island glacier','inset','ne') 
modismoaps
axis off 

%%
% So far, so good. But then you try to add something with a colorbar and the inset gets misaligned: 

measuresps('speed','alpha',0.5) 

cb = colorbar; 
ylabel(cb,' ice speed (m/a) ') 
cmocean amp

%% 
% Above I used <https://www.mathworks.com/matlabcentral/fileexchange/57773 |cmocean|> (Thyng et al., 2016) to make 
% the red colormap. To fix the the inset misalignment, simply call |mapzooms| again, which will delete the old inset and try
% to place a new one in the correct spot: 

mapzoomps('pine island glacier','inset','ne') 

%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info 
% The |mapzoomps| function and supporting documentation were written for Antarctic Mapping Tools
% by Chad A. Greene of the University of Texas at Austin's Institute for Geophysics (UTIG) in 
% February 2016. http://www.chadagreene.com