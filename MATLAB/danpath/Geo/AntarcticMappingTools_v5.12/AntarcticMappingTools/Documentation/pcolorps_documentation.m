%% |pcolorps| documentation
% |pcolorps| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017). Click <List_of_functions.html here>
% for a complete list of functions in AMT. 
% 
% Don't have Matlab's Mapping Toolbox?  No problem.  |pcolorps| works just like Matlab's |pcolor| or
% |pcolorm| functions, but plots georeferenced data in Antarctic polar stereographic coordinates (true latitude 71°S).
% 
%% Syntax
% 
%  pcolorps(lat,lon,Z)
%  pcolorps(...,'PropertyName',PropertyValue,...)
%  pcolorps(...,'km') 
%  h = pcolorps(...)
% 
%% Description 
% 
% |pcolorps(lat,lon,Z)| constructs a surface to represent the data grid |Z| 
% corresponding to a georeferenced |lat,lon| grid in South polar stereographic
% eastings and northings. 
% 
% |pcolorps(...,'PropertyName',PropertyValue,...)| specifies any number of
% surface properties. 
% 
% |pcolorps(...,'km')| plots in polar stereographic kilometers instead of the
% default meters. 
%     
% |h = pcolorps(...)| returns a column vector of handles to surface objects.
% 
%% Example
% Plot gridded georeferenced <http://www.mathworks.com/matlabcentral/fileexchange/42353 Bedmap2> data as a pseudocolor plot in cartesian polar stereographic coordinates:   

[lat,lon,bed] = bedmap2_data('bed','resolution','5 km'); 

pcolorps(lat,lon,bed)
xlabel('eastings (m)') 
ylabel('northings (m)')

%% 
% Enhance the appearance of topographic relief with <shadem_documentation.html |shadem|>.  

shadem(3) 

%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info
% This function was written by <http://www.chadagreene.com Chad Greene> of the University of Texas
% Institute for Geophysics (UTIG), February 2015, for inclusion in the
% Antarctic Mapping Tools package. Updated July 2015 to allow plotting in polar stereographic kilometers.