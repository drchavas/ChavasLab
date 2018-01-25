%% 3D movie code for SAM

%-------------------------------------------
%*** do_3d_movie, ..._ffmpeg, movie_fps *** 
%-------------------------------------------
% do_3d_movie: flag (1 or 0) to make movie frames (or not)
% do_3d_movie_ffmpeg: flag (1 or 0) to compile movie frames (or not) into
%   .mp4 movie using ffmpeg
% movie_fps: string, number of frames per second to use for movie
%   compilation with ffmpeg (only relevant if do_3d_movie_ffmpeg =1)
%
%
do_3d_movie=1;
do_3d_movie_ffmpeg=0;
movie_fps = '15';

% figdir: wherever you want the movie frames to be saved
figdir='XXX';

% runid: dtring identifier for simulation -- used below to make frame names
runid='XXX';

%(N.B. -- I set figdir and runid dynamically, based on the title of the 3D
%files that I'm reading in)

%% Step 1: Initialization of movie frames
if do_3d_movie==1
    % size of movie frames, in pixels    
    xpix_size = 960;
    ypix_size = 640;
    
    % can modify these position settings -- they should give a figure
    % that is 100 pixels up and to the right of the lower left corner
    % of the monitor, and which is xpix_size pixels wide, ypix_size high
    hmv3d=figure('Position',[100 100 xpix_size ypix_size],'Color','k');
    
    % use the clock function to get the time at which the movie is
    % initiated; then process to get a string that will serve as a unique
    % id for the movie frames folder
    date_id=clock;
    % format YYYY-MM-DD_hh'h'mm'm'
    date_str=sprintf('%04d-%02d-%02d_%02dh%02dm',date_id(1:5));
    % set up the name of the directory based on figdir and the additional
    % date string -- also could be modified so that the folder name 
    framedir=[figdir 'mframes_3d_' date_str '\'];
    
    if ~exist(framedir,'dir')
        mkdir(framedir);
    end
end

%% Step 2: loop to read in data over multiple timesteps
%-- put frame-generating code below in this loop
%
% non-local variables used within the block are:
%       x3,y3,z3    = 1D arrays of x,y,z (km)
%       qn3         = 3D array of nonprecipitating water (g/kg)
%       NX_GL       = length(x3)
%       NY_GL       = length(y3)
%       tabs3       = 3D array of air temperatures (K)
%       DX,DY       = grid spacing in x and y (m)
%       lst_3d_dd   = time in days
%       lst_3d_hh   = time in hours
%       lst_3d_mm   = time in minutes
%       i           = integer, frame number (or timestep), must start from 0 for ffmpeg to work
%

% Make 3D movie frames
if do_3d_movie==1
    % clear the existing frame
    clf;
    
    % *** CLOUD GRAPHICS BLOCK ***
    % set the isolevel of qn (cloud water/ice) that will be plotted
    qn_isolevel=0.3;
    % set pgfx3d as the handle for the graphics object given by
    % patch(isosurface(...)) -- this is the line that makes the cloud
    % graphics object
    pgfx3d = patch(isosurface(x3,y3,z3,permute(qn3,[2 1 3]),qn_isolevel));
    % set the view angle to have -26 degrees azimuth, and 20 degrees
    % elevation
    view(-26,20);
    % set the axis background to black -- allows for visualization of clouds
    whitebg('k');
    % set lighting properties of the cloud graphics object
    % FaceColor = color of cloud polygon faces; EdgeColor= color of cloud edges
    % AmbientStrength = strength of illumination of clouds by ambient (non-directional) light
    % FaceAlpha = sets the transparency of the cloud surface (1=opaque; 0=transparent)
    set(pgfx3d,'FaceColor','white','EdgeColor','none','AmbientStrength',0.5,'FaceAlpha',0.4);
    % camlight sets a point light source to exist right and up from the camera position
    camlight;
    % lighting phong sets the reflection model for the cloud patch
    % object and other surfaces that will respond to illumination
    % (e.g., the surface contour below)
    lighting phong;
    
    % *** SURFACE AIR T BLOCK ***
    hold on;
    % plot surface air temperature (lowest model layer) on z=0 surface --
    % (might be better to convert to deg. C or F for a general audience??)
    surfc(x3,y3,zeros(NX_GL,NY_GL)',double(tabs3(:,:,1)'),'EdgeColor','none','AmbientStrength',0.5);
    %surfc(x3,y3,zeros(NX_GL,NY_GL)',double(w3(:,:,2)'),'EdgeColor','none','AmbientStrength',0.5);
    
    % *** AXES BLOCK ***
    % x, y, and z-limits
    xlim([0 max(x3)+DX/1000]);
    ylim([0 max(y3)+DX/1000]);
    zlim([0 16]);
    % x, y, and z- labels
    xlabel(' x, km ','fontsize',16);
    ylabel(' y, km ','fontsize',16);
    zlabel(' height, km ','fontsize',16);
    % limits of colorscale
    caxis([288 294]);
    %caxis([-0.2 0.2]);
    % colorbar properties (position is normalized by overall graphics frame dimensions)
    colorbar('Position',[0.92 0.2 0.03 0.6]);
    % data aspect ratio -- vertical scale exaggeration is equal to
    % daspect(1)/daspect(3), so for daspect([1 1 0.25]), this is 4
    daspect([1 1 0.25]);
    
    % *** TITLE ***
    %titlestring={['Isosurface of Condensate=' sprintf('%2.1f',qn_isolevel) 'g/kg and w at 100 m (m/s)'];
    titlestring={['Isosurface of Condensate=' sprintf('%2.1f',qn_isolevel) 'g/kg and Surface Air T (K)'];...
        ['Day ' sprintf('%03d',lst_3d_dd) ', Local Solar Time : ' sprintf('%02d:%02d ',lst_3d_hh,lst_3d_mm)]};
    
    title(titlestring,'fontsize',16);
    
    % *** ISLAND RECTANGLE ***
    % define island outline with 'line' commands
    line([192 192],[192 288],[0 0],'Color','k','LineWidth',1.5); % west edge
    line([192 288],[192 192],[0 0],'Color','k','LineWidth',1.5); % south edge
    line([288 288],[192 288],[0 0],'Color','k','LineWidth',1.5); % east edge
    line([192 288],[288 288],[0 0],'Color','k','LineWidth',1.5); % north edge
    
    % *** NOTES AT BOTTOM ***
    % make fontsize of axes slightly larger than default
    set(gca,'fontsize',12);
    text(-150, -150, {'Note: Vertical Scale Stretched \times 4 Relative to Horizontal Scale';...
        'Island Outlined by Black Rectangle'});
    
    % *** SAVE FRAME ***
    %frameid=sprintf('_3d_mframe_%08d.png',ninit3d+(i-1)*nsave3d);
    frametmpl='_3d_frame_no-%05d.png';
    frameid=sprintf(frametmpl,i);
    framename=[framedir runid frameid];
    set(gcf,'inverthardcopy', 'off');
    print(hmv3d,'-dpng',framename);
    
    % *** COMPILING FRAMES ***
    % This is done at the end of the file, after the main loop -- search
    % below for 'movie'
    
end

%% Step 3: Compile Movie Frames with ffmpeg
if (do_3d_movie==1)&&(do_3d_movie_ffmpeg==1)
    % Here, I use ffmpeg to compile the movie frames into an .mp4 movie
    % this invokes the system shell (via |system(command)|, 
    % rather than only using MATLAB commands, and requires that ffmpeg
    % exists on your system
    
    % Set movie names for large and small versions
    moviename1 = [runid '_3d_cloudmovie_' sprintf('%03dmin-%sfps.mp4',round(24*60/nfilesperday),movie_fps)];
    moviename2 = [runid '_3d_cloudmovie_' sprintf('%03dmin-res1-%sfps.mp4',round(24*60/nfilesperday),movie_fps)];
    
    % cd to the directory with the movie frames
    cd(framedir);
    
    % set the first ffmpeg command - compiles the frames into a large
    % movie, ideally with no quality loss (-sameq)
    ffmpeg_command1 = ['ffmpeg -r ' movie_fps ' -i ' runid frametmpl ' -sameq -vcodec mpeg4 ' moviename1];
    system(ffmpeg_command1);

    % then rescale/filter this first movie, which should shrink the file size a lot, 
    % but not lose much apparent quality  
    ffmpeg_command2 = ['ffmpeg -i ' moviename1 ' -vf scale=iw/1.5:-1 ' moviename2];
    system(ffmpeg_command2);
end



