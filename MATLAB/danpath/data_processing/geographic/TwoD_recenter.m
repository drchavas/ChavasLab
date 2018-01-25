%%TwoD_recenter.m

%%Created: 14 Dec 2012, Dan Chavas

%%Purpose: Takes in a 2d field of data and recenters it to the input (x,y) coordinates


function [xvals_out,yvals_out,data_out] = TwoD_recenter(xvals_in,yvals_in,data_in,x_cent,y_cent);
%{
%%FOR TESTING
clear
clc
close('all')
%}

%%USER INPUT -- FOR TESTING
%{
dx=1;
dy=1;
xvals_in = -10:dx:10;
yvals_in = -10:dy:20;
xvals_in = xvals_in+.5;   %off-center it for testing
yvals_in = yvals_in+.5;   %off-center it for testing
x_cent = 5; %(assume I've done this already)
y_cent = 3;
xmat = repmat(xvals_in,length(yvals_in),1);
ymat = repmat(yvals_in,length(xvals_in),1)';
%data_in=(10./(1+(xmat-x_cent).^2+(ymat-y_cent).^2)+randn(length(yvals_in),length(xvals_in)))'; %need to transpose -- input data is NOT geographic!
data_in=(10./(1+(xmat-x_cent).^2)+randn(length(yvals_in),length(xvals_in)))';   %need to transpose -- input data is NOT geographic!
%}
%%%%%%%%%%


%%Basic stuff
dx = xvals_in(2)-xvals_in(1);
dy = yvals_in(2)-yvals_in(1);
nx = length(xvals_in);
ny = length(yvals_in);

%%TESTING
%{
figure(1)
pcolor(xvals_in,yvals_in,data_in);
hold on
plot(x_cent,y_cent,'gx')
colorbar
%}

%% Find closest grid points to center point
i_closest_east = find(xvals_in>=x_cent,1);
i_closest_west = i_closest_east-1;
j_closest_north = find(yvals_in>=y_cent,1);
j_closest_south = j_closest_north-1;
if(i_closest_west>0)
    dx_closest_west = x_cent - xvals_in(i_closest_west);
else    %it is on the edge of the domain
    dx_closest_west = dx - (x_cent - xvals_in(i_closest_east));
end
if(j_closest_south>0)
    dy_closest_south = y_cent - yvals_in(j_closest_south);
else    %it is on the edge of the domain
    dy_closest_south = dy - (y_cent - yvals_in(j_closest_north));
end
        

%%plot(xvals_in(i_closest_west),yvals_in(j_closest_south),'g*')

%% Replicate the data in all directions
data_rep = repmat(data_in,3,3);

%% Re-center the data on (x_cent,y_cent)
nx_west = floor(nx/2);
ny_south = floor(ny/2);
i_west_rep = i_closest_west-(nx_west-1)+nx; %account for data added to the west
j_south_rep = j_closest_south-(ny_south-1)+ny; %account for data added to the south
i_east_rep = i_west_rep+nx-1;
j_north_rep = j_south_rep+ny-1;
data_out = data_rep(i_west_rep:i_east_rep,j_south_rep:j_north_rep);

%%TESTING
%{
xvals_add = xvals_in(end)-xvals_in(1)+(xvals_in(2)-xvals_in(1));
yvals_add = yvals_in(end)-yvals_in(1)+(yvals_in(2)-yvals_in(1));
xvals_rep = [xvals_in-xvals_add xvals_in xvals_in+xvals_add];
yvals_rep = [yvals_in-yvals_add yvals_in yvals_in+yvals_add];
xmat_rep = repmat(xvals_rep,length(yvals_rep),1);
ymat_rep = repmat(yvals_rep,length(xvals_rep),1)';
%nx_rep = length(xvals_rep);
%ny_rep = length(yvals_rep);
figure(3)
contourf(xmat_rep,ymat_rep,data_rep');
hold on
plot(xvals_rep(i_west_rep),yvals_rep(j_south_rep),'g*')
colorbar
%}

assert(size(data_in,1)==size(data_out,1))
assert(size(data_in,2)==size(data_out,2))

x_west_out = -1*((nx_west-1)*dx+dx_closest_west);
xvals_out = x_west_out:dx:x_west_out+(nx-1)*dx;
y_south_out = -1*((ny_south-1)*dy+dy_closest_south);
yvals_out = y_south_out:dy:y_south_out+(ny-1)*dy;

%%TESTING
%{
figure(2)
pcolor(xvals_out,yvals_out,data_out);
hold on
plot(0,0,'gx')
colorbar
%}

end
