%%TC_center_find.m

%% Created 18 Dec 2012, Dan Chavas

%% Purpose: To find the TC center using two methods based on the (perturbation) pressure field
%% 1) Simple: location of minimum pressure
%% 2) Centroid: center of mass (adjusted to focus on region near TC)
%% The code also quantifies the difference between the two, giving an error message if large

function [x_cent_COM,y_cent_COM,x_cent_SLP,y_cent_SLP] = TC_center_find(xvals,yvals,pp_pert);

dx = xvals(2)-xvals(1);
dy = yvals(2)-yvals(1);

%% Method #1: Minimum pressure point approach
[i_cent j_cent] = find(pp_pert == min(min(pp_pert)));
i_cent=i_cent(1);
j_cent=j_cent(1);
x_cent_SLP = xvals(i_cent);
y_cent_SLP = yvals(j_cent);

%% Method #2: Centroid method #1
%%Option 1
%%(comment in http://www.mathworks.com/matlabcentral/fileexchange/5457-centroid-calculation-function)
%{
    X_hist=sum(pp_pert,1);
    Y_hist=sum(pp_pert,2);
    X=1:length(X_hist); Y=1:length(Y_hist);
    centX_grid=sum(X.*X_hist)/sum(X_hist);
    centY_grid=sum(Y'.*Y_hist)/sum(Y_hist);
    x_cent_COM = xvals(floor(centX_grid))+dx*(centX_grid - floor(centX_grid));
    y_cent_COM = yvals(floor(centY_grid))+dy*(centY_grid - floor(centY_grid));
%}
%%Option 2
%%(http://www.mathworks.com/matlabcentral/newsreader/view_thread/253215)
A = pp_pert;
temp = -A;
A = temp-min(min(temp));    %need highest value to be at center of storm! use only positive values for simplicity
A(pp_pert>min(min(pp_pert))/4)=0;   %ignore pp_pert values less (in magnitude) than 1/4 the domain-minimum

%These next 4 lines produce a matrix C whose rows are
% the pixel coordinates of the image A
C=cellfun(@(n) 1:n, num2cell(size(A)),'uniformoutput',0);
[C{:}]=ndgrid(C{:});
C=cellfun(@(x) x(:), C,'uniformoutput',0);
C=[C{:}];

%This line computes a weighted average of all the pixel coordinates.
%The weight is proportional to the pixel value.
CenterOfMass=A(:).'*C/sum(A(:),'double');   %[row col]
x_cent_COM = xvals(floor(CenterOfMass(1)))+dx*(CenterOfMass(1) - floor(CenterOfMass(1)));
y_cent_COM = yvals(floor(CenterOfMass(2)))+dy*(CenterOfMass(2) - floor(CenterOfMass(2)));

%%TESTING
%{
    figure(9)
    imagesc(xvals,yvals,pp_pert')
    %imagesc(xvals,yvals,A')
    colorbar
    hold on
    plot(x_cent_COM,y_cent_COM,'gx')
    plot(x_cent_SLP,y_cent_SLP,'g+')
    set(gca,'YDir','normal')
    title('+ = min(p); x = centroid')
    axis([-100 100 -100 100])
        %NOTE: for data (but not xmat or ymat), x=rows, y=columns (i.e. NOT geographic!), so need transpose!
%}

%%Quantify difference between centroid and min(slp) estimates
x_cent_diff = x_cent_COM - x_cent_SLP;
y_cent_diff = y_cent_COM - y_cent_SLP;

%%If the two estimates are very large, warn you!
min(min(pp_pert))
if(min(min(pp_pert))<-1000)  %only apply if there exists a decent TC
    assert(abs(x_cent_diff)<20 & abs(y_cent_diff)<20,sprintf('WARNING! Centroid and SLP methods differ by %3.2f km',sqrt(x_cent_diff^2+y_cent_diff^2)))
end


end