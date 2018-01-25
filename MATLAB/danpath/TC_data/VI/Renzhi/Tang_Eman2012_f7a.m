% TRYING TO RE-CREATE THE TANG AND EMANUEL (2012) FIGURE 7A USING DV 
% FOR 24-HR (AS IN ORIGINAL) AND 6-HR (AS IN OUR STUDY)
% DATA: 
% NORTH ATLANTIC, WESTERN NORTH PACIFIC VMAX [KTS, SYMM PART ONLY], MPI, VI
% DATE: 6/11/15
% AUTHOR: Emmi Yonekura

set(0,'defaultaxesfontsize',18);
set(0,'defaultaxesfontsize',18);

clear;
% OPTIONS:
basin = 'both'; %'NA'; %'WNP'; %'both';
cut = 'full'; %'short'; % 
vi_time = 't'; %'avg'; %'t+1'; % 
intens = 'V/MPI'; %'V'; %'V/MPI'; % 
cnt_thresh = 5; %count threshold for showing mean dv 

% LOAD DATA
% load fig_dat.mat *_vm *_relv %already divided by mpi
% dimensions: NA [137x613]; WNP [137x1037] IBTrACS WMO 1970-2010

% NEWLY (6/30/15) CALCULATED VI AND PI VARIABLES
load na_vivars7914_bt.mat  %vectors of trk pts
load wnp_vivars7914_bt.mat
na_VI = na_vi; wnp_VI = wnp_vi;
%Normalize_vi = (-na_NormalIntensity.^3+na_NormalIntensity).*sqrt(3);
% na_vi = Normalize_vi;

load wmo_na_1979_2014.mat fna
load wmo_wnp_1979_2014.mat fwnp
na_vm = fna.wind;  wnp_vm = fwnp.wind;

na_pi = nan(size(na_vm)); na_pi(isfinite(fna.lon)) = na_PI;
wnp_pi = nan(size(wnp_vm)); wnp_pi(isfinite(fwnp.lon)) = wnp_PI;

na_relv = nan(size(na_vm));
na_relv(isfinite(fna.lon)) = na_NormalIntensity; 

wnp_relv = nan(size(wnp_vm));
wnp_relv(isfinite(fwnp.lon)) = wnp_NormalIntensity; 


na_vi = nan(size(na_vm)); na_vi(isfinite(fna.lon)) = na_VI;
wnp_vi = nan(size(wnp_vm)); wnp_vi(isfinite(fwnp.lon)) = wnp_VI;

% remove negative v
na_vm(na_vm<0)=nan; wnp_vm(wnp_vm<0)=nan;
na_relv(na_vm<0) = nan; wnp_relv(wnp_vm<0)=nan;

% load vi_unfilt.mat na_vi wnp_vi

if strcmp(cut,'full')
    na_vi(na_vi<0)=nan; wnp_vi(wnp_vi<0)=nan;
elseif strcmp(cut,'short')
% further restrict ranges
na_vi(na_vi<0 | na_vi>3)=nan; 
wnp_vi(wnp_vi<0 | wnp_vi>3)=nan;
end

% SPECIFY WHICH BASIN BY ASSIGNING VARIABLES
if strcmp(basin,'NA')
%  - North Atlantic
    v = na_vm; 
    relv = na_relv;
    vi = na_vi; % at time t
elseif strcmp(basin,'WNP')
% %  - western North Pacific
    v = wnp_vm; 
    relv = wnp_relv;
    vi = wnp_vi; 
elseif strcmp(basin,'both')
    v = [na_vm wnp_vm];
    relv = [na_relv wnp_relv];
    vi = [na_vi wnp_vi];
end

% CALCULATE DV(6h) and DV(24h), allign matrices w/ time t
if strcmp(intens,'V')
    ntrk = size(v,2);
    dv6 = [v(2:end,:) - v(1:end-1,:); nan(1,ntrk)];
    dv24 = [v(5:end,:) - v(1:end-4,:); nan(4,ntrk)];
elseif strcmp(intens,'V/MPI')
    ntrk = size(relv,2);
    dv6 = [relv(2:end,:) - relv(1:end-1,:); nan(1,ntrk)];
    dv24 = [relv(5:end,:) - relv(1:end-4,:); nan(4,ntrk)];
end

% % USE VI AT WHICH TIME? 
if strcmp(vi_time,'t+1')
    vi_tplus1 = [vi(2:end,:); nan(1,ntrk)];
    vi = vi_tplus1;
elseif strcmp(vi_time, 'avg')
    vi_tplus1 = [vi(2:end,:); nan(1,ntrk)];
    vi = (vi + vi_tplus1)./2; 
end

% PARCEL INTO PHASE SPACE BINS
if strcmp(cut,'short')
%     vi_bins = [0:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:3];
    vi_bins = [0:0.01:0.99,1:0.5:3];
elseif strcmp(cut,'full')
    %vi_bins = [0:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:10, 20:10:100, 200];
%     vi_bins = [0:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:3];
    %vi_bins = [0:0.01:1];
    vi_bins = [0:0.01:0.99,1:0.5:3];
end

if strcmp(intens,'V/MPI')
    relv_bins = 0:0.05:1.5;
    nx = length(relv_bins);
    ny = length(vi_bins);

    cnt6 = zeros(nx,ny); %keep count of obs in bin
    cnt24 = zeros(nx,ny); %keep count of obs in bin
    mn_dv6 = zeros(nx,ny); %calc avg dv within bin
    mn_dv24 = zeros(nx,ny);

    n = numel(relv);
    for p = 1:n
        if isfinite(relv(p)*vi(p)*dv6(p)) %CHK NAN IN DV6!
            k = dsearchn(vi_bins', vi(p));
            j = dsearchn(relv_bins',relv(p));
            cnt6(j,k) = cnt6(j,k) + 1;
            mn_dv6(j,k) = (mn_dv6(j,k)*(cnt6(j,k)-1) + dv6(p))/ cnt6(j,k);
        end

        if isfinite(relv(p)*vi(p)*dv24(p)) %CHK NAN IN DV24!
            k = dsearchn(vi_bins', vi(p));
            j = dsearchn(relv_bins',relv(p));
            cnt24(j,k) = cnt24(j,k) + 1;
            mn_dv24(j,k) = (mn_dv24(j,k)*(cnt24(j,k)-1) + dv24(p))/ cnt24(j,k);
        end
    end
  
elseif strcmp(intens,'V')
    v_bins = 0:5:160;
    relv_bins = v_bins;
    nx = length(v_bins);
    ny = length(vi_bins);

    cnt6 = zeros(nx,ny); %keep count of obs in bin
    cnt24 = zeros(nx,ny); %keep count of obs in bin
    mn_dv6 = zeros(nx,ny); %calc avg dv within bin
    mn_dv24 = zeros(nx,ny);
    n = numel(v);
    for p = 1:n
        if isfinite(v(p)*vi(p)*dv6(p)) %CHK NAN IN DV6!
            k = dsearchn(vi_bins', vi(p));
            j = dsearchn(v_bins',v(p));
            cnt6(j,k) = cnt6(j,k) + 1;
            mn_dv6(j,k) = (mn_dv6(j,k)*(cnt6(j,k)-1) + dv6(p))/ cnt6(j,k);
        end
        

        if isfinite(v(p)*vi(p)*dv24(p)) %CHK NAN IN DV24!
            k = dsearchn(vi_bins', vi(p));
            j = dsearchn(v_bins',v(p));
            cnt24(j,k) = cnt24(j,k) + 1;
            mn_dv24(j,k) = (mn_dv24(j,k)*(cnt24(j,k)-1) + dv24(p))/ cnt24(j,k);
        end
    end
end

% REMOVE MEANS WITH N<cnt_thresh
low6 = logical(cnt6<cnt_thresh);
mn_dv6(low6) = nan;
cnt6(low6) = nan;
low24 = logical(cnt24<cnt_thresh);
mn_dv24(low24) = nan;
cnt24(low24) = nan;
% remove zeros
mn_dv6(mn_dv6==0) = nan;
mn_dv24(mn_dv24==0) = nan;


% MAKE PLOTS - use surf or patch

% figure; 
% surf(log10(vi_bins'),relv_bins',mn_dv24,'facealpha',1,'EdgeColor','none');
% colormap('bluewhitered(256)'); colorbar;
% 
% axis xy; axis([-2 0 0 1.2]) %1.4]) %
% xlabel(['logVI_{' vi_time '}']); ylabel(intens); 
% title([basin ' 24H OBS COUNT > ' num2str(cnt_thresh) ', VI at ' vi_time]);
% 
% ax = subplot(121);
% 
% surf(log10(vi_bins'),relv_bins',mn_dv24,'facealpha',1,'EdgeColor','none');
% colormap('bluewhitered(256)'); colorbar;
% 
% axis xy; axis([-2 0 0 1.2]) %1.4]) %
% xlabel(['logVI_{' vi_time '}']); ylabel(intens); 
% title([basin ' 24H OBS COUNT > ' num2str(cnt_thresh) ', VI at ' vi_time]);
% set(ax,'CLim',[-0.3 0.3]);

% ax = subplot(122); 
% 
% surf(log10(vi_bins'),relv_bins',mn_dv6,'facealpha',1,'EdgeColor','none');
% colormap('bluewhitered(256)'); colorbar;
% % set(ax,'CLim',[-0.3 0.3]);
% 
% % contour(log10(vi_bins'),relv_bins',cnt6,5:10:200); 
% axis xy; axis([-2 0 0 1.2]) %1.4]) %
% xlabel(['logVI_{' vi_time '}']); ylabel(intens); 
% title([basin ' 6H OBS COUNT >' num2str(cnt_thresh) ', VI at ' vi_time]);

% subplot(222); surf(log10(vi_bins),relv_bins,mn_dv24); colorbar; caxis ([-0.3 0.3]); 
% axis xy;view([0 90]);axis([-2 0 0 1.2]) %1.4]) %
% xlabel(['logVI_{' vi_time '}']); ylabel(intens); title([basin ' 24H D(' intens ')']);
% 
% subplot(224); surf(log10(vi_bins),relv_bins,mn_dv6); colorbar; caxis ([-0.3 0.3]); 
% axis xy;view([0 90]);axis([-2 0 0 1.2]) %1.4]) %
% xlabel(['logVI_{' vi_time '}']); ylabel(intens); title([basin ' 6H D(' intens ')']);


figure; 
subplot(121); contour(vi_bins,relv_bins',cnt24,[20,50,90,120,150,300],'k','LineWidth',2,'ShowText','on'); 
hold on;
surf(vi_bins,relv_bins,mn_dv24,'facealpha',0.7,'EdgeColor','none'); colorbar; caxis ([-0.22 0.22]); 
colormap('bluewhitered(256)');
set(gca,'xscale','log');
% set(gca,'XTick',vi_bins); 
% Xticklabel = cell(1,length(vi_bins));
% for i = 1:length(vi_bins)
%    if(vi_bins(i)==0.001||vi_bins(i)==0.01||vi_bins(i)==0.1||vi_bins(i)==1)
%        Xticklabel{i} = num2str(vi_bins(i));
%    else
%        Xticklabel{i} = [];
%    end
% end      
set(gca,'XTick',[0 0.001 0.01 0.1 1])
% set(gca,'XTicklabel',Xticklabel); 
axis([0.0006 3 0 1.2]) %1.4]) %
xlabel(['VI_{' vi_time '}']); ylabel(intens); 
title([basin ' 24H OBS COUNT > ' num2str(cnt_thresh) ', VI at ' vi_time]);

subplot(122); contour(vi_bins,relv_bins',cnt6,[20,50,90,120,150,300],'k','LineWidth',2,'ShowText','on'); 
hold on;
surf(vi_bins,relv_bins,mn_dv6,'facealpha',0.7,'EdgeColor','none'); colorbar; caxis ([-0.06 0.06]); 
set(gca,'xscale','log');
set(gca,'XTick',[0 0.001 0.01 0.1 1]);
axis([0.0006 3 0 1.2]) 
colorbar;xlabel(['VI_{' vi_time '}']); ylabel(intens); 
title([basin ' 6H OBS COUNT >' num2str(cnt_thresh) ', VI at ' vi_time]);


% subplot(222); surf(vi_bins,relv_bins,mn_dv24,'facealpha',0.1,'EdgeColor','none');alpha 0.2; colorbar; caxis ([-0.25 0.25]); 
% % hold on;
% % contour(vi_bins,relv_bins',cnt6,10:20:200,'k'); 
% colormap('bluewhitered(256)');
% set(gca,'xscale','log');
% set(gca,'XTick',[0.01 0.1 1]);
% axis xy;view([0 90]);axis([0 1 0 1.2]) %1.4]) %
% xlabel(['VI_{' vi_time '}']); ylabel(intens); title([basin ' 24H D(' intens ')']);
% 
% subplot(224); surf(vi_bins,relv_bins,mn_dv6); colorbar; caxis ([-0.06 0.06]); 
% set(gca,'xscale','log');
% set(gca,'XTick',[0.01 0.1 1]);
% axis xy;view([0 90]);axis([0 1 0 1.2]) %1.4]) %
% xlabel(['VI_{' vi_time '}']); ylabel(intens); title([basin ' 6H D(' intens ')']);
% 
% 
% return;
% NO MORE LOGVI
% figure; 
% subplot(221); contour(vi_bins',relv_bins',cnt24,5:10:200); 
% colorbar; axis xy; axis([0 0.2 0 1.55]) %1.2]) %
% xlabel(['VI_{' vi_time '}']); ylabel(intens); 
% title([basin ' 24H OBS COUNT > ' num2str(cnt_thresh) ', VI at ' vi_time]);
% 
% subplot(223); contour(vi_bins',relv_bins',cnt6,5:10:200); 
% colorbar; axis xy; axis([0 0.2 0 1.55]) % 1.2]) %
% xlabel(['VI_{' vi_time '}']); ylabel(intens); 
% title([basin ' 6H OBS COUNT >' num2str(cnt_thresh) ', VI at ' vi_time]);
% 
% subplot(222); surf(vi_bins,relv_bins,mn_dv24); colorbar; caxis ([-0.3 0.3]); 
% axis xy;view([0 90]);axis([0 0.2 0 1.55]) % 1.2]) %
% xlabel(['VI_{' vi_time '}']); ylabel(intens); title([basin ' 24H D(' intens ')']);
% 
% subplot(224); surf(vi_bins,relv_bins,mn_dv6); colorbar; caxis ([-0.3 0.3]); 
% axis xy;view([0 90]);axis([0 0.2 0 1.55]) %1.2]) %
% xlabel(['VI_{' vi_time '}']); ylabel(intens); title([basin ' 6H D(' intens ')']);

return;
% % SO THAT FAILED. LET'S ALSO PLOT HISTOGRAMS OF VI FOR POS/NEG CHANGE
% 
% % PARCEL INTO PHASE SPACE BINS
% if strcmp(cut,'short')
%     vi_bins = [0:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:3];
% elseif strcmp(cut,'full')
%     vi_bins = [0:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:10, 20:10:100, 200];
% end
% ny = length(vi_bins);
% 
% cnt6_pos = zeros(ny,1); %keep count of obs in bin
% cnt6_neg = zeros(ny,1); %keep count of obs in bin
% cnt24_pos = zeros(ny,1); %keep count of obs in bin
% cnt24_neg = zeros(ny,1); %keep count of obs in bin
% 
% n = numel(relv);
% for p = 1:n
%     if isfinite(vi(p)*dv6(p)) %CHK NAN IN DV6!
%         k = dsearchn(vi_bins', vi(p));
%         if dv6(p) >=0
%             cnt6_pos(k) = cnt6_pos(k) + 1;
%         else
%             cnt6_neg(k) = cnt6_neg(k) + 1;
%         end
%     end
%     
%     if isfinite(vi(p)*dv24(p)) %CHK NAN IN DV24!
%         k = dsearchn(vi_bins', vi(p));
%         if dv24(p) >=0
%             cnt24_pos(k) = cnt24_pos(k) + 1;
%         else
%             cnt24_neg(k) = cnt24_neg(k) + 1;
%         end
% 
%     end
% end
% figure;
% subplot(211)
% bar([cnt6_neg cnt6_pos],'grouped')
% legend('NEG','POS')
% xlabel('VI')
% set(gca, 'XTick',1:4:48,'XTickLabel',vi_bins(1:4:48))
% title(['Histogram 6hr d(' intens ')'])
% 
% subplot(212)
% bar([cnt24_neg cnt24_pos],'grouped')
% legend('NEG','POS')
% xlabel('VI')
% set(gca, 'XTick',1:4:48,'XTickLabel',vi_bins(1:4:48))
% title(['Histogram 24hr d(' intens ')'])


% Multiple linear regression of just these variables
stats = regstats(y,X,'linear'); 
% performs a multilinear regression of the responses in y on 
% the predictors in X. X is an n-by-p matrix of p predictors 
% at each of n observations. y is an n-by-1 vector of observed 
% responses.
s=regstats(reshape(dv6,n,1),[reshape(relv,n,1) reshape(vi,n,1)],'linear');
% using both basins, run out of memory
% trying NA...

