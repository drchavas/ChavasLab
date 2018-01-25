%sounding_drc.m -- extract mean environmental profile from ax simulation
%
% Syntax:  [sounding_all,pressure_all,zz00_RCE,th00_RCE,qv00_RCE,u00,v00,...
%       p_sfc,qv_sfc_RCE,th_sfc_RCE,pi00_RCE] = sounding_drc(dir_in,...
%       simulation_name,t0,tf,x0,xf,moist,override_instab,...
%       override_nonequil,override_dampisothermal)
%
% Inputs:
%   dir_in - path to simulation subdirectory is located
%   simulation_name - name of simulation subdirectory
%   t0 [day] - time of start of desired averaging period
%   tf [day] - time of end of desired averaging period
%   x0 [-] - inner radial gridpoint number
%   xf [-] - outer radial gridpoint number
%   moist [0/1] - 1 = moist, 0 = dry sounding
%   override_instab [0/1] - override error detecting large instabillity
%   override_nonequil [0/1] - override error detecting non-equilibration
%   override_dampisothermal [0/1] - override error detecting non-isothermal
%       damping layer
%
% Outputs:
%   sounding_all - matrix of values ready for input into sounding file
%   pressure_all - matrix of values ready for input into sounding pressure file
%   zz00_RCE [m] - RCE heights (does not include top-layer buffer in sounding_all)
%   th00_RCE [K] - RCE potential temp (does not include top-layer buffer in sounding_all)
%   qv00_RCE [kg/kg] - RCE water vapor mixing ratio (does not include top-layer buffer in sounding_all)
%   u00 [m/s] - zonal wind, SET TO ZERO (does not include top-layer buffer in sounding_all)
%   v00 [m/s] - meridional wind, SET TO ZERO (does not include top-layer buffer in sounding_all)
%   p_sfc [Pa] - surface pressure; taken from initial sounding
%   qv_sfc_RCE [kg/kg] - near-surface water vapor mixing ratio; set to same
%       value as 1st model level
%   th_sfc_RCE [K] - near-surface potential temp; set to same
%       value as 1st model level
%   pi00_RCE [-] - non-dim pressure on same levels as zz00_RCE
%
% Example:
%
% Other m-files required: nc_extract_axisym, snd_extract
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 23 Feb 2015; Last revision:

%------------- BEGIN CODE --------------


function [sounding_all,pressure_all,zz00_RCE,th00_RCE,qv00_RCE,u00,v00,p_sfc,qv_sfc_RCE,th_sfc_RCE,pi00_RCE] = sounding_drc(dir_in,simulation_name,t0,tf,x0,xf,moist,override_instab,override_nonequil,override_dampisothermal)

if(nargin<7)
    moist = 1;  %default is moist
end
if(nargin<8)
    override_instab = 0;    %default is to give error if there is a large instability (indicates some sort of problem)
end
if(nargin<9)
    override_nonequil = 1;    %default is to ignore error for large temporal changes in vertical profile during period
end
if(nargin<10)
    override_dampisothermal = 0;    %default is to give error if damping layer is not isothermal
end

%% Domain bounds
y0=0;
yf=10000;
z0=0;
zf=10000; %large number if want all levels

%% Constants (values taken from CM1 model)
c_CM1 = constants_CM1(); %c_CM1: [g rd cp rv p00 xlv cpv]

%g=c_CM1(1); %[m/s2]
Rd=c_CM1(2);  %[J/kg/K]
%Cpd=c_CM1(3); %[J/kg/K]; spec heat of dry air
Rv=c_CM1(4);   %[J/K/kg]
%p0 = c_CM1(5); %[Pa]
%Lv=c_CM1(6);   %[J/kg]
%Cpv=c_CM1(7); %[J/kg/K]; spec heat of water

eps=Rd/Rv;

assert(abs(eps-.622) < .001,'CM1_constants has an error!')

%% Extract environmental mean sounding
    
%%EXTRACT RUN_TYPE AND simulation_name NAME
dir_simulation = sprintf('%s/%s',dir_in,simulation_name);

%% Extract T_sst from file name
i_sst=strfind(simulation_name,'SST');
if(isempty(i_sst))
    sst_str='300.00';
else
    i_sstK=strfind(simulation_name,'K');
    i_sstK=i_sstK(find(i_sstK>i_sst,1));
    sst_str=simulation_name(i_sst+3:i_sstK-1);
end
SST = str2num(sst_str);

%% Extract data from netcdf file
run_type_str='axisym';

%%EXTRACT TIMESTEP SIZE FROM ARBITRARY FILE
var = 'thpert'; %doesn't matter, just need any variable
clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
numfiles=length(dir(sprintf('%s/cm1out_*.nc',dir_simulation)));
[~,xmin_sub,xmax_sub,~,~,zmin_sub,zmax_sub,dx,~,dz,~,~,nz_sub,~,~,~,~,~,time,~] = ...
    nc_extract_axisym(dir_in,simulation_name,'cm1out_000002.nc',var,x0,xf,y0,yf,z0,zf);

zz00_RCE = zmin_sub:dz:zmax_sub;    %heights [m]

dt = time; %[s]; time of cm1out_0001.nc is defined as zero
%i_tf = numfiles-1;    %last file of simulation
%tf = (i_tf-1)*dt/(24*60*60);
%i_t0 = max(0,round(tf-Tmean_RCE*24*60*60/dt)+1); %timestep corresponding to start time for averaging
%t0 = (i_t0-1)*dt/(24*60*60);
i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
i_tf = min(numfiles-1,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging

nt = i_tf-i_t0+1;

assert(nt>0,sprintf('Those timesteps dont work for %s; max available time is %3.3f day',simulation_name,i_tf*dt/86400))

%% EXTRACT DATA FROM input_sounding
snd_file = 'input_sounding';
%load initial sounding variables
[~,pp00,th00,qv00,u00,v00,~,~,~,~,~,~,~,p_sfc,~,~] = ...
    snd_extract(dir_simulation,snd_file,dz,nz_sub);


%% LOOP OVER FILES TO OBTAIN TIME SERIES DATA
%%initialize the output vectors arbitrarily large
clear data_hmean_qv data_hmean_th data_hmean_p data_hmean_pi
data_hmean_qv=zeros(1000,1);
data_hmean_th=zeros(1000,1);
data_hmean_p=zeros(1000,1);
data_hmean_pi=zeros(1000,1);
vert_prof_temp_dth=zeros(1000,1);

for ii=1:i_tf-i_t0+1

    t_file=i_t0+ii-1;

    %%NC FILENAME (default)
    if(t_file<10)
        nc_file = sprintf('cm1out_00000%i.nc',t_file);
    elseif(t_file<100)
        nc_file = sprintf('cm1out_0000%i.nc',t_file);
    elseif(t_file<1000)
        nc_file = sprintf('cm1out_000%i.nc',t_file);
    elseif(t_file<10000)
        nc_file = sprintf('cm1out_00%i.nc',t_file);
    elseif(t_file<100000)
        nc_file = sprintf('cm1out_0%i.nc',t_file);
    else
        nc_file = sprintf('cm1out_%i.nc',t_file);
    end

    %%EXTRACT qv DATA: note that qvpert is in [kg/kg]
    if(moist == 1)
        var = 'qvpert';
        clear data
        [data,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = nc_extract_axisym(dir_in,simulation_name,nc_file,var,x0,xf,y0,yf,z0,zf);

        %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
%         v_def_qv = v_def;
%         v_units_qv = v_units;
        data_hmean_qv=data_hmean_qv(1:nz_sub);   %NOTE qvpert IS IN [kg/kg]!!!
        xvals = xmin_sub:dx:xmax_sub;
        data=squeeze(data);
        vert_prof_temp_qv = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius

        data_hmean_qv=data_hmean_qv+vert_prof_temp_qv/(i_tf-i_t0+1); %[kg/kg]
    end

    %%EXTRACT theta DATA
    var = 'thpert';
    clear data
    [data,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = nc_extract_axisym(dir_in,simulation_name,nc_file,var,x0,xf,y0,yf,z0,zf);

    %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
%     v_def_th = v_def;
%     v_units_th = v_units;
    data_hmean_th=data_hmean_th(1:nz_sub);
    vert_prof_temp_dth=vert_prof_temp_dth(1:nz_sub);
    if(ii>=2)
        vert_prof_temp_th_save = vert_prof_temp_th;
    end

    xvals = xmin_sub:dx:xmax_sub;
    data=squeeze(data);
    vert_prof_temp_th = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius

    data_hmean_th=data_hmean_th+vert_prof_temp_th/(i_tf-i_t0+1); %[g/kg]

    if(ii>=2)
        vert_prof_temp_dth = vert_prof_temp_dth + (vert_prof_temp_th - vert_prof_temp_th_save); %integrated change over time
    end

    %%EXTRACT pressure DATA
    var_temp = 'prspert';
    %        var_temp = 'vinterp'
    clear data
    [data,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = nc_extract_axisym(dir_in,simulation_name,nc_file,var_temp,x0,xf,y0,yf,z0,zf);

    %        vv = data;

    %        var_temp = 'uinterp'
    %        clear data
    %        [data,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = nc_extract_axisym(dir_in,simulation_name,nc_file,var_temp,x0,xf,y0,yf,z0,zf);

    %        uu = data;

    %        data = sqrt(uu.^2 + vv.^2);

    %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
%     v_def_p = v_def;
%     v_units_p = v_units;
    data_hmean_p=data_hmean_p(1:nz_sub);
    xvals = xmin_sub:dx:xmax_sub;
    data=squeeze(data);
    vert_prof_temp_p = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius
    data_hmean_p=data_hmean_p+vert_prof_temp_p/(i_tf-i_t0+1); %[Pa]
    %        wspd = data_hmean_p;

    clear data

    %%EXTRACT pi DATA
    var_temp = 'pi';
    clear data
    [data,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = nc_extract_axisym(dir_in,simulation_name,nc_file,var_temp,x0,xf,y0,yf,z0,zf);

    %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
%     v_def_pi = v_def;
%     v_units_pi = v_units;
    data_hmean_pi=data_hmean_pi(1:nz_sub);
    xvals = xmin_sub:dx:xmax_sub;
    data=squeeze(data);
    vert_prof_temp_pi = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius
    data_hmean_pi=data_hmean_pi+vert_prof_temp_pi/(i_tf-i_t0+1); %nondim

    clear data
    %}

end

vec_length = min(nz_sub,length(zz00_RCE)-z0);


qv00_RCE = qv00(z0+1:z0+vec_length)' + [data_hmean_qv(z0+1:z0+vec_length)];   %qv' [kg/kg]
th00_RCE = th00(z0+1:z0+vec_length)' + [data_hmean_th(z0+1:z0+vec_length)];   %th' [K]
p00_RCE = pp00(z0+1:z0+vec_length)' + [data_hmean_p(z0+1:z0+vec_length)];   %p' [Pa]
pi00_RCE = [data_hmean_pi(z0+1:z0+vec_length)];   %FULL pi [-]

%alternate method
%    p200_RCE = 100000*(data_hmean_pi(z0+1:z0+vec_length)).^(1004/287);

%%Fix any vertical instabilities in th00_RCE (just set value to mean of levels it lies between)
%%MOIST ONLY! For dry case, this needs to be updated since dth/dz<0 in lowest several model levels
%TEST    th00_RCE(10)=th00_RCE(9)-.001;
%TEST    th00_RCE(9)=th00_RCE(8)-.01;
if(moist==1)

    RCE_instab_remove = 0;

    thv00_RCE=th00_RCE.*(1+qv00_RCE/eps)./(1+qv00_RCE);

    instab_check = thv00_RCE(2:end)-thv00_RCE(1:end-1);
    max_instab = min(instab_check);
    z_instab = find(instab_check == max_instab);

    if(max_instab<0)    %there is an instability to remove

        if(override_instab~=1)
            assert(z_instab>=2 || abs(max_instab)<.3,'WARNING: THERE IS A LARGE NEAR-SURFACE INSTABILITY IN MOIST RCE PROFILE')
            assert(abs(max_instab)<.3,'WARNING: THERE IS A LARGE INSTABILITY IN MOIST RCE PROFILE')
        end

        RCE_instab_remove = 1;   %the profile will be adjusted

        %%Iterate upwards through the RCE profile and apply mass-weighted averaging ("mixing") over
        %%each statically-unstable layer and repeat until the entire profile is stable
        thv00_RCE_adj = thv00_RCE;  %adjusted profile begins as true RCE profile
        dthv = thv00_RCE_adj(2:end)-thv00_RCE_adj(1:end-1);  %change in theta with height

        while(sum(dthv<-.0001)>0)   %there are unstable layers
            for j = 1:length(thv00_RCE_adj)-1
                dthv_j = thv00_RCE_adj(j+1)-thv00_RCE_adj(j);
                if(dthv_j<0)   %unstable
                    thv_ave = (p00_RCE(j)*thv00_RCE_adj(j) + p00_RCE(j+1)*thv00_RCE_adj(j+1))/(p00_RCE(j)+p00_RCE(j+1));
                    thv00_RCE_adj(j) = thv_ave;
                    thv00_RCE_adj(j+1) = thv_ave;
                end
            end

            %%NOTE: I do NOT mix qv as well, as this is not a true
            %%instability but instead a result of numerical noise, unlike
            %%in the dry case (though in that case, qv=0 anyways)

            %update change in theta with height
            dthv = thv00_RCE_adj(2:end)-thv00_RCE_adj(1:end-1);  %change in theta with height

        end

        %%Adjust th00_RCE according to adjustmnet in thv00_RCE assuming qv00_RCE fixed
        th00_RCE_adj=thv00_RCE_adj.*(1+qv00_RCE)./(1+qv00_RCE/eps);
        th00_RCE_temp = th00_RCE;
        th00_RCE = th00_RCE_adj;    %only need the final adjusted profile

    end

end
%}

%%Surface theta and qv
%    SST_C = SST-273.15;
%    es_SST = 6.112*exp((17.67*SST_C)./(SST_C+243.5))*100;  %[Pa]
%    qvs_SST = eps*(es_SST./(p_sfc-es_SST));  %[kg/kg]

qv_sfc_RCE = qv00_RCE(1);
th_sfc_RCE = th00_RCE(1);

%%Save input soundings for plotting
%{
qv_sfc_init = qv_sfc;
th_sfc_init = th_sfc;
zz00_init = zz00(z0+1:z0+vec_length)'; %heights [m]
qv00_init = qv00(z0+1:z0+vec_length)';   %[kg/kg]
th00_init = th00(z0+1:z0+vec_length)';   %[K]
%}



%%sounding header: %1) p_sfc (mb); 2) th_sfc (K); 3) qv_sfc (g/kg)
%%sounding columns: 1) zz00 (m); 2) th00_RCE (K); 3) qv00_RCE (g/kg); 4) u00 (m/s); 5) v00 (m/s)

%%Ensure that qv00_RCE > .00001 g/kg
qv00_RCE(qv00_RCE<.00001/1000)=.00001/1000;

%%need extra layers at top above model top
zz00_ext = zz00_RCE(end)+(zz00_RCE(end)-zz00_RCE(end-1))*[1];
th00_ext = th00_RCE(end)+(th00_RCE(end)-th00_RCE(end-1))*[1];
qv00_ext = qv00_RCE(end)*[1];
u00_ext = u00(end)*[1];
v00_ext = v00(end)*[1];

sounding_all=[1000*[zz00_RCE';zz00_ext'] [th00_RCE;th00_ext'] 1000*[qv00_RCE;qv00_ext'] [u00(1:vec_length)';u00_ext'] [v00(1:vec_length)';v00_ext']];
sounding_all=double(sounding_all);

pressure_all=[p00_RCE];
pressure_all=double(pressure_all);

%%Make adjusted dry RCE sounding + perturbation (dry RCE is statically unstable)
if(moist == 0)

    %%Iterate upwards through the RCE profile and apply mass-weighted averaging ("mixing") over
    %%each statically-unstable layer and repeat until the entire profile is stable
    th00_RCE_adj = th00_RCE;  %adjusted profile begins as true RCE profile
    dth = th00_RCE_adj(2:end)-th00_RCE_adj(1:end-1);  %change in theta with height

    sprintf('Producing adjusted dry RCE sounding that is statically stable; also perturbation')
    while(sum(dth<0)>0)   %there are unstable layers
        for j = 1:length(th00_RCE_adj)-1
            dth_j = th00_RCE_adj(j+1)-th00_RCE_adj(j);
            if(dth_j<0)   %unstable
                th_ave = (p00_RCE(j)*th00_RCE_adj(j) + p00_RCE(j+1)*th00_RCE_adj(j+1))/(p00_RCE(j)+p00_RCE(j+1));
                th00_RCE_adj(j) = th_ave;
                th00_RCE_adj(j+1) = th_ave;
            end
        end

        %update change in theta with height
        dth = th00_RCE_adj(2:end)-th00_RCE_adj(1:end-1);  %change in theta with height

    end

    %magnitude of the adjustment from the true RCE
    th00_RCE_pert = th00_RCE-th00_RCE_adj;

    %%Save soundings for the adjusted profile and the perturbation

    %ADJUSTED PROFILE
    sounding_all=[1000*[zz00_RCE;zz00_ext'] [th00_RCE_adj;th00_ext'] 1000*[qv00_RCE;qv00_ext'] [u00(1:vec_length)';u00_ext'] [v00(1:vec_length)';v00_ext']];
    sounding_all=double(sounding_all);

end

%% TESTING FOR INSTABILITY/EQUILIBRATION

vert_prof_dth_max = max(vert_prof_temp_dth);
z_vert_prof_dth_max = find(vert_prof_temp_dth == vert_prof_dth_max);
if(override_nonequil~=1)
    assert(vert_prof_dth_max < dth_max_equil,'WARNING: NOT EQUILIBRATED, max(dth_equil) = %5.5f at level %i',vert_prof_dth_max,z_vert_prof_dth_max)
end

%% Calculate absolute temperature RCE profile
T00_RCE = th00_RCE.*pi00_RCE;

T_RCE_test = T00_RCE;
if(override_dampisothermal~=1)
    assert(std(T_RCE_test(end-5:end))<.1,'WARNING: DAMPING LAYER MAY NOT BE ISOTHERMAL!')
end

if(RCE_instab_remove == 1)  %an instability was removed
    sprintf('%s',simulation_name)
    sprintf('A vertical static instability was removed from the RCE sounding')
    sprintf('Maximum instability = %5.4f K (th_v) at level %i',max_instab,z_instab)
else
    sprintf('%s',simulation_name)
    sprintf('There were no instabilities in the RCE sounding')
end

fclose('all');

%%Adjust height units to [m]
zz00_RCE = 1000*zz00_RCE;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: PLOT %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
clf     %clear figures without closing them

%%PLOT VERTICAL PROFILES
%PLOT DATA
%vec_length=min(sum(~isnan(qv00_RCE(:,rr))),length(zz00));

figure(1)
subplot(1,2,1)
hold on
%plot(1000*[qv_sfc_RCE;qv00_RCE],[0;zz00_RCE/1000],pl_clrs,'LineWidth',1)
hpl1 = plot(1000*[qv00_RCE],[zz00_RCE/1000],pl_clrs,'LineWidth',1);

set(gca,'fontweight','bold','fontsize',11)

figure(1)
subplot(1,2,2)
hold on
plot([T00_RCE],[zz00_RCE/1000],pl_clrs,'LineWidth',1,'LineStyle','--');
hpl2 = plot([th00_RCE],[zz00_RCE/1000],pl_clrs,'LineWidth',1);

set(gca,'fontweight','bold','fontsize',11)

%%Plot initial soundings
%{
figure(1)
subplot(1,2,1)
hold on
plot(1000*[qv_sfc_init;qv00_init],[0;zz00_init/1000],pl_clrs{2*rr},'LineWidth',1)

figure(1)
subplot(1,2,2)
hold on
plot([th_sfc_init;th00_init],[0;zz00_init/1000],pl_clrs{2*rr},'LineWidth',1)
%}



figure(1)
if(moist == 1)
subplot(1,2,1)
set(gca,'fontweight','bold','fontsize',11)
axis([0 1.01*1000*max(cellfun(@(x) max(x(:)),qv00_RCE)) 0 zz00(end)/1000])
input_xlabel=sprintf('q_v [g/kg]',v_def_qv);
xlabel(input_xlabel);
input_ylabel=sprintf('Height AGL [km]');
ylabel(input_ylabel);
input_legend=strrep([simulation_names_plot],'_','\_');
legend(hpl1,input_legend)
input_title=strrep(sprintf('water vapor mixing ratio [g/kg]'),'_','\_');
title(input_title)

figure(1)
subplot(1,2,2)
end
set(gca,'fontweight','bold','fontsize',11)
axis([.9*min(min(cellfun(@(x) min(x(:)),T00_RCE))) max(max(cellfun(@(x) max(x(:)),th00_RCE))) 0 zz00(end)/1000])
input_xlabel=sprintf('theta / T [%s]',v_units_th);
xlabel(input_xlabel);
input_ylabel=sprintf('Height AGL [km]');
ylabel(input_ylabel);
input_legend=strrep([simulation_names_plot],'_','\_');
legend(hpl2,input_legend)
input_title=strrep(sprintf('theta / T (dash) [%s]',v_units_th),'_','\_');
title(input_title)

input_title=strrep(sprintf('vert profs: horiz-mean days %i-%i; i0-if = %i-%i (dr = %2.0f km)',t0,tf,x0,xf,dx),'_','\_');
annotation('textbox','String',input_title,'Position',[.1 .5 .9 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')

%    title({input_title1,input_title2})
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%------------- END OF CODE --------------