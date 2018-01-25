%TC_wind_fourier_decomp_test.m

clear
clc
close('all')

load wind_data_ex.mat
    %vars: rr, th_grid, VV_TC

r1 = 100*1000;  %[m]
r2 = 110*1000;  %[m]

ii_good = find(rr>=r1 & rr<=r2);

rr_good = rr(ii_good);
th_good = th_grid(ii_good);
VV_good = VV_TC(ii_good);
size(VV_good)

[th_good,i_sort] = sort(th_good);
rr_good = rr_good(i_sort);
VV_good = VV_good(i_sort);

figure(1)
polar(th_good,rr_good/1000,'rx')
figure(2)
plot(th_good,VV_good)
xlabel('azimuth [deg]')
ylabel('wind speed [ms-1]')
axis([0 360 0 1.1*max(VV_good)])

%% Perform fourier decomposition
p = abs(fft(VV_good))/(N/2);    %% absolute value of the fft
p = p(1:N/2).^2;    %% take the positve frequency half, only
freq = [0:N/2-1]/T; %[Hz]; frequency in real time

semilogy(freq,p);   %plot power spectrum
axis([0 30 0 1]);