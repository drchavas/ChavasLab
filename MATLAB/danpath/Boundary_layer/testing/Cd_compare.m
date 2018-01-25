%Cd_compare.m

clear
clc
close all
set(0,'DefaultFigureVisible', 'on');

addpath(genpath('~/Dropbox/Research/MATLAB/'));

figure(4)
clf(4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%From Donelan2004.m -- email from Mark Donelan 2014-06-17
% Donelan et al 2004 GRL: CD vs U10; 
% Data lifted off Figure 2 by Alex Soloviev.
%Squares
VD1 =  ([263 251 225 223 198 194 174 166 152 139 116 130 111 103   97  92  85  74  69  64   63  60  58]-46)*60/(560-46);
CD1 = -([289 296 309 310 324 326 338 342 350 358 370 362 371 376  379 380 382 382  379 372 367 358 341]-431)*5e-3/(431-26);
%Asterisks
VD2 = ([180 226 276 334 412 400 376]-46)*60/(560-46);
CD2 = -([351 324 278 235 246 246 237]-431)*5e-3/(431-26);
%Circles
VD3 = ([209 502 476 384 254 307 429 356]-46)*60/(560-46);
CD3 = -1.0*([340 244 231 255 305 265 241 241]-431)*5e-3/(431-26);
%Diamonds
VD4 = ([396 433 315 267 183 148 354 220 84]-46)*60/(560-46);
CD4 = -([237 237 248 274 324 341 243 307 381]-431)*5e-3/(431-26);

figure(4);clf;plot(VD1,CD1,'s');
figure(4);hold on;plot(VD2,CD2,'*');
figure(4);hold on;plot(VD3,CD3,'o');
figure(4);hold on;plot(VD4,CD4,'d');
axis([0 60 0 .0045])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



V_in = 1:50;
[Cd_D04] = Cd_Donelan04(V_in);
[Cd_S08] = Cd_Smith2008Black07(V_in);
[Cd_S14] = Cd_Smith2014CBLAST07(V_in);

hold on
hpl(1) = plot(V_in,Cd_D04,'r')
hold on
hpl(2) = plot(V_in,Cd_S08,'g')
hpl(3) = plot(V_in,Cd_S14,'c')
legend(hpl,{'D04','S08','S14'},'Location','SouthEast')
title('Blue symbols = Donelan et al. 2004 data')

%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('Cd_compare.pdf')
saveas(gcf,plot_filename,'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
