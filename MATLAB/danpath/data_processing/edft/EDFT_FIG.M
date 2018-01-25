% This program calculate and plot output of spectral estimation methods discussed in [4]:
%	- DFT,
%	- Extended DFT (EDFT),
%	- High-Resolution DFT (HRDFT) [2],
%	- Iterative algorithms based on the Capon filter approach and GWLS solution
% EDFT spectral estimate is tested on the three 64-point sequences: 
% x_uni - uniform complex-value test sequence used in [3];
% x_non - nonuniform complex-value test sequence used in [3];
% x_mak - Well known Marple&Kay 64-point real sample sequence [1].
% The spectral plots are not intended to demonstrate relative performance of approaches, but only
% to point out features of each technique. 
%
% References: 
% [1] - S.M.Kay, S.L.Marple. Spectrum analysis - a modern perspective. Proceedings IEEE 69, 1981.
% [2] - Mauricio D. Sacchi, Tadeusz J. Ulrych, Colin J. Walker. Interpolation and Extrapolation Using a 
%       High-Resolution Discrete Fourier Transform. IEEE TRANS. ON SIGNAL PROC., VOL.46, NO.1, JANUARY 1998.
% [3] - Vilnis Liepins. High-resolution spectral analysis by using basis function adaptation approach.
%       Doctoral Thesis for Scientific Degree of Dr.Sc.Comp., University of Latvia, 1997. 
%       Abstract available on http://www.opengrey.eu/item/display/10068/330816
% [4] - Vilnis Liepins. Extended Fourier analysis of signals. CoRR abs/1303.2033 (2013).
%       Available online at http://arxiv.org/abs/1303.2033 
%
% E-mail:	vilnislp@gmail.com

	clear
	it=10;				% Number of iterations
	N=1000;				% Number of frequencies
	fn=[-ceil((N-1)/2):floor((N-1)/2)]/N;	% Uniform frequency set [-0.5 ... 0.5[ (sampling frequency=1)
	fne=[-ceil((2*N-1)/2):floor((2*N-1)/2)]/N;	% Extended uniform frequency set [-1 ... 1[
	K=64;				% Length of test sequences 1 and 2	
        	t_u=0:K-1;              			% Uniform time set (sampling period=1)
	rand('seed',137);			% uncomment to generate always the same input sequences
        	t_n=t_u+rand(1,K)*0.8;  		% Nonuniform time set (mean sampling period=1)
	Eu=exp(-i*2*pi*t_u.'*fn);		% Exponents matrix (Uniform time), dim (KxN)
	En=exp(-i*2*pi*t_n.'*fn);		% Exponents matrix (Nonuniform time), dim (KxN)
	Eue=exp(-i*2*pi*t_u.'*fne);		% Extended Exponents matrix (Uniform time), dim (Kx2N)
	Ene=exp(-i*2*pi*t_n.'*fne);		% Extended Exponents matrix (Nonuniform time), dim (Kx2N)

% Uniform 64-point complex-value sequence: x_uni
	PHT=2*pi*(rand-0.5);			% Random initial phase for the complex exponent
	x_s1=exp(i*(2*pi*0.35*t_u+PHT));	% Complex exponent at normalized frequency 0.35
	x_i=20*sin(pi*(t_u-K/2)*0.25)./(pi*(t_u-K/2)+eps).*exp(i*pi*0.25*(t_u-K/2));	% Rectangular pulse [0 ... 0.25]
	x_i(K/2+1)=5;
	f_x=-0.5:1/1024:-0.25;			% Set 257 frequencies for simulation of a band-limited noise
	rand_ph=2*rand(257,1)*ones(1,K);			
	x_n=sum(exp(i*2*pi*(f_x'*t_u+rand_ph)))/sqrt(K); % Band-limited noise [-0.5 ...-0.25]
	x_uni=x_s1+x_i+x_n;			% x_uni - uniform composite sequence
	xmax=ceil(max(abs(x_uni)));
	lim=2^9+0.5;
	sigm_x=x_uni*x_uni'/K;					
	x_uni=xmax*round(x_uni/xmax*lim)/lim; 	% Simulate 10-bit Analog-to-Digital Converter
	sigm_n=sqrt((xmax/lim)^2/6/K); 					
	WT=[20*ones(1,N/4)/K sigm_n*ones(1,N/4) 20*ones(1,N/4)/K sigm_n*ones(1,N/4)];
	WT(find(fn==0.35))=1;			% WT - true spectrum of signals x_uni and x_non
	WE=[sigm_n*ones(1,N/2) WT sigm_n*ones(1,N/2)]; % WE - extended true spectrum of signals x_uni and x_non
	sigm_p=(xmax/lim)^2*6/N;
	WP=[(20*ones(1,N/4)).^2/N sigm_p*ones(1,N/4) (20*ones(1,N/4)).^2/N sigm_p*ones(1,N/4)];
	WP(find(fn==0.35))=N;		% WP - true power spectal density of signals x_uni and x_non

% Nonuniform 64-point complex-value sequence: x_non
	x_ns=exp(i*(2*pi*0.35*t_n+PHT));	% complex exponent at normalized frequency 0.35
	x_ni=20*sin(pi*(t_n-K/2)*0.25)./(pi*(t_n-K/2)).*exp(i*pi*0.25*(t_n-K/2));	% Rectangular pulse [0 ... 0.25]
	x_nn=sum(exp(i*2*pi*(f_x'*t_n+rand_ph)))/sqrt(K); % Band-limited noise  [-0.5 ...-0.25]
	x_non=x_ns+x_ni+x_nn;		% x_non - nonuniform composite sequence
	x_non=xmax*round(x_non/xmax*lim)/lim; 	% Simulate 10-bit Analog-to-Digital Converter

% Sparse sequences (NaNs inserted): 
% x_u48 - sequence x_uni with 16 missing samples
% x_u40 - sequence x_uni with 24 missing samples 
% x_u32 - sequence x_uni with 32 missing samples
% Generate missing samples numbers for sequences x_u48, x_u40, x_u32	
	k16=[];k24=[];k32=[];
	x_u48=x_uni;x_u40=x_uni;x_u32=x_uni;
	x_48=x_uni;x_40=x_uni;x_32=x_uni;
	while length(k16)~=16, k16=find(sprandn(64,1,0.3));, end
	while length(k24)~=24, k24=find(sprandn(64,1,0.4));, end
	while length(k32)~=32, k32=find(sprandn(64,1,0.8));, end
% Insert zeros or NaNs where sample is missing
	x_48(k16)=zeros(1,16);		% x_48 - sparse data with 16 zeros (FFT input) 
	x_u48(k16)=NaN*ones(1,16);		% x_u48 - sparse data with 16 NaNs (EDFT input)
	x_40(k24)=zeros(1,24);		% x_40 - sparse data with 24 zeros (FFT input) 
	x_u40(k24)=NaN*ones(1,24);		% x_u40 - sparse data with 24 NaNs (EDFT input)
	x_32(k32)=zeros(1,32);		% x_32 - sparse data with 32 zeros (FFT input) 
	x_u32(k32)=NaN*ones(1,32);		% x_u32 - sparse data with 32 NaNs (EDFT input)	

% Marple&Kay 64-point data set x_mak /apmplified two times/ 
	x_mak=2*[1.291061 -2.086368 -1.691316 1.243138 1.641872 -0.008688 -1.659390 -1.111467 0.985908 1.991979 -0.046613...
	-1.649269 -1.040818 1.054665 1.855816 -0.951182 -1.476495 -0.212242 0.780202 1.416003 0.199282 -2.027026... 
	-0.483577 1.664913 0.614114 -0.791469 -1.195311 0.119801 0.807635 0.895236 -0.012734 -1.763842 0.309840...
	1.212892 -0.119905 -0.441686 -0.879733 0.306181 0.795431 0.189598 -0.342332 -0.328700 0.197881 0.071179...
	0.185931 -0.324595 -0.366092 0.368467 -0.191935 0.519116 0.008328 -0.425946 0.651478 -0.639978 -0.344389...
	0.814130 -0.385168 0.064218 -0.380008 -0.163008 1.180961 0.114206 -0.667626 -0.814997];
	KK=64;
	Emk=exp(-i*2*pi*(0:KK-1).'*fn);	

% Calculate DFTs and Power spectrums for plots
	disp('Calculating DFT output...');	
	dft_x_uni=fftshift(fft(x_uni,N))/K;
	dft_x_ue=x_uni*Eue/K;
	dft_x_non=x_non*En/K;
	dft_x_ne=x_non*Ene/K;
	fft_xmak=fft(x_mak,N);
	dft_x_mak=fftshift(fft_xmak)/length(x_mak);
	dft_x_48=fftshift(fft(x_48,N))/48;
	dft_x_40=fftshift(fft(x_40,N))/40;
	dft_x_32=fftshift(fft(x_32,N))/32;

       	dft_f1=20*log10(abs(dft_x_uni));
	dft_f1c=10*log10(abs(K*dft_x_uni).^2/N);
       	dft_f2=20*log10(abs(dft_x_non));
       	dft_f2c=10*log10(abs(K*dft_x_non).^2/N);
	dft_f3a=20*log10(abs(dft_x_ue));
	dft_f3b=20*log10(abs(dft_x_ne));
      	dft_f4a=20*log10(abs(dft_x_48));
      	dft_f4b=20*log10(abs(dft_x_40));
      	dft_f4c=20*log10(abs(dft_x_32));
	dft_f5a=20*log10(abs(dft_x_mak));
	ifft_f6a=real(ifft(fft_xmak));

% Calculate EDFTs and Power spectrums for plots
	disp('Calculating Extended DFT (EDFT) output...');
	[F1,S1]=edft(x_uni,N,it);
	[F10,S10]=edft(x_uni,N,1,fftshift(WT));
	[F11,S11]=nedft(x_uni,t_u,fne,10);
	S1=fftshift(S1);F1=fftshift(F1);
	[F2,S2]=nedft(x_non,t_n,fn,it);
	[F20,S20]=nedft(x_non,t_n,fn,1,WT);
	[F21,S21]=nedft(x_non,t_n,fne,15);
	[F3,S3]=edft(x_mak,N,it);
	[F4a,S4a]=edft(x_u48,N,it);
	[F4a,S4b]=edft(x_u40,N,it);
	[F4a,S4c]=edft(x_u32,N,it);

	edft_f1a=20*log10(abs(fftshift(S10)));
       	edft_f1b=20*log10(abs(S1));
       	edft_f1c=10*log10(abs(F1).^2/N);
	edft_f1d=real(F1./S1/K);
	edft_f2a=20*log10(abs(S20));
       	edft_f2b=20*log10(abs(S2));
       	edft_f2c=10*log10(abs(F2).^2/N);
	edft_f2d=real(F2./S2/K);
	edft_f3a=20*log10(abs(S11));
	edft_f3b=20*log10(abs(S21));
	edft_f3c=real(F21./S21/K/2);
      	edft_f4a=20*log10(abs(fftshift(S4a)));
      	edft_f4b=20*log10(abs(fftshift(S4b)));
      	edft_f4c=20*log10(abs(fftshift(S4c)));
       	edft_f5b=20*log10(abs(fftshift(S3)));
	ifft_f6b=real(ifft(F3));

% Calculate High-Resolution DFT of Marple&Kay sequence and Power spectrum for plot
	disp('Calculating High-Resolution DFT (HRDFT) output...');
	W=ones(1,N);			% Initial conditions for W
	for l=1:it
	fh3=edft(x_mak,N,1,W);
       	W=fh3/N;				% Weight for the next HRDFT iteration
	end 
       	hrf_f5c=20*log10(abs(fftshift(fh3/N))+eps);
	ifft_f6c=real(ifft(fh3));

% Figure 1: Plot True spectrums (red), DFT (blue) and non-iterative estimate by EDFT              
	figure(1)
	clf
	db_min=-80;			% Axis min [db]
	db_max=10;			% axis max [db]
	%True spectrum for test sequences in Subplot 1 and 2
	true_sp=20*log10(WT);		% True Power spectrum of test sequences 1 and 2
	true_psd=10*log10(WP);		% True PSD of test sequences 1 and 2	
	subplot(411)
	plot(fn,true_sp,'r-',fn,dft_f1,'b-',fn,edft_f1a)% Plot True, DFT and non-iterative EDFT Power Spectrum 
	axis([-0.5 0.5 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('UNIFORM sequence: (a) - True Power Spectrum [red], DFT [blue] and non-iterative EDFT')
	subplot(412)
	plot(fn,true_sp,'r-',fn,dft_f1,'b-',fn,edft_f1b)	% Plot True, DFT and EDFT Power Spectrum
	axis([-0.5 0.5 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('(b) - True Power Spectrum [red], DFT [blue] and EDFT 10th iteration')
	subplot(413)
	plot(fn,true_psd,'r-',fn,dft_f1c,'b-',fn,edft_f1c) % Plot True, DFT and EDFT Power Spectrum
	axis([-0.5 0.5 db_min db_max+30])
	xlabel('Frequency')
	ylabel('PSD [dB]')
	title('(c) - True PSD [red], DFT [blue] and EDFT 10th iteration')
	subplot(414)
	plot(fn,ones(1,N),'b-',fn,edft_f1d)		% Plot DFT and EDFT relative resolution
	axis([-0.5 0.5 0 ceil(N/K)])
	xlabel('Frequency')
	ylabel('Resolution')
	title('(d) - Relative frequency resolution of DFT [blue] and EDFT 10th iteration')
	
% Figure 2: Plot True spectrums (red) and estimates by DFT (blue) and EDFT 10th iteration           
	figure(2)
	clf
	subplot(411)
	plot(fn,true_sp,'r-',fn,dft_f2,'b-',fn,edft_f2a)% Plot True, DFT and non-iterative EDFT Power Spectrum
	axis([-0.5 0.5 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('NONUNIFORM sequence: (a) - True Power Spectrum [red], DFT [blue] and non-iterative EDFT')
	subplot(412)
	plot(fn,true_sp,'r-',fn,dft_f2,'b-',fn,edft_f2b)	% Plot True, DFT and EDFT Power Spectral Density
	axis([-0.5 0.5 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('(b) - True Power Spectrum [red], DFT [blue] and EDFT 10th iteration')
	subplot(413)
	plot(fn,true_psd,'r-',fn,dft_f2c,'b-',fn,edft_f2c) % Plot True, DFT and EDFT Power Spectral Density
	axis([-0.5 0.5 db_min db_max+30])
	xlabel('Frequency')
	ylabel('PSD [dB]')
	title('(c) - True PSD [red], DFT [blue] and EDFT 10th iteration')
	subplot(414)
	plot(fn,ones(1,N),'b-',fn,edft_f2d)		% Plot DFT and EDFT relative resolution
	axis([-0.5 0.5 0 ceil(N/K)])
	xlabel('Frequency')
	ylabel('Resolution')
	title('(d) - Relative frequency resolution of DFT [blue] and EDFT 10th iteration')

% Figure 3: Plot Extended True spectrums (red) and estimates by DFT (blue) and EDFT              
	figure(3)
	clf
	true_e=20*log10(WE);		% Extended True spectrum of test sequences 1 and 2
	subplot(311)
	plot(fne,true_e,'r-',fne,dft_f3a,'b-',fne,edft_f3a) % Plot True, DFT and EDFT Power Spectrum
	axis([-1 1 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('UNIFORM sequence: (a) - True Power Spectrum [red], DFT [blue] and EDFT 10th iteration')
	subplot(312)
	plot(fne,true_e,'r-',fne,dft_f3b,'b-',fne,edft_f3b) % Plot True, DFT and EDFT Power Spectrum
	axis([-1 1 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('NONUNIFORM sequence: (b) - True Power Spectrum [red], DFT [blue] and EDFT 15th iteration')
	subplot(313)
	line([0.35; 0.35],[N/K; 0],'linestyle',':','color','r')
	axis([-1 1 0 ceil(N/K)])
	hold on
	plot([-0.5 0;-0.25 0.25],[1 1;1 1],'r:');
	plot(fne,ones(1,2*N)/2,'b-',fne,edft_f3c)	% Plot DFT and EDFT relative resolution
	hold off
	xlabel('Frequency')
	ylabel('Resolution')
	title('(c) - Relative frequency resolution of DFT [blue] and EDFT 15th iteration')

% Figure 4: Plot True spectrums (red), DFT (blue) and EDFT estimate for sequences wit NaNs inserted             
	figure(4)
	clf
	subplot(311)
	plot(fn,true_sp,'r-',fn,dft_f4a,'b-',fn,edft_f4a)	% Plot True, DFT and EDFT Power Spectrum for sequence with 16 NaNs 
	axis([-0.5 0.5 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('True Power Spectrum [red], DFT [blue] and EDFT 10th iteration. UNIFORM sequence: (a) - 16 lost samples' )
	subplot(312)
	plot(fn,true_sp,'r-',fn,dft_f4b,'b-',fn,edft_f4b)	% Plot True, DFT and EDFT Power Spectrum for sequence with 24 NaNs
	axis([-0.5 0.5 db_min db_max])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('UNIFORM sequence: (b) - 24 lost samples')
	subplot(313)
	plot(fn,true_psd,'r-',fn,dft_f4c,'b-',fn,edft_f4c) % Plot True, DFT and EDFT Power Spectrum for sequence with 32 NaNs
	axis([-0.5 0.5 db_min db_max+30])
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('UNIFORM sequence: (b) - 32 lost samples')

% Figure 5: Plot True spectrums (red) and estimates by DFT, EDFT and HRDFT of Marple&Kay sequence      
	figure(5)
	clf
	subplot(311)
	% True spectrum for Marple&Kay test sequence
	fr_true=[0.1 0.2 0.21; 0.1 0.2 0.21;];
	amp_true=[-20 0 0; db_min db_min db_min];
	col_n=15*log10(0.1*(sin(pi*(eps:0.001:1))).^2); % Colored noice as described in [1]
	line(fr_true,amp_true,'linestyle','-','color','r')
	axis([0 0.5 db_min db_max])
	hold on
	plot(0.2:0.0003:0.5,col_n,'-r')  
	plot(fn,dft_f5a)			% Plot DFT Power Spectrum
	hold off
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('Marple&Kay sequence: (a) - True Power Spectrum [red] and DFT')
	subplot(312)
	line(fr_true,amp_true,'linestyle','-','color','r')
	axis([0 0.5 db_min db_max])
	hold on
	plot(0.2:0.0003:0.5,col_n,'-r')  
	plot(fn,edft_f5b)			% Plot EDFT Power Spectrum
	hold off
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('(b) - True Power Spectrum [red] and EDFT 10th iteration')
	subplot(313)
	line(fr_true,amp_true,'linestyle','-','color','r')
	axis([0 0.5 db_min db_max])
	hold on
	plot(0.2:0.0003:0.5,col_n,'-r')  
	plot(fn,hrf_f5c)			% Plot HRDFT Power Spectrum
	hold off
	xlabel('Frequency')
	ylabel('Power [dB]')
	title('(c) - True Power Spectrum [red] and HRDFT 10th iteration')

% Figure 6: Plot Marple&Kay sequence and reconstructed data by DFT, EDFT and HRDFT.
	figure(6)
	clf
	subplot(311)
	x_max=ceil(max(abs(ifft_f6c)));
	plot(ifft_f6a)			% Plot Inverse DFT sequence
	axis([0 N+1 -x_max x_max])
	hold on
	plot(x_mak,'b-')			% Plot Marple&Kay sequence
	hold off
	xlabel('Sample number')
	ylabel('Magnitude')
	title('Marple&Kay sequence: (a) - True (blue) and reconstructed by Inverse DFT')
	subplot(312)
	plot(ifft_f6b)			% Plot Inverse EDFT sequence
	axis([0 N+1 -x_max x_max])
	hold on
	plot(x_mak,'b-')			% Plot Marple&Kay sequence
	hold off
	xlabel('Sample number')
	ylabel('Magnitude')
	title('(b) - True (blue) and reconstructed by Inverse EDFT')
	subplot(313)
	plot(ifft_f6c)			% Plot Inverse HRDFT sequence
	axis([0 N+1 -x_max x_max])
	hold on
	plot(x_mak,'b-')			% Plot Marple&Kay sequence
	hold off
	xlabel('Sample number')
	ylabel('Magnitude')
	title('(c) - True (blue) and reconstructed by Inverse HRDFT')

%% NOTE: Code below is commented because of the power spectrums coincide with the EDFT output (Figure 5b).
%% Calculate Capon filter output and GWLS solution of Marple&Kay sequence and Power spectrums for plot
%	disp('Calculating Capon filter output...');
%	x_mr=(flipud(x_mak.')).';		% Re-ordering Marple&Kay sequence for Capon filter input
%	W=ones(1,N);			% Initial conditions for iterative algorithm
%	for l=1:it
%	r=ifft(fftshift(W));			% Calculate autocorrelation function by applying ifft
%	RT=toeplitz(r(1:KK)).';			% Compose and transpose the autocorrelation matrix
%	ER=inv(RT)*conj(Emk);
%	sc1=(x_mr*ER)./sum(Emk.*ER);		% Calculate amplitude spectrum for iteration (l) 
%	W=sc1.*conj(sc1);		 	% Weight for the next iteration 
%	end
%	cap_f6a=20*log10(abs(sc1));
%	disp('Calculating GWLS solution output...');
%	W=ones(1,N);			% Initial conditions for iterative algorithm
%	for l=1:it
%	r=ifft(fftshift(W));			% Calculate autocorrelation function by applying ifft
%	RT=toeplitz(r(1:KK)).';			% Compose and transpose the autocorrelation matrix
%	ER=Emk.'*inv(RT);
%	sc2=((ER*x_mak.').')./sum(conj(Emk).*ER.'); % Calculate amplitude spectrum for iteration (l) 
%	W=sc2.*conj(sc2);		 	% Weight for the next iteration 
%	end
%	gwls_f6b=20*log10(abs(sc2));
%% Figure 7: Plot True spectrum (red) and estimates by DFT (blue) and Capon filter approach       
%	figure(7)
%	clf
%	subplot(211)
%	line(fr_true,amp_true,'linestyle','-','color','r')
%	axis([0 0.5 db_min db_max])
%	hold on
%	plot(0.2:0.0003:0.5,col_n,'-r')  
%	plot(fn,dft_f4a,'b-',fn,cap_f6a)		% Plot DFT and Capon filter Power Spectrums
%	hold off
%	xlabel('Frequency')
%	ylabel('Power [dB]')
%	title('Marple&Kay sequence: (a) - True Power Spectrum [red], DFT [blue] and Capon filter 10th iteration')
%	subplot(212)
%	line(fr_true,amp_true,'linestyle','-','color','r')
%	axis([0 0.5 db_min db_max])
%	hold on
%	plot(0.2:0.0003:0.5,col_n,'-r')  
%	plot(fn,dft_f4a,'b-',fn,gwls_f6b)		% Plot DFT and GWLS solution Power Spectrums
%	hold off
%	xlabel('Frequency')
%	ylabel('Power [dB]')
%	title('(b) - True Power Spectrum [red], DFT [blue] and GWLS solution 10th iteration')