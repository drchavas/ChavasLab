% Demonedft.m demonstrates performance of NEDFT algorithm by iterations in case of nonuniform sampling and
% provide comparison with traditional nonuniform discrete Fourier transform (DFT) output. NEDFT call line
% is [F,S,Stopit]=nedft(X,tk,fn,1,W), where DEMO data X are calculated as sum of the following components:
% - real sinusoid at normalized frequency 0.137, magnitude 2,
% - complex exponent at frequency 0.137+1/(4*K), magnitude 1 (K=64), where K=64 is length of sequence X,
% - mean value (frequency 0), magnitude 1
% - pulse in frequency ranges [-0.137 0] and
% - noise generated by randn function (SNR~30 dB).
% Number of frequencies in fn (length of Fourier transform) N=512.
% Number of samples K=64 (uniform and nonuniform sequences) or KK=32 (nonuniform sequence, 1/2 of samples missing).
% So, the true spectrum of DEMO sequence (without noise part) looks like greeting 'HI', where 'I' actually
% are two close peaks at normalized frequencies 0.137 and 0.137+1/(4*K).
%
% In the first, demonedft.m plots real parts of the input sequences:
% - Subplot(311) display DEMO data X1 consisting of 64 samples, tk=[1:K];
% - Subplot(312) outputs nonuniformly sampled DEMO data X2, tkn2=[1:K]+0.5*(rand(1,K)-0.5);
% - Subplot(313) showing sparse (32 samples) nonuniformly sampled data X3, tkn3=tkn2(1:2:K); 
% To proceed you should select one of these sequences - enter [1] to see DEMO for 64 uniformly sampled data or
% select [2] or [3], where 64 or 32 nonuniformly sampled data will be processed by DFT and NEDFT algorithms.	
%
% In the next, fifteen (I=15) NEDFT iterations are performed. You will be asked to strike any key after 1,2,3,5
% and 10 iteration. Plots in blue colour are equal to DFT analysis and can be used for comparison with NEDFT output
% (green colour). To resolve two close sinusoids (at 1/4 of Nyquist), NEDFT should increase resolution at least 
% 4 times (in comparison with DFT, see Subplot(413)).
%
% NEDFT and DFT outputs are displayed in four subplots:
% - Subplot(221) show Power Spectral Density as 10*log10(abs(F).^2/N)) in normalized frequencies range [-0.5 0.5[.
% - Subplot(222) display Power Spectrum as 10*log10(abs(S).^2). The plot indicate power of sinusoids in X.
% - Subplot(413) plot division F./S/K titled as Relative frequency resolution and expose how the frequency resolution
% of NEDFT algorithm is changed along the frequency axes in respect to DFT analysis. The first iteration showing
% that DFT and NEDFT analysis have almost equal resolution for all frequencies in range [-0.5 0.5[. During the next
% NEDFT iterations we can discover that the relationship F./S/K for NEDFT analysis is not that simple as for DFT.
% Although sum(F./S/K)=N and remains constant for all iterations, NEDFT have ability to increase resolution 
% around the powerful narrowband components (sinusoids) and decrease resolution at frequencies where data X have
% weak power components. NEDFT is called as high resolution method and that's true, but with the following remark-
% NEDFT keeping the same 'summary' resolution as DFT or, in other words, squares under curves for DFT (blue colour)
% and NEDFT (green colour) analysis are equal. The maximum frequency resolution is limited by value of division N/K.
% For example, if K=64 and N=512, then NEDFT can potentially improve frequency resolution 512/64=8 times. 
% - Subplot(414) plot real parts of reconstructed sequences obtained as result of inverse fast Fourier transform
% to outputs of DFT and NEDFT - ifft(fftshift(dft_out)) and ifft(fftshift(F)), correspondingly. It is well known
% that reconstruction of data by applying ifft(fft(X,N)) for N>length(X), return sequence where initial data X are
% padded with zeros to length N. That result is true also for DFT and for the first iteration of NEDFT. The next
% iterations are showing ability of NEDFT to obtain reconstructed data consisting of input sequence X values plus
% non-zero extrapolation of X to length N. As shown in subplot 414, sequence X are extrapolated in both directions
% - backward and forward. Also, for nonuniformly sampled input sequences [X2] and [X3], the reconstructed data are
% re-sampled at uniform time set tn=1:N.	
%
% E-mail:	Vilnislp@gmail.com
% Reference: 	V.Liepin'sh, "A spectral estimation method of nonuniformly sampled band-limited signals. 
%		Automatic Control and Computer Sciences, Vol.28, No.2, pp.66-77, 1994.

	clear
	disp('NEDFT DEMO program started...');

% Set the length for input sequence X (K), DFT (N) and number of iterations (I).
	
	K=64;
	N=512;
	I=15;

% Relative central frequency of the test signal (value in range ]0,0.5]) and SNR [dB].
	fc=0.137;			% frequency of sinusoids	
	SNR=30;			% signal to noise ratio [dB]
	sigma=sqrt((2+1+1)/10^(SNR/10));% sigma for given SNR
			
% Uncomment to generate always the same signal and/or noise
	%rand('seed',777);
	%randn('seed',777);

% Generate uniformly and nonuniformly spaced time, frequency, noise vectors and initial phases of sinusoids.

	tk1=1:K;			% uniform time vector (K samples, sampling period T=1)
	tkn2=[1:K]+0.5*(rand(1,K)-0.5);	% nonuniform time vector (K samples, mean sampling period Ts=1)
	tkn3=tkn2(1:2:K); 		% sparse nonuniform time vector (K/2 samples, mean sampling period Ts=2)
	fn=[-ceil((N-1)/2):floor((N-1)/2)]/N;% uniform normalized frequency set (sampling frequency- 1)
	tn=1:N;			% uniform time vector (N samples, sampling period T=1)
	ti1=tk1-K/2-0.5;		% uniform centerred time vector for pulse
	tin2=tkn2-K/2-0.5;		% nonuniform centerred time vector for pulse
	cno=sigma*(randn(1,K)+i*randn(1,K))/2;	% complex noise
	ph1=2*pi*rand;		% initial phase of real sinusoid at frequency fc
	ph2=2*pi*rand;		% initial phase of complex exponent at frequency fc

% Generate uniformly sample input sequence X1 (K samples).
	XC1=2*sin(2*pi*fc*tk1+ph1);				% sinusoid
	XC2=exp(i*(2*pi*(fc+1/4/K)*tk1+ph2));			% complex exponent
	XC3=ones(1,K);					% mean value
	XC4=sin(pi*fc*ti1)./(ti1+eps).*exp(-i*pi*fc*ti1);		% rectangular pulse
	X1=XC1+XC2+XC3+XC4;				% signal w/o noise			
	SNR1=10*log10((X1*X1')/(cno*cno'));			% calculated SNR
	X1=X1+cno;					% signal+noise

% Generate nonuniformly sampled sequence X2 (K samples)
       	XN1=2*sin(2*pi*fc*tkn2+ph1);				% sinusoid
	XN2=exp(i*(2*pi*(fc+1/4/K)*tkn2+ph2));			% complex exponent
	XN3=ones(1,K);					% mean value
	XN4=sin(pi*fc*tin2)./(tin2+eps).*exp(-i*pi*fc*tin2);		% rectangular pulse
	X2=XN1+XN2+XN3+XN4;				% signal w/o noise		
	SNR2=10*log10((X2*X2')/(cno*cno'));			% calculated SNR
	X2=X2+cno;					% signal+noise

% Generate nonuniformly sampled sequence X3 (taking half of X2 samples)
	
	X3=X2(1:2:K);

% Plot the real part of data X1, X2 and X3.

	figure(1)
	clf
	Xmin=min(real([X1 X2]))-1;
	Xmax=max(real([X1 X2]))+1;
	subplot(311)
	stem(tk1,real(X1))
	axis([0 K+1 Xmin Xmax])
	xlabel('Sample number')
	ylabel('real(X1)')
	title(['Input sequence X1 - uniform sampling: ',int2str(K),' samples, SNR=',int2str(round(SNR1))])
	subplot(312)
	stem(tkn2,real(X2))
	axis([0 K+1 Xmin Xmax])
        	xlabel('Sample number')
	ylabel('real(X2)')
	title(['Input sequence X2 - nonuniform sampling: ',int2str(K),' samples, SNR=',int2str(round(SNR2))])
	subplot(313)
	stem(tkn3,real(X3))
	axis([0 K+1 Xmin Xmax])
        	xlabel('Sample number')
	ylabel('real(X3)')
	title(['Input sequence X3 - sparse nonuniform sampling: ',int2str(K/2),' samples'])

% Selecting of sequence for NEDFT input and calculating DFT output

	disp('DEMO data are calculated as sum of the following components:');
        	disp(['   - real sinusoid at normalized frequency ',num2str(fc),', magnitude 2;']);
        	disp(['   - complex exponent at frequency ',num2str(fc),'+1/(4*',int2str(K),'), magnitude 1;']);
        	disp('   - mean value (frequency 0), magnitude 1;');
	disp(['   - rectangular pulse in frequency ranges [-',num2str(fc),' 0];']);
        	disp(['   - noise generated by randn function (SNR~',int2str(SNR),' dB).']);
	disp(' ');
	disp('Input [1] to select uniformly sampled sequence X1');
	disp('Input [2] to select nonuniformly sampled sequence X2');
	disp('Input [3] to select sparse nonuniformly sampled sequence X3');
	DEMO=input('Select input sequence for DEMO:');  
	if DEMO==1,
	    X=X1;tk=tk1;KK=K;		% X=X1 - uniform data (K samples)
	    E=exp(-i*2*pi*tk1.'*fn);		% E - complex exponents matrix
	    dft_out=X*E; 			% Calculate DFT output		
        	elseif DEMO==2,
	    X=X2;tk=tkn2;KK=K;		% X=X2 - nonuniform data (K samples)
	    E=exp(-i*2*pi*tkn2.'*fn);		% E - complex exponents matrix
	    dft_out=X*E; 			% Calculate DFT output							
	elseif DEMO==3,
	    X=X3;tk=tkn3;KK=K/2;		% X=X3 - sparse nonuniform data (K/2 samples)
	    E=exp(-i*2*pi*tkn3.'*fn);		% E - complex exponents matrix
	    dft_out=X*E; 			% Calculate DFT output		
	else
	   disp('...end of DEMO')		% End DEMO if other symbol entered
	   return
	end

% Calculate outputs for DFT used in 4 plottings.

	sub1_PWD=10*log10(abs(dft_out).^2/N);
	sub2_PS=10*log10(abs(dft_out/KK).^2);
	sub3_FR=ones(1,N)*KK/K;
	sub4_RD=real(ifft(fftshift(dft_out)));

% Set values for input argument W used for the first NEDFT iteration.

	W=ones(1,N);

% Get NEDFT outputs F and S for iteration it.

	for it=1:I,

	[F,S,Stopit]=nedft(X,tk,fn,1,W);

% Break NEDFT iterations if get ill conditioned matrix or results may be inaccurate.
	if Stopit(2)==1,it=it-1;break,end 

% Calculate outputs for NEDFT iteration it used in 4 plottings.

	sub1=10*log10(abs(F).^2/N);
	sub2=10*log10(abs(S).^2);
	sub3=real(F./S/K);
	sub4=real(ifft(fftshift(F)));
	frsig=[-fc 0 fc fc+1/4/K];

% Plots Power Spectral Densities in subplot221.

	figure(2)	
	clf
	subplot(221)
	plot([frsig;frsig],[10*log10(N)*[1 1 1 1];-100*[1 1 1 1]],'r:');
	hold on
	plot([-fc;0],[10*log10(pi^2/N);10*log10(pi^2/N)],'r:');
	plot(fn,sub1,'g-',fn,sub1_PWD,'b-')
	hold off
	axis([-0.5 0.5 -100 50])
	xlabel('Normalized frequency')
	ylabel('10*log[abs(F).^2/length(F)]')
	title(['Power Spectral Density (iter.',int2str(it),')'])
	
% Plots Power Spectrums in subplot222.

	subplot(222)
	plot([frsig;frsig],[[0 0 0 0];-100*[1 1 1 1]],'r:');
	hold on
	plot([-fc;0],[20*log10(pi/K);20*log10(pi/K)],'r:');
        	plot(fn,sub2,'g-',fn,sub2_PS,'b-')
	hold off
	axis([-0.5 0.5 -100 50])
	xlabel('Normalized frequency')
	ylabel('10*log[abs(S).^2]')
	title(['Power Spectrum (iter.',int2str(it),')'])

% Plots Frequency Resolutions in subplot413.

	subplot(413)
	plot([frsig;frsig],[[0 0 0 0];N/K*[1 1 1 1]],'r:');
        	hold on
	plot(fn,sub3,'g-',fn,sub3_FR,'b-')
        	plot([-fc;0],[1;1],'r:');
	hold off
	axis([-0.5 0.5 0 N/K])
	xlabel('Normalized frequency')
       	ylabel(['(F./S)/length(X',int2str(DEMO),')'])
	title(['Relative frequency resolution (iter.',int2str(it),')'])

% Plots True and Reconstructed data in subplot414.

	subplot(414)
	plot(tk1+1,real(X1),'r-',tn,sub4,'g-',tn,sub4_RD,'b-')
	axis([1 N Xmin Xmax])
	xlabel('Sample number')
	ylabel('real(ifft(F))')
	title(['True (red) and reconstructed sequences by NEDFT (green, iter.',int2str(it),') and DFT (blue)'])

% Waiting for keyboard action

	drawnow
        	if it==1
	disp(' ');
        	disp('Plots True Spectrum [red color], NEDFT [green colour] and DFT [blue colour].');
        	end
        	if it==1|it==2|it==3|it==5|it==10,
	disp(['Calculate [F,S]=nedft(X',int2str(DEMO),',tk,fn,',int2str(it),') and ifft(F). Strike any key to continue.']);   
	pause
	end

% Calculate NEDFT input argument W for next iteration.

	W=S;

	end
	disp(['Calculate [F,S]=nedft(X',int2str(DEMO),',tk,fn,',int2str(it),') and ifft(F). ...end of DEMO.']);