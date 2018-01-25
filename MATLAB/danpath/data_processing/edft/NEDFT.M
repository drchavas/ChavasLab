function [F,S,Stopit]=nedft(X,tk,fn,I,W)

% NEDFT - Nonuniform Extended Discrete Fourier Transform.
%
% SYNTAX
% a. Mandatory inputs/outputs	
%    F=nedft(X,tk,fn)		
%	Function NEDFT returns discrete Fourier transform F of input sequence X sampled at arbitrary
%	selected time moments tk: X(tk) >>> F(fn), where frequencies fn, in general, also may selected
%	arbitrary. If fn is less than X, input sequences X and tk will be truncated.
% b. Mandatory and optional inputs/outputs
%    [F,S,Stopit]=nedft(X,tk,fn,I,W)		
%    I	Optional input parameter I can be used for limiting maximum number of iterations. If I is not 
%	specified in input arguments, default value for I is set by parameter 'Miteration', that is, 
%	nedft(X,tk,fn)=nedft(X,tk,fn,Miteration). To complete iteration process faster, the value for 
%	'Miteration' should be decreased.
%    W	Input weight vector W, if specified, override the default values W=ones(size(fn)). W must have
%	at least length(X) nonzero elements.	
%    S	The second output argument S represents the Amplitude spectrum. Peak values of abs(S) can be
%	used for estimate amplitudes of sinusoids in the input sequence X.
%  Stopit is an informative output parameter. The first row of Stopit showing the number of performed iteration,
%	the second row indicate breaking of iteration reason and may have the following values:
%	0- Maximum number of iteration performed.
%	1- Sum of outputs division sum(F./S) is not equal to K*N within Relative deviation 'Rdeviat'. 
%	the calculations is interrupted because of results could be inaccurate. If this occur in the first
%	NEDFT iteration, then outputs F and S are zeros.
%	2- Relative threshold 'Rthresh' reached. To complete iteration process faster, the value for
%	'Rthresh' should be increased.
% ALGORITHM
%    Input: 
%	X- input sequence
%	E- complex exponents matrix (Fourier transform basis) - E=exp(-i*2*pi*tk.'*fn);
%	I- (optional) number of maximum iteration.
%	W- (optional) weight vector W. If not specified, W = ones(1,size(fn)) used for the first iteration.
%    Output F and S for each NEDFT iteration are calculated by following formulas:
%	1. R=E*diag(W/N)*E';
%	2. F=W.*(X*inv(R)*E); 
%	   S=(X*inv(R)*E)./diag(E'*inv(R)*E).';
%	3. W=S.*conj(S); - the weight vector W for the next iteration.
%    A special case: if length(X) is equal to length(fn), the NEDFT output do not depend on selected weight 
%	vector W and is calculated in non-iterative way.   
% Tips for selection of mandatory NEDFT inputs X(tk) and fn:
%	1. Input sequence X(tk) for NEDFT can be sampled uniformly or nonuniformly. Uniform sampling
%	can be considered as a special case of nonuniform sampling, where tk=[0,1,...,K-1]*T and T is 
%	sampling period. Nonuniform sampling can be realized in many different ways, like as:
%	- uniform sampling with randomly missed samples (known as sparse data);
%	- uniform sampling with missed data segments (known as gapped data);
%	- uniform sampling with jitter: tk=([0,1,...,K-1] + jitter*rand(1,K))*Ts, where value for jitter is selected
%	in range [0...1[ and Ts is the mean sampling period;   
%	- additive nonuniform sampling: tk=tk-1 + (1+jitter*(rand-0.5))*Ts, k=1,...K-1, t0=0;
%	- signal dependent sampling, e.g, level-crossing sampling, etc... .
%	2. Frequencies for fn can be selected arbitrary. This mean, that user can choose not only the length
%	of NEDFT (number of frequencies in fn), but also the way how to distribute frequencies along the 
%	frequency axis. On other hand, to get adequate sequence X representation, frequencies fn should
%	be selected to cover overall range, where the input sequence X spectrum is supposed to be found,
%	otherwise, in result of NEDFT, all components having spectra outside fn will be incorporated.
%	Note that fn should contain negative frequencies too, and for a real value X(tk) analysis each positive
%	frequency in fn should have corresponding negative one. 
%	3. Frequencies for vector fn can be added in any order. Therefore it is possible to combine different 
%	frequency sets in one or just add individual frequencies of interest to fn, e.g, fn=[fn1 fn2 f1 f2], where 
%	fn1 and fn2 are different frequency sets, f1,f2 - specific frequencies. NEDFT outputs will be calculated
%	accordingly- F(fn)=[F(fn1) F(fn2) F(f1) F(f2)], S(Fn)=[S(Fn1) S(fn2) S(f1) S(f2)]. 
% FEATURES
%	1. NEDFT output F(fn) is the discrete Fourier transform of sequence X(tk).
%	The Power Spectral Density function of nonuniform sequence X(tk) can be estimated by the following
%	formula: abs(F).^2/(N*Ts), Ts - mean sampling period.
%	2. In general, the function Y=inedft(F,fn,tn) (see attached program) is used to calculate the reconstructed 
%	sequence Y(tn). If frequencies fn are selected on the same grid as used by FFT algorithm, then ifft(F) can
%	be applied to get uniformly re-sampled and extrapolated to length(fn) version of input sequence X(tk).
%	3. NEDFT output S(fn) estimate amplitudes and phases of sinusoidal components in sequence X(tk). 
%	4. NEDFT can increase frequency resolution length(fn)/length(X) times. Division of outputs 1/(Ts*(F./S))
%	demonstrate the frequency resolution of NEDFT. he following is true for any NEDFT iteration: 
%		0<F./S<=length(fn),
%		sum(F./S)=length(fn)*length(X).	
%	5. If input arguments are matrixes, the NEDFT operation is applied to each column.
%
%	See also FFT, IFFT, FFTSHIFT, EDFT, INEDFT.
%
% E-mail: 		vilnislp@gmail.com
% Reference: 	V.Liepin'sh, "A spectral estimation method of nonuniformly sampled band-limited signals.
%		Automatic Control and Computer Sciences, Vol.28, No.2, pp.66-77, 1994.

%======================= Set default parameters for NEDFT ============================
Miteration=30;				% Limit for maximum number of iteration (Stopit 0). 
Rdeviat=0.0005;				% Value for relative deviation (Stopit 1).
Rthresh=0.0001;				% Value for relative threshold (Stopit 2).
%======================= Check NEDFT input arguments ==============================
if nargin<3,error('Not enough input arguments. See help nedft.'),end
if sum(any(isinf(X)))sum(any(isnan(X))), error('Input argument X contain Inf or NaN.  See help nedft.'), end
if size(X,1)==1,				% Check size of input sequence X.
    trf=0;
else
    X=X.'; tk=tk.'; fn=fn.'; trf=1;
end
[L K]=size(X);				% K - length of input sequence X.
if size(tk,1)~=L | size(tk,2)~=K, error('Size of input arguments X and tk must be equal. See help nedft.'), end
if size(fn,1)~=L, error('Incorrect size of input argument fn. See Help nedft.'), end
N=size(fn,2);				% N - length of DFT.
if N<K,					% Truncate sequence X if N<K.
    X=X(:,1:N); tk=tk(:,1:N); K=N;
end
if nargin<4,				% Set value for maximum number of iterations.
    I=Miteration;				% Default value for I.
else
    if isempty(I),I=Miteration;end, I=floor(I(1));	% Check input argument I.
end
if nargin>4,				% Check of input argument W. 
    if trf==1,W=W.';end
    if (size(W,2)~=N)|(size(W,1)~=L),error('Incorrect size of input argument W. See help nedft.'), end   
    W=W.*conj(W); 
    if any(find(sum(W>0)<K)), error('Too many zeros in input argument W. See help edft.'), end    
    else
    W=ones(L,N);				% Default values for W.
end
%======================= Check for a special cases ==============================
if K==N, I=1; W=ones(L,K); end			% If K=N, perform just one NEDFT iteration.
%======================= Set default values for NEDFT output arguments ===============
F=zeros(L,N); S=zeros(L,N);			% Fill zeros in output matrixes F and S.
Stopit=[I*ones(1,L); zeros(1,L)];			% Stopit 0: Set values for default Stopit.
%======================= Calculate NEDFT for each X column l ======================
for l=1:L,
E=exp(-i*2*pi*tk(l,:).'*fn(l,:));			% Calculate the complex exponents matrix E.
    for it=1:I,				% Start iterations...
% Calculate the correlation matrix R by using a loop structure.
	for n=1:K, 
	    for k=n:K,
		R(k,n)=sum(W(l,:).*conj(E(n,:)).*E(k,:))/N;
		if n~=k, 
		    R(n,k)=conj(R(k,n));
		else
		    R(n,n)=real(R(n,n));
		end
	    end
	end			     
% Calculate the correlation matrix R by using vectorized form and RE=R\E (an alternative approach).
%	R=E*diag(W(l,:)/N)*E';
%     	RE=R\E; 
% Calculate RE=inv(R)*E and ERE=diag(E'*inv(R)*E).'=sum(conj(E).*RE).
	RE=inv(R)*E; 
	ERE=sum(conj(E).*RE);
% Stopit 1: Break iterations if sum(F./S) is not equal to N*K.
	if abs(ERE*W(l,:).'/N/K-1)>Rdeviat, Stopit(:,l)=[it-1; 1]; break, end
% Calculate outputs for iteration (it): N-point NEDFT (F) and Amplitude Spectrum (S).
	F(l,:)=X(l,:)*RE;
	S(l,:)=F(l,:)./ERE;
	F(l,:)=F(l,:).*W(l,:);
% Calculate weight (W) for the next iteration.
        	W(l,:)=S(l,:).*conj(S(l,:));
% Stopit 2: Break iterations if relative threshold reached.
        	SW(it)=sum(W(l,:));
        	if it>1,
            	    thit=abs(SW(it-1)-SW(it))/SW(1);
            	    if thit<=Rthresh, Stopit(:,l)=[it; 2]; break, end
        	end		
    end					% ... end iterations.
end
%======================= Adjust size of NEDFT output ==============================
if trf==1,F=F.';S=S.';end			% Adjust size of NEDFT outputs.