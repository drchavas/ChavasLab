%FourierDecomposition_example2.m

%Source: http://faculty.olin.edu/bstorey/Notes/Fourier.pdf

clear;close('all');clc

N = 10^6;

%%Sample a simple sine wave with amplitude 1
T = 1;  %[s]; length of time interval
t = [0:N-1]'/N; %[-]; time in fractions of interval
t = t*T;    %[s]; actual time
freq1 = 10; %[Hz]; frequency
amp1 = 1;   %[-]; amplitude
f = amp1*sin(2*pi*freq1*t); %sine wave with frequency [freq1]
c_n = fft(f);   %MATLAB FFT coefficients

%Note: The ordering of the corresponding MATLAB frequencies is as follows [0 1 2 3 4 -3 -2 -1]. 
n = [0:N/2 -N/2+1:-1]';
%Note that 8 discrete data points yields 8 Fourier coe±cients and the
%  highest frequency that will be resolved is N/2.
%  The order that the coefficients come in is often called
%  reverse-wrapa-around order. The first half of the list of numbers is the
%  positive frequenciesa nd the second half is the negative frequencies.
%  The real part of the FFT corresponds to the cosines series and the
%  imaginary part corresponds to the sine.

%When taking an FFT of a real number data set (i.e. no complex numbers in
%the original data), the positive and negative frequencies turn out to be
%complex conjugates.

%%MATLAB FFT returns data that needs to be divided by N/2 to get the
%%coefficients that we used in the sin and cos series. To obtain the
%%coefficients a_n and b_n we use the following formulas, where c_n is the
%%coefficient of the MATLAB FFT:
a_n = -imag(c_n)/(N/2); %0<n<N/2; imaginary part = sine
b_n = real(c_n)/(N/2); %0<n<N/2; real part = cosine
% [n a_n b_n]


%% Power Spectrum
%Note: Usually we only care how much information is contained at a
%   particular frequency and we don?t really care whether it is part of
%   the sine or cosine series. Therefore, we are interested in the
%   absolute value of the FFT coefficients. The absolute value will
%   provide you with the total amount of information contained at a given
%   frequency, the square of the absolute value is considered the power
%   of the signal. Remember that the absolute value of the Fourier
%   coefficients are the distance of the complex number from the origin
p = abs(fft(f))/(N/2);    %% absolute value of the fft
p = p(1:N/2).^2;    %% take the positve frequency half, only
freq = [0:N/2-1]/T; %[Hz]; frequency in real time

semilogy(freq,p);   %plot power spectrum
axis([0 30 0 1]);