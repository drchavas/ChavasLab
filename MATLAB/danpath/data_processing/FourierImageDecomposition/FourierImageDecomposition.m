%% Demonstration of fourier decomposition. 
% Goal: to visualize that (2D) images can be decomposed 
% into sinusoidal waves. 
% ----------
% Implementation: At each 'iteration' the most prominent 
% wave* is taken away from the input image, and set aside.
% Where the residual image gradually loses contrast, the removed 
% waves together start forming a new reconstructed 
% image (initially low pass filtered). After 20 
% waves the rest of the waves is removed in groups of 
% increasing size.
% (*together with its 'conjugate' wave). 
% ----------
% Author: Job Bouwman
% date: 20 july 2013 
% ----------

close all; clear all; clc;

%% user settings:
% speed in which the demo is carried out:
demoSpeed = 20;  % somewhere between 1 and 10

% The selected image:
% I = imread('moon.tif');
I = imread('cameraman.tif');

% The scale in which it is performed;
scale = 0.5;   % either {1 , 0.5 or 0.25}  


%% The demo:
I = imresize(mean(double(I),3), scale);
I = I - mean(I(:));
N = size(I);
minI = min(I(:));
maxI = max(I(:));
viewRange = [(1.25*minI - 0.25*maxI), (1.25*maxI - 0.25*minI)];

figure(1);
subplot(2,1,1); imagesc(I, viewRange); axis equal tight; colormap gray;
title('2D Fourier decomposition', 'FontSize',18); 
xlabel('original image', 'FontSize',12); 

FT_I = fftn(I);
FT_new = zeros(size(FT_I));
ABS_FT_I = abs(FT_I);
FT_cur = zeros(size(FT_new));
newIm = zeros(size(FT_new));

n = 1; nWaves = 1;
while n < numel(I)
    FT_cur = 0*FT_cur;
    if n > 20;  nWaves = 10 ; end
    if n > 200; nWaves = 100; end
    if n > 2000; nWaves = 1000; end
    
    for p = 1:nWaves*2
        [a,b] = find(ABS_FT_I == max(ABS_FT_I(:)), 1, 'first');
        ABS_FT_I(a,b) = 0;
        FT_cur(a,b) = FT_I(a,b);
    end
    I_cur = ifftn(FT_cur);

    canvas = cat(2, real(I - newIm - I_cur), zeros(size(I)), newIm);
    canvasShow = canvas;
    canvasShow(:,1:N(2)) = canvasShow(:,1:N(2)) + I_cur;
    subplot(2,1,2); imagesc(real(canvasShow), viewRange); axis equal tight; colormap gray;
    
    if n<20; 
        title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
    else
        title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);  
    end
    
    xlabel('residual image (left)             extracted wave (middle)             accumulated image (right)', 'FontSize',12); 
    pause((1/demoSpeed)*exp(-n/2));

    for L = [(min(max(1, round((1 + sin(linspace(-pi/2, pi/2, 40/min(20, n))))/2*N(2))), N(2)+1))   (N(2)+1)]
        canvasShow = canvas;
        canvasShow(:,L:(L+(N(2)-1))) = canvasShow(:,L:(L+(N(2)-1))) + I_cur;
        subplot(2,1,2); imagesc(real(canvasShow), viewRange); axis equal tight; colormap gray;
        if n<20; 
            title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
        else
            title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);  
        end
        xlabel('residual image (left)             extracted wave (middle)             accumulated image (right)', 'FontSize',12);
        pause((1/demoSpeed)*exp(-(n+4)));
        drawnow
    end
    if n<20; 
        title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
    else
        title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);  
    end
    xlabel('residual image (left)             extracted wave (middle)             accumulated image (right)', 'FontSize',12);
    pause((1/demoSpeed)*exp(-n/2));
    for L = N(2) + [(min(max(1, round((1 + sin(linspace(-pi/2, pi/2, 40/min(20, n))))/2*N(2))), N(2)+1))   (N(2)+1)]
        canvasShow = canvas;
        canvasShow(:,L:(L+(N(2)-1))) = canvasShow(:,L:(L+(N(2)-1))) + I_cur;
        subplot(2,1,2); imagesc(real(canvasShow), viewRange); axis equal tight; colormap gray;
        if n<20; 
            title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
        else
            title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);  
        end
        xlabel('residual image (left)             extracted wave (middle)             accumulated image (right)', 'FontSize',12);
        pause((1/demoSpeed)*exp(-(n+4)));
        drawnow
    end
    newIm = newIm + I_cur;
    n = n + nWaves;
end
