clear all;clc;close all;

N = 250; % size of the diffraction pattern image
m = 30; % size of the sample image
m2 = m/2;
alpha = 1.0;

sample = rgb2gray(imread('lena.jpg')); 
sample = imresize(sample,[m m]);
sample = MatMap(sample,0,1); % the sample matrix is normalized

%figure;imshow(sample,'InitialMagnification',200);axis on;title('sample image');
%%
Ss = zeros(N,N);
% put the sample image into a zero background
Ss(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1) ,floor(N/2)-(m2-1):floor(N/2)+1+(m2-1)) = sample; 
% support
sup = squarMask(N,m,floor(N/2),floor(N/2));
% sup = circle_mask(N+1,m,ceil(N/2),ceil(N/2));sup=sup(1:353,1:353);

Ss = Ss.*sup;

% S = round(abs(fftshift(fft2(Ss)))); % generate the modulus of the diffraction pattern
Sfile = load('emc_cen.mat');
S = double(Sfile.pat(51:300,51:300));

% add noise if necessary
%   S = abs(S + wgn(N,N,0.2*mean(S(:))));

% add mask if necessary
%   mask=0;
%   mask = zeros(N,N)+1;
%   mask(145:155,145:155)=0;
%   S=S.*mask;
    mask = Sfile.mask(51:300,51:300);

figure;
imagesc(log(1+S));
axis square;
title('The modulus of diffraction pattern (log)');
%%
itnum = 2000; % iteration number
jmax = 10;
newg = rand(N,N);
figure;
for j=1:jmax
    g = newg;
    al = double(N/2*(jmax-j+1)/jmax);
    Rf = 1e10;
for i = 1:itnum
    %=================HIO========================
     g = hio(g,S,sup,0.3,mask);
     
    %=================DM========================
%        g = DM(g,S,sup,mask);
%   

    %=================ASR========================
%          g = asr(g,S,sup,0.1,mask);

    %==============filtering===================
    g = filtering(g,double(al),sup);
    
    %==============display the reconstruct sample image===================
    imshow(real(g(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1))),'InitialMagnification',200);
    %imshow(abs(g));
    title(strcat('iternum',num2str(i)));
    pause(0.01); % 
    %==============calculate new g================
    temp_g = zeros(N,N);
    temp_g(sup==1) = g(sup==1);
    if length(mask)==1
        real_g = abs(fftshift(fft2(temp_g)));
    else
        real_g = abs(fftshift(fft2(temp_g))).*mask;
    end
    real_S = S;
    gamma = mean(mean(real_g))/mean(mean(real_S));
    score = sum(sum(abs(real_g-gamma*real_S)))/sum(S(:));
    if score<Rf
        Rf = score;
        newg = real(g);
        disp(i);
    end
    g = real(g);
end
end
%=========================
temp=fftshift(fft2(newg));
abstemp=abs(temp);
rad_new=abstemp(N/2,:);
rad=S(N/2,:);
figure;plot(1:N,rad,'k-');hold on;plot(1:N,rad_new,'r-');