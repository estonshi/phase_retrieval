%------------------------------------------
% Note: (for simulated data) If central data are missed, please use small
% support (30~50)
%------------------------------------------
N = 250; % size of the diffraction pattern image
m = 40; % size of the sample image
m2 = m/2;
alpha = 1.0;

%%
% support
sup = squarMask(N,m,floor(N/2),floor(N/2));
% sup = circle_mask(N+1,m,ceil(N/2),ceil(N/2));sup=sup(1:353,1:353);

Sfile = load('pattern.mat');
% S = double(Sfile.pat(51:300,51:300));  % if it is 'emc_cen_1.mat'

% add noise if necessary
%   S = abs(S + wgn(N,N,0.2*mean(S(:))));

% add mask if necessary
%   mask=1;
%   mask = zeros(N,N)+1;
%   mask(145:155,145:155)=0;
%   S=S.*mask;
  mask = Sfile.mask(51:300,51:300);
  S = S.*mask;

figure;
imagesc(log(1+S));
axis square;
title('The modulus of diffraction pattern (log)');
%%
itnum = 200; % iteration number
jmax = 5;
newg = rand(N,N);
figure;
for j=1:jmax
    g = newg;
    al = double(N/2*(jmax-j+1)/jmax);
    Rf = 1e10;
for i = 1:itnum
    %=================HIO========================
     g = hio(g,S,sup,0.5,mask);
     
    %=================DM========================
%        g = DM(g,S,sup,mask);
%   

    %=================ASR========================
%          g = asr(g,S,sup,0.1,mask);

    %==============filtering===================
    g = filtering(g,double(al),sup);
    
    %==============display the reconstruct sample image===================
    imshow(abs(g(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1))),'InitialMagnification',200);
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

% phase figure
phase = acos(real(temp)./abstemp);
third = find(real(temp)<0 & imag(temp)<0);
phase(third) = abs(2*3.1416-phase(third));
fourth = find(real(temp)>0 & imag(temp)<0);
phase(fourth) = 3.1415+phase(fourth);
figure;imagesc(log(1+phase));