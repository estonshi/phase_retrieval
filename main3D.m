% 3D phase retrieval
% clear all;

%%
% Parameters. Carefully check before reconstruction
N = 250; 
mode = 'exp';  % 'simu' or 'exp'
data_path = '../data/exp_3D.mat';
jmax = 5;
jjmax = 1;
iternum = [200];
m_series = [150];
init_hio_factor = 0.3;

init_model = load('init_model.mat');
newg = init_model.exp_pat;
% newg = init_model.simu_pat;
wantsave = 0;

%%
is_OK = check(m_series,iternum,jjmax);
m = m_series(1);
m2 = m/2;
Support = squarMask3D(N,m,floor(N/2),floor(N/2),floor(N/2));

sample = load(data_path);
if strcmp(mode,'simu')
    sampleData = DownSample3(sample.data,[100,100,100]);
    sampleData = MatMap(sampleData,0,1);
    Window = zeros(N,N,N);
    Window(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1)) = sampleData;
    Window = Window.*Support;
    pattern = abs(fftshift(fftn(Window)));
else
    pattern = double(sample.data);
end

mask = sample.mask;

if exist('newpattern')
    pattern = newpattern;
    mask = 1;
end

figure;imagesc(log(1+squeeze(pattern(floor(N/2),:,:))));axis square;title('The modulus of diffraction pattern (log)');

%%
figure;
for jj=1:jjmax
for j=1:jmax
    g = newg;
    al = double(N/2*(jmax-j+1)/jmax);
    Rf = 1e10;
    hio_factor = init_hio_factor;
    for i=1:iternum(jj)
        %========= HIO ==========
%         hio_factor = hio_factor*(iternum(jj)-i)/iternum(jj);
        g = hio(g,pattern,Support,hio_factor,mask);
        %========= Filter =========
        g = filtering(g,al,Support);
        %========= Disp =========
        if strcmp(mode,'simu')
            subplot(1,2,1);
            imshow(real(squeeze(g(floor(N/2),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1)))),'InitialMagnification',200);
        else
            subplot(1,2,1);
            q_pattern = fftshift(fftn(squeeze(g(floor(N/2),:,:))));
            imagesc(log(1+abs(q_pattern)));
            title('q space');
            subplot(1,2,2);
            imagesc(log(1+abs(squeeze(g(floor(N/2),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1))))));
            title('real space');
        end
        suptitle(strcat('iternum : ',num2str(jj),'->',num2str(j),'->',num2str(i)));
        pause(0.01);
        %==============calculate new g================
        temp_g = zeros(N,N,N);
        temp_g(Support==1) = g(Support==1);
        if length(mask)==1
            real_g = real(fftshift(fftn(temp_g)));
        else
            real_g = real(fftshift(fftn(temp_g))).*mask;
        end
        real_S = pattern;
%         gamma = mean(real_g(:))/mean(real_S(:));
        score = sum(sum(sum(abs(real_g-real_S))))/sum(real_S(:));
        if score<Rf
            Rf = score;
            newg = real(g);
            disp(i);
        end
        g = real(g);
    end
end
    m = m_series(min(jj+1,jjmax));
    m2 = m/2;
    Support = squarMask3D(N,m,floor(N/2),floor(N/2),floor(N/2));
end

%%
if length(size(newg))==2
    temp=fftshift(fft2(newg));
else
    temp=fftshift(fftn(newg));
end
abstemp=abs(squeeze(temp(floor(N/2),:,:)));
rad_new=abstemp(floor(N/2),:);
rad=squeeze(pattern(floor(N/2),floor(N/2),:));
figure;plot(1:N,rad,'k-');hold on;plot(1:N,rad_new,'r-');
ylim([0 700]);

if wantsave==1
    save matlab init_hio_factor iternum jmax jjmax mode newg pattern Rf;
end