% 3D phase retrieval
% clear all;
%------------------------------------------
% Note: input data should in form : {'data':pattern,'mask':mask}. 'mask'
%       is optional.
%------------------------------------------
%%
% Parameters. Carefully check before reconstruction
N = 250;    % Your martrix size N*N*N
mode = 'exp';  % 'simu' or 'exp'
data_path = '../data/exp_3D_5k.mat';
jmax = 1;    % update OSS filtering parameter while j changes
save_all = 0;
repeat_max = 5;

jjmax = 3;   % update iter, m, hio_factor and threshold while jj changes
iternum = [50,50,50];   % NOTE : real total iter numbers are iter*j
m_series = [40,40,40];    % support area constraint
init_hio_factor = [0.3,0.3,0.1];
init_threshold = [0.5,0.707,0.707];    % support intensity constraint

% initiate reconstruction
% init_model = load('init_model.mat');
% newg = init_model.exp_pat;
% newg = init_model.simu_pat;
newg = rand(N,N,N);

%%
repeat_times = 0;
if repeat_max~=0
    patt_ave = zeros(N,N,N);
end

for repeat_times=1:repeat_max
    close all;
    g = rand(N,N,N);
    is_OK = check(init_hio_factor,m_series,iternum,jjmax);
    m = m_series(1);
    m2 = m/2;
    Support = squarMask3D(N,m,floor(N/2),floor(N/2),floor(N/2));
    newSupport = Support;

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

    try
        mask = sample.mask;
    catch
        mask = 1;
    end

    if exist('newpattern')
        pattern = newpattern;
        mask = 1;
    end

    figure;imagesc(log(1+squeeze(pattern(floor(N/2),:,:))));axis square;title('The modulus of diffraction pattern (log)');

    %%
    figure;
    for jj=1:jjmax
        hio_factor = init_hio_factor(jj);
        threshold = init_threshold(jj);
        for j=1:jmax
        %     g = newg;
            al = double(N/4.0*(jmax-j+1)/jmax);
            Rf = 1e10;
            for i=1:iternum(jj)
                %========= HIO ==========
                g = hio(g,pattern,newSupport,hio_factor,mask);
                %========= Filter =========
                g = filtering(g,al,newSupport);
                %========= Disp =========
                if strcmp(mode,'simu')
                    subplot(1,2,1);
                    imshow(real(squeeze(g(floor(N/2),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1)))),'InitialMagnification',200);
                else
                    subplot(1,2,1);
                    gg = abs(g);
                    q_pattern = fftshift(fftn(gg(:,:,:)));
                    imagesc(log(1+abs(squeeze(q_pattern(floor(N/2),:,:)))));
                    title('q space');
                    subplot(1,2,2);
                    imagesc(log(1+abs(squeeze(gg(floor(N/2),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1))))));
                    title('real space');
                end
                suptitle(strcat('iternum : ',num2str(jj),'->',num2str(j),'->',num2str(i)));
                pause(0.01);
                %==============calculate new g================
                temp_g = zeros(N,N,N);
                temp_g(Support==1) = g(Support==1);
                if length(mask)==1
                    real_g = abs(fftshift(fftn(temp_g)));
                else
                    real_g = abs(fftshift(fftn(temp_g))).*mask;
                end
                real_S = pattern.*mask;
                score = sum(sum(sum(abs(real_g-real_S))))/sum(real_S(:));
                if score<Rf
                    Rf = score;
                    newg = abs(g);
                    disp(i);
                end
                g = abs(g);
                % =================update support================= %
                if threshold==0
                    newSupport = Support;
                else
                    newSupport = zeros(N,N,N);
                    thres_g = threshold*(median(g(:))+max(g(:)))/2.0;  % %%%%%%%%%%
                    newSupport(g>thres_g) = 1;
                    newSupport = newSupport.*Support;
                end
            end
        end
        m = m_series(min(jj+1,jjmax));
        m2 = m/2;
        Support = squarMask3D(N,m,floor(N/2),floor(N/2),floor(N/2));
    end

    if repeat_max~=0
        [cx,cy,cz] = ind2sub(size(newSupport),find(newSupport==1));
        cx = round(mean(cx)); cy = round(mean(cy)); cz = round(mean(cz));
        temp_g = newSupport.*g;
        temp_g = temp_g(cx-m2:cx+m2,cy-m2:cy+m2,cz-m2:cz+m2);
        patt_ave(floor(N/2)-m2:floor(N/2)+m2,floor(N/2)-m2:floor(N/2)+m2,floor(N/2)-m2:floor(N/2)+m2) = ...
            patt_ave(floor(N/2)-m2:floor(N/2)+m2,floor(N/2)-m2:floor(N/2)+m2,floor(N/2)-m2:floor(N/2)+m2) + temp_g;
    end
end
g = patt_ave / repeat_max;
%%
if length(size(g))==2
    temp=fftshift(fft2(g));
else
    temp=fftshift(fftn(g));
end
abstemp=abs(squeeze(temp(floor(N/2),:,:)));
rad_new=abstemp(floor(N/2),:);
rad=squeeze(pattern(floor(N/2),floor(N/2),:));
figure;plot(1:N,rad,'k-');hold on;plot(1:N,rad_new,'r-');
ylim([0 700]);

if save_all==1
    save matlab init_hio_factor iternum jmax jjmax mode newg pattern Rf;
end