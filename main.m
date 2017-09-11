%------------------------------------------
% Note: (for simulated data) If central data are missed, please use small
% support (30~50)
%------------------------------------------
% Parameters. Carefully check before reconstruction
N = 401;    % Your martrix size N*N*N
mode = 'exp';  % 'simu' or 'exp'
data_path = '../data/shell.mat';
jmax = 3;    % update OSS filtering parameter while j changes
save_all = 0;
repeat_max = 1;

jjmax = 3;   % update iter, m, hio_factor and threshold while jj changes
iternum = [100,100,100];   % NOTE : real total iter numbers are iter*j
m_series = [50,50,50];    % support area constraint
init_hio_factor = [0.9,0.9,0.7];
init_threshold = [0.0,0.00,0.00];    % support intensity constraint

%newg = rand(N,N);
newg = zeros(N,N);
newg(195:205,195:205) = 1;

%%
repeat_times = 0;


for repeat_times=1:repeat_max
    close all;
    %
    g = zeros(N,N);
    g(195:205,195:205) = 1;
    %
    is_OK = check(init_hio_factor,m_series,iternum,jjmax);
    m = m_series(1);
    m2 = m/2;
    Support = squarMask(N,m,floor(N/2),floor(N/2));
    newSupport = Support;

    sample = load(data_path);
    if strcmp(mode,'simu')
%         sampleData = DownSample3(sample.data,[100,100,100]);
%         sampleData = MatMap(sampleData,0,1);
%         Window = zeros(N,N,N);
%         Window(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1)) = sampleData;
%         Window = Window.*Support;
%         pattern = abs(fftshift(fftn(Window)));
        assert(0~=0,'Do not support simulation yet.');
    else
        pattern = double(sample.data);
    end

    mask = sample.mask;

    if exist('newpattern')
        pattern = newpattern;
        mask = 1;
    end

    figure;imagesc(log(1+pattern));axis square;title('The modulus of diffraction pattern (log)');

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
                    imshow(real(squeeze(g(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1)))),'InitialMagnification',200);
                else
                    subplot(1,2,1);
                    gg = abs(g);
                    q_pattern = fftshift(fftn(gg(:,:,:)));
                    imagesc(log(1+abs(q_pattern)));
                    title('q space');
                    subplot(1,2,2);
                    imagesc(log(1+abs(squeeze(gg(floor(N/2)-(m2-1):floor(N/2)+1+(m2-1),floor(N/2)-(m2-1):floor(N/2)+1+(m2-1))))));
                    title('real space');
                end
                suptitle(strcat('iternum : ',num2str(jj),'->',num2str(j),'->',num2str(i)));
                pause(0.01);
                %==============calculate new g================
                temp_g = zeros(N,N);
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
                    newSupport = zeros(N,N);
                    thres_g = threshold*(median(g(:))+max(g(:)))/2.0;  % %%%%%%%%%%
                    newSupport(g>thres_g) = 1;
                    newSupport = newSupport.*Support;
                end
            end
        end
        m = m_series(min(jj+1,jjmax));
        m2 = m/2;
        Support = squarMask(N,m,floor(N/2),floor(N/2));
    end

    if repeat_max~=0
        [cx,cy] = ind2sub(size(newSupport),find(newSupport==1));
        cx = round(mean(cx)); cy = round(mean(cy));
        temp_g = newSupport.*g;
        temp_g = temp_g(cx-m2:cx+m2,cy-m2:cy+m2);
        patt_ave(floor(N/2)-m2:floor(N/2)+m2,floor(N/2)-m2:floor(N/2)+m2) = ...
            patt_ave(floor(N/2)-m2:floor(N/2)+m2,floor(N/2)-m2:floor(N/2)+m2) + temp_g;
    end
end
g = patt_ave / repeat_max;
%%
if length(size(g))==2
    temp=fftshift(fft2(g));
else
    temp=fftshift(fftn(g));
end
abstemp=abs(temp);
rad_new=abstemp(floor(N/2),:);
rad=squeeze(pattern(floor(N/2),:));
figure;plot(1:N,log(1+rad),'k-');hold on;plot(1:N,log(1+rad_new),'r-');
ylim([0 10]);

if save_all==1
    save matlab init_hio_factor iternum jmax jjmax mode newg pattern Rf;
end

% phase figure
phase = acos(real(temp)./abstemp);
third = find(real(temp)<0 & imag(temp)<0);
phase(third) = abs(2*3.1416-phase(third));
fourth = find(real(temp)>0 & imag(temp)<0);
phase(fourth) = 3.1415+phase(fourth);
figure;imagesc(log(1+phase));