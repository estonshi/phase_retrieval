% Compressive Sensitivity
global CS_mask;
global ori;
sim = 1;
% set up
if sim==1
    pic = imread('../lena.jpg');
    pic = pic(:,:,1);
    pic = imresize(pic,[64,64]);
    pic_fre = fftshift(fft2(pic));
    pic_fre(28:36,28:36) = 0;
    ori = pic_fre;
    CS_mask = zeros(64,64)+1;
    CS_mask(28:36,28:36) = 0;
    x0 = abs(ifft2(ifftshift(pic_fre)));
else
    Sfile = load('../emc_cen.mat');
    ori = double(Sfile.pat(51:300,51:300));
    CS_mask = Sfile.mask(51:300,51:300);
    x0 = newg;
end
% optimize
options=optimset('largescale','on','display','iter','MaxFunEvals',10000); 
[x,fval,exitflag,output,lambda] = fmincon(@obj,x0,[],[],[],[],[],[],@cons,options);
