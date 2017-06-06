function newsample = DownSample3( sample,newsize )
%UNTITLED3 
ori_size = size(sample);
center = ori_size/2;
temp1 = fftshift(fftn(sample));
lins = ceil(center-newsize/2):ceil(center-newsize/2)+newsize-1;
temp2 = temp1(lins,lins,lins);
newsample = abs(ifftn(ifftshift(temp2)));

end

