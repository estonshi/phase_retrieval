function g2 = projectM(g,S,mask)

if length(size(g))==3
    G  = fftshift(fftn(g));
else
    G  = fftshift(fft2(g));
end
G2 = S .* (G ./ abs(G));
if length(mask)~=1
    G2(mask==0)=G(mask==0);
end
if length(size(g))==3
    g2 = ifftn(ifftshift(G2));
else
    g2 = ifft2(ifftshift(G2));
end
end