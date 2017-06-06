function newg = filtering( g, alpha, sup )
%   use gaussian to filter g
s = size(g);
G = zeros(s);
newg = zeros(s);
center = s/2;
if length(s)==3
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                G(i,j,k) = exp(-0.5*((i-center(1))^2+(j-center(2))^2+(k-center(3))^2)/alpha^2);
            end
        end
    end
    temp1 = g.*(~sup);
    temp2 = fftshift(fftn(temp1));
    temp3 = temp2.*G;
    temp4 = ifftn(ifftshift(temp3));
    newg(sup==1) = g(sup==1);
    newg(sup~=1) = temp4(sup~=1);
else
    for i=1:s(1)
        for j=1:s(2)
            G(i,j) = exp(-0.5*((i-center(1))^2+(j-center(2))^2)/alpha^2);
        end
    end
    temp1 = g.*(~sup);
    temp2 = fftshift(fft2(temp1));
    temp3 = temp2.*G;
    temp4 = ifft2(ifftshift(temp3));
    newg(sup==1) = g(sup==1);
    newg(sup~=1) = temp4(sup~=1);
end

end

