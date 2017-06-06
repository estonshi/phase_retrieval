function mask = squarMask3D( N,m,x,y,z )
%SQUARMASK3D 

mask = zeros(N,N,N);
mask(x-m/2+1:x+m/2,y-m/2+1:y+m/2,z-m/2+1:z+m/2) = 1;

end

