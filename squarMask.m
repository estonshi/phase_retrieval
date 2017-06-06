function mask = squarMask( N,m,x,y )

mask = zeros(N,N);
mask(x-m/2+1:x+m/2,y-m/2+1:y+m/2) = 1;

end

