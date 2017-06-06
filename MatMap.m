function Mout = MatMap(M,ymin,ymax)
%  scale M to [ymin ymax] 

data = M(:);
data = double(data);
mapdata = (ymax - ymin)*((data - min(data))/(max(data) - min(data)))...
    + ymin;

Mout = reshape(mapdata,size(M));