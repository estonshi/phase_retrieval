function g = hio( g,S,sup,alpha,mask )
%HIO Summary of this function goes here
%   Detailed explanation goes here
g2 = projectM(g,S,mask);
g = ((g2>=0)&(sup==1)).*g2 + ((g2<0)|(sup==0)).*(g-alpha*g2);  % HIO
% g = (sup==1).*g2 + (sup==0).*(g-alpha*g2);
end

