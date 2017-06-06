function mask = circle_mask(N,m,x,y)

mOKflag = 1;
if mod(m,2)==1
    disp('m is not an oven');
    mOKflag = 0;
end

r = m/2;

xyOKflag = 1;
if (x>=r)&(x<=(N-r))&(y>=r)&(y<=(N-r))
%     disp('x y ok');
    xyOKflag = 1;
else
    disp('x y error. Area overflow.');
    xyOKflag = 0;
end

if (xyOKflag == 1)&(mOKflag == 1)
    [xx yy] = meshgrid(-N/2:N/2-1);
    z = sqrt(xx.^2 + yy.^2);
    clear xx yy xyOKflag;
    z = (z<r);      
    z = circshift(z,[x-N/2,y-N/2]);
else
    disp('PIEmask something is wrong!');
    z = [];
end

mask = double(z);