% Constraint function for CS
% Unlinear and equal
function [C,Ceq] = cons(o)
global CS_mask;
global ori;
% calculate constraint
Ceq = abs(CS_mask.*fftshift(fft2(o))-ori);
C = -o;
end