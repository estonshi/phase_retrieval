% Objective function for CS
% Unlinear and equal
function TV = obj(o)
size_o = size(o);
uy = o(2:size_o(1),1:size_o(2))-o(1:size_o(1)-1,1:size_o(2));
ux = o(1:size_o(1),2:size_o(2))-o(1:size_o(1),1:size_o(2)-1);
TV = sqrt(sum(sum(uy(:,1:size_o(2)-1).^2+ux(1:size_o(1)-1,:).^2)));
end