function g = DM( g,S,sup,mask )
%DM Summary of this function goes here
%   Detailed explanation goes here
g = g + projectM(2*projectSup(g,sup)-g,S,mask) - projectSup(g,sup);

end

