function g = er( g,S,sup,mask )
%ER Summary of this function goes here
%   Detailed explanation goes here
g = projectM(projectSup(g,sup),S,mask);

end

