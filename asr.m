function g = asr( g,S,sup,alpha,mask )
%ASR Summary of this function goes here
%   Detailed explanation goes here
g = alpha.*rejectM(rejectSup(g,sup),S,mask) + (1-alpha).*g;   % ASR

end

