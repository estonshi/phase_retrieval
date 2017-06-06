function g2 = rejectM(g,S,mask)

g2 = 2.*projectM(g,S,mask) - g;