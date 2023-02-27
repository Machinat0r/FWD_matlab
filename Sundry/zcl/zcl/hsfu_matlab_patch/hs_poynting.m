function J = hs_poynting(m)


it=0:.02:1;  it=it(:);

J=[ [it it it*0+1];[it*0+1 flipud(it) flipud(it)]];

end