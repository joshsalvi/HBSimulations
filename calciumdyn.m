function f=calciumdyn(t,C)
vpmca = 7.5e-7;
vmet = 2.8e-4;
kpmca = 6e-7;
Po=0.15;
Cm=50e-9;
f=(vmet*Po-374*vpmca.*C.^2./(C.^2+(5*kpmca)^2))*(1-C/Cm); 
end
