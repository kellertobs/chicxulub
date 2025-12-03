function [Vq] = vapour(T,C,Plith)

Tv = 100 + 275/2.15e7 .* Plith;
Vq = 1./(1+exp(-(T-Tv)./max(4,C.*500))); 

end