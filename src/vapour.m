function [Vq] = vapour(T,C,Plith)

Tv = 100 + 275/2.15e7 .* Plith;
Vq = (1 + tanh((T-Tv)./max(4,C.*400)))/2; 

end