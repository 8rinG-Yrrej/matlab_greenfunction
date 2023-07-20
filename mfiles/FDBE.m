function [f,nB,fPnB,fMf,dfomfun,dnBomfun]=FDBE
%FDBE: collection of Fermi-Dirac and Einstein-Boltzmann distributions
% [f,nB,fPnB,fMf,dfomfun,dnBomfun]=FDBE;

f=str2func('@(xoT) (1-tanh(xoT/2))/2'); 
nB=str2func('@(xoT) if_then_else(xoT==0, 0, (coth(xoT/2)-1)/2)');

% n_f(x/T)+n_B(y/T)
fPnB = str2func('@(x,y) (if_then_else(y==0,0,coth(y/2))-tanh(x/2))/2');
% f(x/T)-f(y/T)
fMf = str2func('@(x,y) (tanh(y/2)-tanh(x/2))/2');

dfomfun=str2func('@(xoT,T) xoT/(4*T)./(cosh(xoT/2).^2)');
dnBomfun=str2func('@(xoT,T) if_then_else(xoT==0, 0, xoT/(4*T)./(sinh(xoT/2).^2))');

end
