% function [HMt,HMv,HMut,HMuv] = block_Htransf_mat(om,doms,iomUh,iomUt)
% Gives block hilbert transform matrices for a signal x = [t; u; v].
% with a uniform grid in the middle u part.
% HILBERT = [ HMt [HMut; @hilbert; HMuv] HMv ].
% x(iomUh:iomUt) == u.
function [HMt,HMv,HMut,HMuv] = block_Htransf_mat(om,doms,iomUh,iomUt)

Lom = length(om);
iomU = (iomUh:iomUt);
indt = (1:iomUh-1);
indv = (iomUt+1:Lom);
HMt = doms(indt)'./(om - om(indt)')/pi;
HMv = doms(indv)'./(om - om(indv)')/pi;

HMut = doms(iomU)'./(om(indt) - om(iomU)')/pi;
HMuv = doms(iomU)'./(om(indv) - om(iomU)')/pi;
HMt = zapdiag(HMt);
HMv = rot90(zapdiag(rot90(HMv,2)),2);

end

%%
function m = zapdiag(h)
m = h;
for ii = 1:min(size(h))
  m(ii,ii) = 0;
end
end
