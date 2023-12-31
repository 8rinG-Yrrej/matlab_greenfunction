function y = block_Hctransf(x,iomUh,iomUt,HMt,HMv,HMut,HMuv)
%BLOCK_HCTRANSF y = block_Hctransf(x,iomUh,iomUt,HMt,HMv,HMut,HMuv)
%   Gives block hilbert transform for a signal x = [t; u; v],
%   with a uniform grid in the middle u part.
%   HILBERT = [ HMt [HMut; @hilbert; HMuv] HMv ].
%   HMt,HMv,HMut,HMuv are generated by block_Htransf_mat
%   x(iomUh:iomUt) == u.
%   reshape to accomodate higher dimensional arrays

y = reshape((HMt*x(1:iomUh-1,:)) + (HMv*x(iomUt+1:end,:)) ...
  + [ (HMut*x(iomUh:iomUt,:))   ; ...
         cached_naiveHc(x(iomUh:iomUt,:))          ; ...
      (HMuv*x(iomUh:iomUt,:)) ; ...
    ], size(x));

end
