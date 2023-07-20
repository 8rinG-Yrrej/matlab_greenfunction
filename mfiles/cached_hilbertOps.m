function [hilbert_c,hilbert_Re,hilbert_Im,hilbert_aGR,hilbert_aGA]=cached_hilbertOps(om,gpu_side_hilbert)
%cached_hilbertOps returns some helper functions for non-periodic Hilbert
%transform useful for frequency domain Green's functions, over uniform or uniform-log sample (frequency) grid. In the
%uniform case, they are all based on `cached_naiveHc`, which is returned in
%`hilbert_c`. It performs the non-periodic version of Matlab signal toolbox
%function `hilbert`.
% hilbert_c(x) is like the Matlab `imag(hilbert(x))`
% hilbert_Re is like the Matlab `hilbert`
% hilbert_Im is like the Matlab `1i*(hilbert(x))`
% hilbert_aGR returns retarded Green's functions from a spectral function
% hilbert_aGA returns advanced Green's functions from a spectral function

if nargin == 1
  gpu_side_hilbert=false;
end

if isuniform(om)
  hilbert_c = @(x) cached_naiveHc(x);
  hilbert_Re = @(x) cached_naiveH(x);
  hilbert_Im = @(x) 1i*hilbert_Re(x);

  hilbert_aGR = @(a) cached_naiveHc(a)-1i*a;
  hilbert_aGA = @(a) cached_naiveHc(a)+1i*a;
else
  % Break up hilbert transform into uniform frequency region and
  % nonuniform parts
  % v = [t;u;v]
  % HM = [HMt | HMut; hilbert; HMuv; HMv]

  iomU=uniformSubindex(om);
  iomUh=iomU(1); iomUt=iomU(end);
  doms=gradient(om);

  [HMt,HMv,HMut,HMuv] = block_Htransf_mat(om,doms,iomUh,iomUt);

  if gpu_side_hilbert
    HMut = gA(HMut);
    HMuv = gA(HMuv);
    HMt = gA(HMt);
    HMv = gA(HMv);
    iomUh_g = gA(iomUh);
    iomUt_g = gA(iomUt);
    hilbert_Re = @(x) block_Htransf(x,iomUh_g,iomUt_g,HMt,HMv,HMut,HMuv);
    hilbert_Im = @(x) 1i*hilbert_Re(x);

    hilbert_c = @(a) block_Hctransf(a,iomUh_g,iomUt_g,HMt,HMv,HMut,HMuv);
    hilbert_aGR = @(a) block_Hctransf(a,iomUh_g,iomUt_g,HMt,HMv,HMut,HMuv)-1i*a;
    hilbert_aGA = @(a) block_Hctransf(a,iomUh_g,iomUt_g,HMt,HMv,HMut,HMuv)+1i*a;
  else
    hilbert_Re = @(x) block_Htransf(gather(x),iomUh,iomUt,HMt,HMv,HMut,HMuv);
    hilbert_Im = @(x) 1i*hilbert_Re(gather(x));

    hilbert_c = @(a) block_Hctransf(a,iomUh,iomUt,HMt,HMv,HMut,HMuv);
    hilbert_aGR = @(a) block_Hctransf(a,iomUh,iomUt,HMt,HMv,HMut,HMuv)-1i*a;
    hilbert_aGA = @(a) block_Hctransf(a,iomUh,iomUt,HMt,HMv,HMut,HMuv)+1i*a;
  end
end

end
