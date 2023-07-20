function z = cached_naiveHc(y)
%CACHED_NAIVEHC z = cached_naiveHc(y) is a naive implementation of Hilbert
%   transform (non-periodic version), using convolution kernel 1/pi/n,
%   1-n<=n<=n-1. At n=0 the kernel is zero. This is a smoother version of
%   hilbert transform with no even-odd effect or oscillatory behavior. Note
%   that naiveH^2 does not give back original function. This is essentially
%   the non-periodic version of Matlab signal toolbox `imag(hilbert(x))`.
%
% arguments
%     y (1,:) {mustBeNumeric}
% end
persistent ker fker

fftThres = 2000;
has_bigger = false;
if ismatrix(y) && ~isvector(y)
  has_mat = true;
  n = size(y,1);
elseif isvector(y)
  has_mat = false;
  n = length(y);
else
  has_bigger = true;
  n = size(y,1);
end

if isempty(ker) || 2*n ~= length(ker)+1
  x = (1-n):(n-1); % use a kernel with zero at the middle

  ker = 1./x/pi;
  ker(n) = 0;
  fker = fft(ker);
end

if has_bigger
  if ~iscolumn(ker)
    ker = ker(:);
    fker = fker(:);
  end
  % get the hilbert transform of y
  z = ifft(fft(y(:,:),2*n-1,1).*fker,[],1);
  z = reshape(z(n:end,:), size(y));

elseif has_mat

  if ~iscolumn(ker)
    ker = ker(:);
    fker = fker(:);
  end
  % get the hilbert transform of y
  if n > fftThres % using FFT
    z = ifft(fft(y,2*n-1).*fker);
    z = z(n:end,:);
  else
    z = zeros(size(y));
    for ii = 1:size(y,2)
      z(:,ii) = conv(y(:,ii),ker,'same');
    end
  end

else
  if isrow(fker) ~= isrow(y)
    ker = ker.';
    fker = fker.';
  end

  % get the hilbert transform of y
  if n > fftThres % using FFT
    z = ifft(fft(y,2*n-1).*fker);
    z = z(n:end);
  else % using plain convolution on vector
    z = conv(y,ker,'same');
  end

end

end
