function z = cached_naiveH(y)
%CACHED_NAIVEH z = cached_naiveH(y) calculates hilbert transform on real
%    axis (non-periodic version). Returns z = y + 1i*cached_naiveHc(y).
%
%    The returned value is the same form as Matlab built-in. z = y + im yi,
%    which is `hilbert` in the signal toolbox.
%    The 1/x/pi kernel and its FFT is cached by `cached_naiveHc`.
% 
%  arguments
%      y (1,:) {mustBeNumeric}
%  end

z = y + 1i*cached_naiveHc(y);

end
