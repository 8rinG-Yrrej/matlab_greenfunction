% function [I] = uniformSubindex(A,c)
% Get the indices of a uniformly spaced subarray of A,
% centered about c. c defaults to be center of A.
function [I] = uniformSubindex(A,c)

if nargin < 2
  c = round(length(A)/2);
end
if c > length(A) || c < 1
  warning('uniformSubindex(A,c): c beyond indices of A.');
  I = [];
  return
end

w = diff(A);
rg = abs(w-w(c)) > 10*eps(max(abs(A(1:end-1)),abs(A(2:end))));
r = find(rg(c:end));
l = find(rg(c-1:-1:1));

if isempty(r)
  ir = length(A);
else
  ir = r(1)+c-1;
end

if isempty(l)
  il = 1;
else
  il = c-l(1)+1;
end

I = il:ir;

end
