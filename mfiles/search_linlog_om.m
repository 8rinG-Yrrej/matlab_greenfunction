function [om,LomU,Nom,omMax]=search_linlog_om(domU,omMaxU,omMax_goal,alpha)
%[om,LomU,Nom,omMax]=search_linlog_om(domU,omMaxU,omMax_goal,alpha)
% Returns a logarithmic om grid, with lowest om's linear and symmetric about
% zero. It does not contain zero. The uniform and log regions are connected
% smoothly so that up to 2nd order derivatives, the function om is
% continuous. |alpha| controls the 'radius' of the smooth corner that
% connects the two regions.
%
% Argin:
% - domU: minimal grid spacing
% - omMaxU: end of the uniform region
% - omMax_goal: max value of the grid desired
% - alpha: length scale for derviative matching to get a smoother grid.
% Argout:
% - om: frequency grid
% - LomU: length of the central uniform region of the grid
% - Nom: total length of the grid
% - omMax: max value of the grid

ll=0; rr=100*alpha;

while abs(rr-ll)>1
  mm=round((ll+rr)/2);
  maxom=(domU*alpha*(exp((mm)/alpha)-1)-domU/alpha*(mm).^2/2-domU/alpha^2*(mm).^3/6+omMaxU);
  if maxom>omMax_goal
    rr=mm;
  else
    ll=mm;
  end
end
Nom=mm;
omMax=maxom;
[om,LomU]=linlog_om(domU,mm,alpha,omMaxU);
end


function [om,LomU]=linlog_om(domU,Nom,alpha,omMaxU)
om = (domU/2:domU:omMaxU)';
LomU = length(om)*2;
om = unique([om; (domU*alpha*(exp((0:Nom)/alpha)-1)-domU/alpha*(0:Nom).^2/2-domU/alpha^2*(0:Nom).^3/6+om(end))']);
om = [flip(-om); om];
end

