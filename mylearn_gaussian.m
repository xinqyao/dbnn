function CPD = mylearn_gaussian(CPD, fam, data, ns, cnodes)
%function CPD = learn_params(CPD, fam, data, ns, cnodes)
% This function is borrowed from original learn_params:
%
% LEARN_PARAMS Compute the maximum likelihood estimate of the params of a gaussian CPD given complete data
% CPD = learn_params(CPD, fam, data, ns, cnodes)
%
% data(i,m) is the value of node i in case m (can be cell array).
% We assume this node has a maximize_params method.
%
% Now it isn't hasty to update parameters

ncases = size(data, 2);
%CPD = reset_ess(CPD);
% make a fully observed joint distribution over the family
fmarginal.domain = fam;
fmarginal.T = 1;
fmarginal.mu = [];
fmarginal.Sigma = [];
if ~iscell(data)
  cases = num2cell(data);
else
  cases = data;
end
%hidden_bitv = zeros(1, max(fam));
for m=1:ncases
  % specify (as a bit vector) which elements in the family domain are hidden
  hidden_bitv = zeros(1, max(fmarginal.domain));
  ev = cases(:,m);
  hidden_bitv(find(isempty(ev)))=1;
  CPD = update_ess(CPD, fmarginal, ev, ns, cnodes, hidden_bitv);
end
%CPD = maximize_params(CPD);
