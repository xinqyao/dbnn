function CPD = mylearn_tabular(CPD, fam, data, ns, cnodes)
%function CPD = mylearn_params(CPD, fam, data, ns, cnodes)
% This function is borrowed from original learn_params:
%
% LEARN_PARAMS Compute the ML/MAP estimate of the params of a tabular CPD given complete data
% CPD = learn_params(CPD, local_data)
%
% local_data(i,m) is the value of i'th family member in case m (can be cell array).
%
% Now it doesn't update parameters immediately.

local_data = data(fam, :); 
if iscell(local_data)
  local_data = cell2num(local_data);
end
assert(length(CPD.sizes) == size(local_data, 1));
P = prod(CPD.sizes);
indices = subv2ind(CPD.sizes, local_data'); % each row of data' is a case 
count = histc(indices, 1:P); %modified by yao
count = myreshape(count, size(CPD.counts));
CPD.counts = CPD.counts+count;
%switch CPD.prior_type
% case 'none', CPD.CPT = mk_stochastic(counts); 
% case 'dirichlet', CPD.CPT = mk_stochastic(counts + CPD.dirichlet); 
% otherwise, error(['unrecognized prior ' CPD.prior_type])
%end
