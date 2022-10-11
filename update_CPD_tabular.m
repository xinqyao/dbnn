function CPD = update_CPD_tabular(CPD)
% my function to calculate parameters from estimates

count=myreshape(CPD.counts,CPD.sizes);
switch CPD.prior_type
 case 'none', CPD.CPT = mk_stochastic(count); 
 case 'dirichlet', CPD.CPT = mk_stochastic(count + CPD.dirichlet); 
 otherwise, error(['unrecognized prior ' CPD.prior_type])
end
CPD.counts=zeros(prod(CPD.sizes),1); %clear counts
