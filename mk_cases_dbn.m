function cases=mk_cases_dbn(data,bnet)
%CASES=MK_CASES_DBN(DATA,BNET)
% Note: the content of 
% continuous-value node MUST be put behind

cases={[]};

nnds   = bnet.nnodes_per_slice;
cnodes = bnet.cnodes_slice;
dnodes = bnet.dnodes_slice;
nodesz = bnet.node_sizes_slice;
onodes = bnet.observed;
conodes= intersect(cnodes,onodes);
donodes= setdiff(onodes,conodes);
blocksz= length(dnodes)+sum(nodesz(cnodes));

T=length(data)/blocksz;
cases{1}=cell(nnds,T);
for j=[1:T]
    if ~isempty(donodes)
        cases{1}(donodes,j)=num2cell(data((j-1)*blocksz+donodes));
    end
    current=(j-1)*blocksz+length(dnodes);
    for k=conodes
        cases{1}(k,j)={data(current+1:current+nodesz(k))};
        current=current+nodesz(k);
    end   
end
