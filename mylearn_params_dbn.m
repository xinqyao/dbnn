function [bnet,sCPD] = mylearn_params_dbn(bnet, sCPD, data)

% bnet = mylearn_params_dbn(bnet, data)
% This module is borrowed from original BNT function:

% LEARN_PARAM_DBN Estimate params of a DBN for a fully observed model
% bnet = learn_params_dbn(bnet, data)
%
% data{l}(i,t) is the value of node i in slice t (can be a cell array) at series l
%
% We set bnet.CPD{i} to its ML/MAP estimate.

% Now it is capable of parameters typing.

nnds=bnet.nnodes_per_slice;
cnodes=bnet.cnodes_slice;
cnodes2=bnet.cnodes;
ns=bnet.node_sizes_slice;
ns2=bnet.node_sizes(:)';

ncases = size(data,2);

data2=[];
data1=[];
for i=1:ncases
    [ss T] = size(data{i});
    data1 = [data1 data{i}(:,1)];
    data2 = [data2 [data{i}(:,1:T-1);data{i}(:,2:T)]];
end

% slice 1
for j=1:nnds
    e=bnet.equiv_class(j,1);
    if adjustable_CPD(bnet.CPD{e})
        fam = family(bnet.dag,j);
        switch(class(bnet.CPD{e}))
            case 'tabular_CPD'
                sCPD{e} = mylearn_tabular(sCPD{e},fam,data1,ns,cnodes);
            case 'gaussian_CPD'
                bnet.CPD{e} = mylearn_gaussian(bnet.CPD{e},fam,data1,ns,cnodes);
        end
    end
end

       

% slices 2:T
% data2(:,t) contains [data(:,t-1); data(:,t)].
% Then we extract out the rows corresponding to the parents in 
% the current and previous slice.
% data2 = [data(:,1:T-1); 
%          data(:,2:T) ];
for j=1:nnds
  j2 = j+nnds;
  e=bnet.equiv_class(j,2);
  if adjustable_CPD(bnet.CPD{e})
    fam = family(bnet.dag,j2);
    switch(class(bnet.CPD{e}))
        case 'tabular_CPD'
            sCPD{e} = mylearn_tabular(sCPD{e},fam,data2,ns2,cnodes2);
        case 'gaussian_CPD'
            bnet.CPD{e} = mylearn_gaussian(bnet.CPD{e},fam,data2,ns2,cnodes2);
    end
  end
end
