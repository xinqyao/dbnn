function [sCPD,pred,predraw]=dbn(input,laa,lss,dmax,matfile)
%****************************************************************
%*                                                              *
%*     Protein Secondary Structure Prediction Based on DBN      * 
%*                                                              *
%*            written by Xin-Qiu Yao, copyright 2006            *
%*                                                              *
%****************************************************************
%
%[BNET PRED PREDRAW] = DBN(INPUT,LAA,LSS,DMAX,MATFILE)
%
%BNET     return the trained network
%PRED     return the prediction
%PREDRAW  return the raw probability 
%         distribution
%INPUT    input data:
%         1. the file name of the input 
%            data; for training
%         2. the data name in the matlab 
%            workspace; for prediciton
%LAA      the lag in profiles
%LSS      the lag in secondary structure
%DMAX     the maximum allowed segment length 
%         in the multinomial modelling
%MATFILE  file name for the network to be 
%         loaded; used for prediction only 
    
% Initialization %
pred=[];
predraw=[];
sCPD={};

SS=1; D=2; F=3;
d=[F+1:F+lss];
nd=length(d);
ind=[SS D F d]; 
if(laa==0), R=[]; else, R=max(ind)+1; end
nR=length(R);
ind=[ind R];
AA=max(ind)+1; 
ind=[ind AA];

nnds = length(ind);
bTrain=false;
%test=textread(input,'%d',1,'delimiter',',','emptyvalue',NaN);
if(ischar(input))
    bTrain=true;
end

%******************** Net Construction ********************
%disp('Constructing model...')
intra = zeros(nnds);
inter = zeros(nnds);

%intra connections
intra(SS,D)=1;
intra(D,F)=1;
intra(SS,AA)=1;
intra(d,AA)=1;   
intra(R,AA)=1;

%inter connections
inter(SS,SS)=1;
inter(D,D)=1;
inter(F,[SS D])=1;
if(nd>0)
   inter(SS,d(end))=1;
   for i=2:nd
      inter(d(i),d(i-1))=1;
   end
end

%node sizes
ns=zeros(1,nnds); %nodes sizes;
ns(SS)=3; ns(D)=dmax; ns(F)=2;
ns(d)=4;
ns(AA)=20;
ns(R)=laa*21;

dnodes=[SS D F d]; %discrete nodes
onodes=[1:nnds];   %observed nodes

% parameters typing
%% type F, d1(:), d2(1:nd-1), R(:), AA
eclass1 = [1 2 3 3+ones(1,nd)];
eclass1 = [eclass1 max(eclass1)+ones(1,nR)];
eclass1 = [eclass1 max(eclass1)+1];
eclass2 = [max(eclass1)+1 max(eclass1)+2 eclass1(F)]; 
eclass2 = [eclass2 max(eclass2)+ones(1,nd-1)];
eclass2 = [eclass2 max(eclass2)+ones(1,nd>0) eclass1(R) eclass1(AA)];
eclass  = [eclass1 eclass2];
nclasses= max(eclass);
%************************************************

bnet=mk_dbn(intra,inter,ns,'discrete',dnodes,...
            'observed',onodes,'eclass1',eclass1,...
            'eclass2',eclass2);

% CPD definition
bnet.CPD{eclass1(SS)}= tabular_CPD(bnet,SS);
bnet.CPD{eclass1(D)} = tabular_CPD(bnet,D);
bnet.CPD{eclass1(F)} = tabular_CPD(bnet,F);
if ~isempty(d)
    bnet.CPD{eclass1(d(1))} = tabular_CPD(bnet,d(1));
end
bnet.CPD{eclass1(AA)}= gaussian_CPD(bnet,AA);
if ~isempty(R)
    bnet.CPD{eclass1(R)}=root_CPD(bnet,R);
end
bnet.CPD{eclass2(SS)}=tabular_CPD(bnet,SS+nnds);
bnet.CPD{eclass2(D)}=tabular_CPD(bnet,D+nnds);
if ~isempty(d)
    bnet.CPD{eclass2(d(1))}=tabular_CPD(bnet,nnds+d(1));
end
if(nd>1)
   bnet.CPD{eclass2(d(end))}=tabular_CPD(bnet,nnds+d(end));
end

if(bTrain)
% temporary CPD stored
sCPD=cell(1,nclasses);
for i=1:nclasses
    sCPD{i}=struct(bnet.CPD{i});
end
%********************** END Net Construction ************************

%********************** Learning ************************************
disp('Learning...')
data=[];
fp=fopen(input,'r');
dat=readdata(fp); %(fp,'%f64',1,'delimiter',',','emptyvalue',NaN);
ncases=0;
while(~feof(fp))
%    dat=cell2num(dat);
    if(dat==inf)
        cases=mk_cases_dbn(data,bnet);
        [bnet,sCPD] = mylearn_params_dbn(bnet, sCPD, cases); %only statistics
        clear cases;
        data=[];
        ncases=ncases+1;
        if(mod(ncases,20)==0)
           disp(['Finished ' num2str(ncases)])
        end
    else
        data=[data;dat];
    end
    dat=readdata(fp); %textscan(fp,'%f64',1,'delimiter',',','emptyvalue',NaN);
end

for i=1:max(bnet.equiv_class(:))
    if adjustable_CPD(bnet.CPD{i})
        switch(class(bnet.CPD{i}))
            case 'tabular_CPD'
                sCPD{i} = update_CPD_tabular(sCPD{i});
%                j=find(eclass==i);
%                bnet.CPD{i}=tabular_CPD(bnet,j(1),sCPD{i}.CPT);
            case 'gaussian_CPD'
                bnet.CPD{i}=maximize_params(bnet.CPD{i});
                sCPD{i}=struct(bnet.CPD{i});
                sCPD{i}.Wsum=[];
                sCPD{i}.WXsum=[];
                sCPD{i}.WYsum=[];
                sCPD{i}.WXXsum=[];
                sCPD{i}.WYsum=[];
                sCPD{i}.WXYsum=[];
%                j=find(eclass==i);
%                bnet.CPD{i}=gaussian_CPD(bnet,j(1),'mean',sCPD{i}.mean,...
%                            'weights',sCPD{i}.weights,'cov',sCPD{i}.cov);
        end
    end
end

fclose(fp);
%clear cases;
%********** End learning *******************************************

end %bTrain 

if(~bTrain)
%********** Prediction ****************
%disp('Predicting...')
load(matfile);

% set CPD parameters 
bnet.CPD{eclass1(SS)}= tabular_CPD(bnet,SS,sCPD{eclass1(SS)}.CPT);
bnet.CPD{eclass1(D)} = tabular_CPD(bnet,D,sCPD{eclass1(D)}.CPT);
bnet.CPD{eclass1(F)} = tabular_CPD(bnet,F,sCPD{eclass1(F)}.CPT);
if ~isempty(d)
    bnet.CPD{eclass1(d(1))} = tabular_CPD(bnet,d(1),sCPD{eclass1(d(1))}.CPT);
end
bnet.CPD{eclass1(AA)}= gaussian_CPD(bnet,AA,'mean',sCPD{eclass1(AA)}.mean,...
                       'cov',sCPD{eclass1(AA)}.cov,'weights',sCPD{eclass1(AA)}.weights);
bnet.CPD{eclass2(SS)}=tabular_CPD(bnet,SS+nnds,sCPD{eclass2(SS)}.CPT);
bnet.CPD{eclass2(D)}=tabular_CPD(bnet,D+nnds,sCPD{eclass2(D)}.CPT);
if ~isempty(d)
    bnet.CPD{eclass2(d(1))}=tabular_CPD(bnet,nnds+d(1),sCPD{eclass2(d(1))}.CPT);
end
if(nd>1)
   bnet.CPD{eclass2(d(end))}=tabular_CPD(bnet,nnds+d(end),sCPD{eclass2(d(end))}.CPT);
end

ns=bnet.node_sizes_slice;
bnet.observed=[R AA];
bnet.hidden_bitv([SS D F d nnds+SS nnds+D nnds+F nnds+d])=1;
ev=mk_cases_dbn(input,bnet); %data,bnet);
engine=smoother_engine(jtree_2TBN_inf_engine(bnet));
engine=enter_evidence(engine,ev{1});
T=size(ev{1},2);
marg=zeros(ns(SS),T);
predd=zeros(1,T);
for j=1:T
   marg0=marginal_nodes(engine,SS,j);
   marg(:,j)=marg0.T;
   predd(j)=argmax(marg(:,j));
end
predraw=marg;
pred=map(predd,'str');

end %~bTrain;
