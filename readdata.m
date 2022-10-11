function dat=readdata(fp)
buff='';
c=fscanf(fp,'%c',1);
while(~feof(fp) && c~=char(10) && c~=',')
   buff=[buff c];
   c=fscanf(fp,'%c',1);
end
if(length(find(buff=='/'))==2)
   dat=inf;
else
   dat=str2num(buff);
end
if(isempty(dat)) 
    if(c==',') 
       dat=NaN; 
    elseif(~feof(fp))
       dat=readdata(fp);
    end
end
