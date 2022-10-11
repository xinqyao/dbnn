function data=readplain(fid)
%

data=[];
dat=readdata(fid);
while(~feof(fid))
   if(dat==inf)
      return;
   else
      data=[data; dat];
   end
   dat=readdata(fid);
end
