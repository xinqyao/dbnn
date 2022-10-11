function writeraw(fid,data,s)
% write out data in FLAT format

if(nargin<3)
   for i=1:length(data)
      fprintf(fid,'%.3f,%.3f,%.3f\n',data(1,i),data(2,i),data(3,i));
   end
elseif(s=='r')
   for i=length(data):-1:1
      fprintf(fid,'%.3f,%.3f,%.3f\n',data(1,i),data(2,i),data(3,i));
   end
end
fprintf(fid,'//\n');
