function y=map(x,type)
%Y=MAP(X)
%Map amino acids or secondary structures to numbers

str=['H' 'E' 'C'];
aa=['A' 'R' 'N' 'D' 'C' 'Q' 'E' 'G' 'H' 'I' 'L' 'K' 'M' 'F' 'P' 'S' 'T' 'W' 'Y' 'V'];
sz=length(x);
bchar = false;
bnum  = false;
bstr  = false;
baa   = false;
if(strcmp(type,'str'))
   bstr = true;
elseif(strcmp(type,'aa'))
   baa  = true;
else
   error('Unknow type!');
end
if ischar(x)
   bchar  = true;
elseif isnumeric(x)
   bnum = true;
else
   error('Unknow input value!');
end
if bchar & baa
   data=upper(x);
   data(find(data=='B')) = 'D';
   data(find(data=='X')) = 'A';
   data(find(data=='Z')) = 'Q';
   for i=1:sz
      ind=find(aa==data(i));
      if ~isempty(ind)
         y(i)=ind;
      else
         error('Unknow amino acid letter!');
      end
   end
elseif bchar & bstr
   data=x;
   data(find(data=='-')) = 'C';
   data(find(data=='_')) = 'C';
   data(find(data=='?')) = 'C';
   data(find(data=='.')) = 'C';
   data = upper(data);
   for i=1:sz
      ind=find(str==data(i));
      if ~isempty(ind)
         y(i) = find(str==data(i));
      else
         error('Unknow secondary structure letter!');
      end
   end
elseif bnum & baa
   assert(min(x)>=1 & max(x)<=20);
   y = aa(x);
elseif bnum & bstr
   assert(min(x)>=1 & max(x)<=3);
   y = str(x);
end
