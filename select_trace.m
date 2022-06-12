%load /Users/fialko/w/dat/AllQfaults/Ridgecrest.utm;
%x=Ridgecrest;
load /Users/fialko/Downloads/ridgecrest_downsample.txt;
x=ridgecrest_downsample;
fid = fopen('ridge.utm','w+');
for i=1:length(x(:,1))
  if (x(i,1) ~= 370)
   fprintf(fid,'%f %f \n',x(i,1),x(i,2));
  else
   fprintf(fid,'nan nan \n');
  end
end
status = fclose(fid);
