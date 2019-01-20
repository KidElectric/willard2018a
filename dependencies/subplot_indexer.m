function rci=subplot_indexer(r,c)
%  rci=subplot_indexer(r,c)
%  r= number of rows of subplots
%  c= number of columns of subplots
%  rci= rci(row,column) = subplot number to use
%       i.e. subplot(2,2,rci(1,2)) will plot in the first row, second
%       column subplot. 
%  BI
si=1:(r*c);
rci=zeros(r,c);
s=0;
for i=1:r
    rci(i,1:c)=si((1+s):(1+s+(c-1)));
    s=s+c;
end