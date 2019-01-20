function [p,chi2stat,df]=chi2_test_oe(observed,expected,df)
% [p,c]=chisquare_test_oe(observed,expected) counts
% observed is required and should be in matrix format of counts.
% Expected and df will be calculated automaticaly from observed if excluded
% from input
% Chi-square test, by hand

if ~exist('expected','var') || isempty(expected) %observed must be in matrix format
   n=sum(observed(:));
   rn=sum(observed,1);
   cn=sum(observed,2);
   expected=(cn*rn)/n;
end
oe=(observed-expected).^2 ./ expected;
chi2stat = sum(oe(:));
if ~exist('df','var') || isempty(df)
    [r,c]=size(observed);
    if r>1 && c>1
        df=(r-1)*(c-1);
    else
        df=length(observed)-1; %Assumes most simple possibility, but note that
    end
    %df = (R-1)(C-1), R= num rows, C= num cols in contingency table. 
end
p = 1 - chi2cdf(chi2stat,df);
