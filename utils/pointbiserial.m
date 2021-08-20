function [r,h,p,ci] = pointbiserial(d,x,alpha,tail)
%POINTBISERIAL Calculate Point biserial correlation 
%	R = POINTBISERIAL(X,Y) calculates the correlation with
%   D logical (boolean) variable, values 0 or 1. If values are numeric, all
%      values ~=0 are set to 0 and
%   X as continuos variable.
%
%   [R,H,P,CI] = POINTBISERIAL(D,X) 
%   H = 1 means that you can reject the null hypothesis, the means are
%   equal at the 5% significance level. P holds the p-value associated
%   with the t-statistic, CI is a 95% confidence interval for the true
%   difference in means.
%   A t-test is computed to check for the difference of both groups. Note
%   that a normal distribution is desired for the computation of the
%   t-test. For a non-parametric test use the 4-parameter call below.
%
%   [...] = POINTBISERIAL(D,X,ALPHA) can be used to set a different value for
%   alpha. Default is 5%.
%
%   [...] = POINTBISERIAL(D,X,ALPHA,TAIL) performs the test against the alternative
%   hypothesis specified by TAIL: (taken from ttest2)
%       'both'  -- "means are not equal" (two-tailed test)
%       'right' -- "mean of X with D == 1 is greater than mean of X with D == 0" (right-tailed test)
%       'left'  -- "mean of X with D == 1 is less than mean of X with D == 0" (left-tailed test)
%
%   [R,H,P,STATS] = POINTBISERIAL(D,X,ALPHA,'np') performs a non-parametric Wilcoxon rank sum
%   test (2-tailed). STATS contains test result struct (see ranksum).
%
% Example:
%   parm = rand(100,1);
%   gender = rand(100,1);
%   gender(gender<=.5)=0;
%	[r,h] = pointbiserial(gender,parm,.05);
% Uses:
%	Matlab Statistics Toolbox
%
% Frederik Nagel
% Institute of Music Physiology and Musicians' Medicine
% Hanover University of Music and Drama 
% Hannover
% Germany
%
% e-mail: frederik.nagel@hmt-hannover.de
% homepage: http://www.immm.hmt-hannover.de
%
% May 29, 2006.
%
% See also CORR, TTEST2, RANKSUM
error(nargchk(2, 4, nargin))
if(nargin==2)
    alpha=.05;
    tail = [];
elseif(nargin==3)
    tail = [];
end
% Convert numeric values to logicals
d = logical(d);
% Length 
n = length(d);
% Lengths of groups 0 and 1
n1 = sum(d);
n0 = sum(~d);
if(n0==0)
    error('There are no data with x=0!');
elseif(n1==0)
    error('There are no data with x=1!');
elseif (n0+n1 ~= n || sum(isnan(d))>0 || sum(isnan(x))>0)
    error('Data may not contain NANs');
end
% Mean of groups 0 and 1
x1 = mean(x(d));
x0 = mean(x(~d));
% SD of y
sx  = std(x);
%Correlation coefficient
r = (x1-x0)/sx*sqrt (n0*n1/n^2);
if(isempty(tail))
    [h,p,ci] = ttest2(x(d),x(~d));
elseif(strcmp(tail,'np'))
    [p,h,ci] = ranksum(x(d),x(~d),alpha);    
else
    [h,p,ci] = ttest2(x(d),x(~d),alpha,tail);    
end