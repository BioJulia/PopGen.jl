#= MATLAB Code from James Cai's PGEToolbox
%TAJIMA89D - calculates Tajima's D directly
%
% Syntax: [d,theta,pval] = tajima89d(nsam, Sn, sumk)
%
% nsam   - Number of smaples
% Sn     - Number of segregating sites or total number of mutations
% sumk   - Average number of nucleotide differences
% d      - tajima's D value
% theta  - estimated theta (4Neu) from Sn
% pval   - P-value for a double sided test assuming a normal
%          distribution for D

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (Sn==0), d=0; theta=0; pval=1; return; end

nx=1:(nsam-1);
a1 = sum(1./nx);
a2 = sum(1./nx.^2);
b1=(nsam+1)/(3*(nsam-1));
b2=2*(nsam*nsam+nsam+3)/(9*nsam*(nsam-1));
c1=b1-1/a1;
c2=b2-(nsam+2)/(a1*nsam)+a2/(a1^2);
e1=c1/a1;
e2=c2/(a1^2+a2);

d=(sumk-Sn/a1)/sqrt(e1*Sn+e2*Sn*(Sn-1));

if nargout>1
    theta=Sn/a1;
    % P-value for a double sided test assuming a normal distribution for D
    pval = 2*normcdf(-abs(d));
    % [above] This evaluates the tail probability directly, avoiding
    % numerical cancellation issues in the upper tail.
    % better than the following:
    % p = normcdf(z,0,1);
    % p = 2*min(p,1-p);
    % SEE: Matlab Solution ID - 1-1YDCNJ
end
=#

## the julia port ##