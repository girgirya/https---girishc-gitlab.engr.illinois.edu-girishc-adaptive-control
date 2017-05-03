function [s,sp]=sdash(z,a,bw,n2)
%calculates sigma and sigmadash
%[sigma,sigmadash]=sdash(z,a,bw,n2)
%z=V'*xbar
%Author:Girish Chowdhary


ez  = exp(-a.*z);
ez(find(ez==Inf)) = 1e8;        % 'fix' Inf values to get proper sigma and sigma_prime
s   = 1./(1+ez);                %subtract 0.5 to remove bias in the sigmoids
sp  = a.*ez.*s.*s;
sp  = diag(sp);
s = [bw;s];                      %set 1 to 0 to remove the bias term
sp = [zeros(1,length(s)-1);sp];    
   
   