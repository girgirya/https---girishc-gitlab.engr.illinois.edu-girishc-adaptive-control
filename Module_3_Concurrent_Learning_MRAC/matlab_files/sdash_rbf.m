function sigma= sdash_rbf(x,rbf_c,rbf_mu,bw,n2)
%calculates sigma for a RBF NN given the state the rbf centers and rbf
%widths. Note that centers should be picked over a uniform distribution of
%the domain.
%Author:Girish Chowdhary
    sigma(1)=bw;
     for jj=2:n2+1
         sigma(jj,1)=exp(-norm(x-rbf_c(:,jj-1))^2/rbf_mu(jj-1)^2);
     end
     