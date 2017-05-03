function [Xbg,Delta_bg,nbg,I]=data_point_remover(XbgOld,Delta_bg,sigma,old_eig,deltaErr,nbg)
A=XbgOld; %save initial
%nbg=nbg-1;
                for ii =1:nbg
                    A(:,ii)=sigma;%input new point
                    %dets(ii)= det(A*A');%calculate determinant
                    eigs(ii)=min(svd(A*A',0));%min(sqrt(eig(A(:,1:nbg)*A(:,1:nbg)')));%%calculate min eigenvalue
                    A=XbgOld;%recover old
                end
                %[Y,I]=min(eigs);%look for the minimum eigenvalue one (this is the point to throw)
                [Y,I]=max(eigs);%look for the maximum minimum eigenvalue one (this is the point to throw)
               
                if Y>old_eig  %incorporate point only if it is better than what I have
                    A(:,I)=sigma;
                    Delta_bg(I)=deltaErr;
                    
                    Xbg=A; %send out new A
                else 
                   
                    Xbg=A;%send out old A
                end
       