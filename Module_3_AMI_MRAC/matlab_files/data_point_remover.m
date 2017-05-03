function [Xbg,Delta_bg,nbg,I,SaveFlag]=data_point_remover(XbgOld,Delta_bg,sigma,old_eig,deltaErr,nbg,FilterNoiseFlag)
A=XbgOld; %save initial


%if FilterNoiseFlag==0
%nbg=nbg-1;
%end
                for ii =1:nbg
                    A(:,ii)=sigma;%input new point
                    %dets(ii)= det(A*A');%calculate determinant
                    eigs(ii)=min(svd(A*A'));%*A'));%min(sqrt(eig(A(:,1:nbg)*A(:,1:nbg)')));%%calculate min eigenvalue
                    A=XbgOld;%recover old
                end
                %[Y,I]=min(eigs);%look for the minimum eigenvalue one (this is the point to throw)
                [Y,I]=max(eigs);%look for the maximum minimum eigenvalue one (this is the point to throw)
                
                if Y-old_eig>=0.00001%Y>old_eig%abs(Y-old_eig)>=0.01  %incorporate point only if it is better than what I have
                    A(:,I)=sigma;
                    Delta_bg(I)=deltaErr;
                    
                    Xbg=A; %send out new A
                    SaveFlag=1;
                else 
                   
                    Xbg=A;%send out old A
                    SaveFlag=0;
                end
       