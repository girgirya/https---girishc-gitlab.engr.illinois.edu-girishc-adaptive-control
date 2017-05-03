

%% W integrate RK4 based integration for W and V matrices
        function [W,V,r]=W_V_int(gammaW,gammaV,sigma,sigmadash,V,W,xbar,e,P,B,kappa,residual,xback,pp,dt)
      
        
        r=e'*P*B;
      
       
%% propogate W
       xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W,residual,xback,pp,e);
       
%1
 
    rk1=dt*xp;
    x1=W+rk1/2;
%2
    xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,x1,residual,xback,pp,e);
    rk2=dt*xp;
    x1=W+rk2/2;

%3
    xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,x1,residual,xback,pp,e);
    rk3=dt*xp;
    x1=W+2*rk3;

%4
    xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,x1,residual,xback,pp,e);
    rk4=dt*xp;
    W=W+(rk1+2*(rk2+rk3)+rk4)/6;
    
   

%% propogate V
    clear xp;
    xp=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W,residual,xback,pp,e);

     %dt=0.5*dt;
   

    rk1=dt*xp;
    x1=V+rk1/2;
  
   %2
    xp=stat_V(gammaV,sigmadash,x1,xbar,r,kappa,W,residual,xback,pp,e);
    

   rk2=dt*xp;
   x1=V+rk2/2;
  
   %3  
    xp=stat_V(gammaV,sigmadash,x1,xbar,r,kappa,W,residual,xback,pp,e);
    
    rk3=dt*xp;
    x1=V+2*rk3;
       
   %4
    xp=stat_V(gammaV,sigmadash,x1,xbar,r,kappa,W,residual,xback,pp,e);
    rk4=dt*xp;
   V=V+(rk1+2*(rk2+rk3)+rk4)/6;
    
%% W update model (c1)
    function [Wdot]=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W,residual,xback,pp,e)
         global a bw n1 n2 P B Wc Vc gammaW_Full gammaV back_flag sigma1 RESIDUAL_REC index WdotStore mom  XERR_REC XDOTERR_REC Kalr
          gammaW=gammaW_Full;
         sWDOT=W*0; 
        %Wc=(eye(n2+1)-(sigma*sigma'/(sigma'*sigma)));
         Wd=-((sigma-sigmadash*V'*xbar)*r'+kappa*norm(r)*W+Kalr*sigmadash*V'*(sigmadash*V')'*W)*gammaW(1); %in the ALR the V' term might not be required
         Wc=eye(n2+1)-Wd*pinv(Wd'*Wd)*Wd';
        %Vc=(eye(n1+1)-(gammaV(1)*xbar*xbar'*gammaV(1)/(xbar'*gammaV(1)*gammaV(1)*xbar)));
         Vd= -gammaV(1)*(xbar*r'*W'*sigmadash+kappa*norm(r)*V);
         Vc=eye(n1+1)-Vd*pinv(Vd'*Vd)*Vd';
  
        for ii=1:pp
         %if abs(xback(1,ii))+abs(xback(2,ii))+abs(xback(3,ii))>1;
              if back_flag(ii)==1;
              z=V'*xback(:,ii);
              %z=V'*Vc*xback(:,ii);
              
              [sigma1,sigmadash1]=sdash(z,a,bw,n2);
              %rr=W'*sigma1-residual(ii);%nuad(W,V,x(i))-deltaerr(i);%%%residual(pp);
              rr=residual(ii);%when using error dynamics
              
        sWDOT=((sigma1-sigmadash1*V'*xbar)*rr+kappa*norm(rr)*W)*gammaW(1+ii)+sWDOT;
               
             end
        end
        
          gammaW_Full=gammaW;
         Wdot=Wd-Wc*(sWDOT)+(mom*WdotStore);%different adaptation laws for inst. and bkk learning
         WdotStore=Wdot;%-WdotStore;
%% V update model (c1)
        function[Vdot]=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W,residual,xback,pp,e)
         %Vdot = -gammaV*( xbar*e'*P*B*W'*sigmadash + kappa*norm(e)*V );fae
         global a bw n2 n1 Wc K SsVDOT Vc back_flag VdotStore mom
          sVDOT=V*0;
          SsVDOT=sVDOT;
          for ii=1:pp
              %if abs(xback(1,ii))+abs(xback(2,ii))+abs(xback(3,ii))>1;
              if back_flag(ii)==1
               
              z=V'*xback(:,ii);   
              %z=V'*Vc*xback(:,ii);
              [sigma1,sigmadash1]=sdash(z,a,bw,n2);
              rr=residual(ii);%when using error dynamics
                            
              %%working one
              sVDOT=gammaV(1+ii)*(xback(:,ii)*rr*W'*sigmadash1+kappa*norm(rr)*V)+sVDOT;%designed to regulate rr

           
              end
              
          end
          
         

         Vdot=-gammaV(1)*(xbar*r'*W'*sigmadash+kappa*norm(r)*V)-Vc*(sVDOT)+(mom*VdotStore);
         VdotStore=Vdot;%-VdotStore;
          
       