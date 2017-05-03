function [x,x_rm,xDot,deltaErr]=inv_pendu(x,x_rm,v_crm,v_h,delta,dt,controlDT,Wstar)

deltaErr = Wstar'*[sin(pi*x(1)); abs(x(2))*x(2); exp(-x(1)*x(2))]; %model nonlinearity
%% Reference model update
    
    xp=stat_refmodel(x_rm,v_crm,v_h);
  
   %1
   rk1=controlDT*xp;
   x1=x_rm+rk1;
   %2
   xp=stat_refmodel(x1,v_crm,v_h);
   rk2=controlDT*xp;
   x1=x_rm+rk2;
   %3
   xp=stat_refmodel(x1,v_crm,v_h);
   rk3=controlDT*xp;
   x1=x_rm+2*rk3;
   %4
   xp=stat_refmodel(x1,v_crm,v_h);
   rk4=controlDT*xp;
   
   x_rm=x_rm+(rk1+2.0*(rk2+rk3)+rk4)/6;
   
    
%x_rm=x_rm+xp*controlDT;

%% propogate state dynamics
clear xp
xp=state(x,delta,Wstar);
    xDot=xp;
   %1
   rk1=dt*xp;
   x1=x+rk1;
   %2
   xp=state(x1,delta,Wstar);
   rk2=dt*xp;
   x1=x+rk2;
   %3
   xp=state(x1,delta,Wstar);
   rk3=dt*xp;
   x1=x+2*rk3;
   %4
   xp=state(x1,delta,Wstar);
   rk4=dt*xp;
   
   x=x+(rk1+2.0*(rk2+rk3)+rk4)/6;
  

deltaErr = Wstar'*[sin(pi*x(1)); abs(x(2))*x(2); exp(-x(1)*x(2))]; %model nonlinearity

% End main function   
% Begin nested functions
%% State model
    function [xDot] = state(x,delta,Wstar)
        deltaErr =Wstar'*[sin(pi*x(1)); abs(x(2))*x(2); exp(-x(1)*x(2))]; %model nonlinearity
        x1Dot=x(2);
        x2Dot=delta+deltaErr;%delta-(deltaerr)
        xDot=[x1Dot;x2Dot];
        %Xdotdot=delta+sinx-mod(xdot)*xdot

%% reference model state dynamics (continuous := c1)
    function [x_dot_rm] = stat_refmodel(x_rm,v_crm,v_h)
             x1Dot_rm = x_rm(2);
             x2Dot_rm = v_crm - v_h;
             x_dot_rm=[x1Dot_rm;x2Dot_rm];
%             A=[0 1;-.1 0];
%             B=[0;1];
%             global Kp
%             x_dot_rm=(A-B*Kp)*x_rm+B*(v_crm-v_h);