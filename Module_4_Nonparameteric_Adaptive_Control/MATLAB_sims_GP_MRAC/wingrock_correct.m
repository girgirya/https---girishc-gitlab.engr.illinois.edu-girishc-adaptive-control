function [x,x_rm,xDot,deltaErr,v_crm]=wingrock_correct(x,x_rm,v_h,delta,dt,controlDT,Wstar,xref,omegan_rm,zeta_rm)
deltaErr=Wstar'*[1;x(1);x(2);abs(x(1))*x(2);abs(x(2))*x(2);x(1)^3];

%% Reference model update
    
    xp=stat_refmodel(x_rm,v_h,xref,omegan_rm,zeta_rm);
  
%    %1
%    rk1=controlDT*xp;
%    x1=x_rm+rk1;
%    %2
%    xp=stat_refmodel(x1,v_h,xref,omegan_rm,zeta_rm);
%    rk2=controlDT*xp;
%    x1=x_rm+rk2;
%    %3
%    xp=stat_refmodel(x1,v_h,xref,omegan_rm,zeta_rm);
%    rk3=controlDT*xp;
%    x1=x_rm+2*rk3;
%    %4
%    xp=stat_refmodel(x1,v_h,xref,omegan_rm,zeta_rm);
%    rk4=controlDT*xp;
%    
%    x_rm=x_rm+(rk1+2.0*(rk2+rk3)+rk4)/6;
   
x_rm=x_rm+xp*controlDT;

v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);

%% propogate state dynamics
clear xp
xp=state(x,delta,Wstar);
    xDot=xp;
   %1
%    rk1=dt*xp;
%    x1=x+rk1;
%    %2
%    xp=state(x1,delta,Wstar);
%    rk2=dt*xp;
%    x1=x+rk2;
%    %3
%    xp=state(x1,delta,Wstar);
%    rk3=dt*xp;
%    x1=x+2*rk3;
%    %4
%    xp=state(x1,delta,Wstar);
%    rk4=dt*xp;
%    
%    x=x+(rk1+2.0*(rk2+rk3)+rk4)/6;
  
x=x+dt*xp;
deltaErr=Wstar'*[1;x(1);x(2);abs(x(1))*x(2);abs(x(2))*x(2);x(1)^3]+randn(1,1)*0.01^2;;
% End main function   
% Begin nested functions
%% State model
    function [xDot] = state(x,delta,Wstar)
        x1Dot=x(2);
        deltaErr=Wstar'*[1;x(1);x(2);abs(x(1))*x(2);abs(x(2))*x(2);x(1)^3];

        x2Dot=delta+deltaErr;%delta-(deltaerr)
        xDot=[x1Dot;x2Dot];
        %Xdotdot=delta+sinx-mod(xdot)*xdot

%% reference model state dynamics (continuous := c1)
    function [x_dot_rm] = stat_refmodel(x_rm,v_h,xref,omegan_rm,zeta_rm)
             x1Dot_rm = x_rm(2);
             v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);
             x2Dot_rm = v_crm - v_h;
             x_dot_rm=[x1Dot_rm;x2Dot_rm];
%             A=[0 1;-.1 0];
%             B=[0;1];
%             global Kp
%             x_dot_rm=(A-B*Kp)*x_rm+B*(v_crm-v_h);