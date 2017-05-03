close all;
clear all;


%Neural network based control of an inverted pendulum
%
%Author Girish Chowdhary
%date: 18/oct/06

%Code simulates the control of an inverted pendulum with PCH and NN based
%model reference control,

%Distributed with Adaptive Control workshop at the 50th CDC
%% globals
global RESIDUAL_REC a bw n1 n2 n3 gammaV  B P  index  VAD_REC X_REC XDOT_REC  bv 



%% integration parameters
%integration parameters
t0=0;
tf=150
dt=0.05;
%% reference model parameters
omegan_rm = 0.5;        % reference model natural freq (rad/s)
zeta_rm   = 1;        % reference model damping ratio
%% command array

xcmd      = [0:dt:tf]*0;          % commanded position (rad)


stepup=dt;
step=20;
stepup1=step;
stepdown=2*step;

 xcmd(1/dt:10/dt)=0;
 xcmd(10/dt:30/dt)=1;
 xcmd(50/dt:70/dt)=1;
 xcmd(90/dt:110/dt)=1;
 xcmd(130/dt:150/dt)=1;


%% control parameters
%control parameters
controlUpdateDT = 0.02; % controller update rate
Kp = 2.1;                 % proportional gain
Kd = 2.2;                 % derivative gain
%% Neural network parameters
n1=2;                   %n1 no of states(inputs to NN)
n2=8;                   %n2 no of hidden nodes (neurons)
n3=1;                   %n3 no of output layer nodes\
amin = 0.01;
amax = 5;
a = zeros(n2,1);
% assign activation potential (different for all nodes, zero for the first)
for ii=1:n2,
   a(ii) = tan( atan(amin) + ( atan(amax) - atan(amin) )*(ii)/n2 );
end

gammaW=50;              %W learning rate %3
gammaV=20;               %V learning rate%1
kappa=+0.1;                %NN damping (e-mod)
bv=1;
bw=1;
deltaErr=0;

%% Initialization
% solve Lyapunov equation

  A=[0 1;-Kp -Kd];
  B = [0; 1];
  Q = eye(2);
  P = lyap(A',Q);%A' instead of A because need to solve A'P+PA+Q=0
  
  [Lp,PP,E] = lqr(A,B,eye(2)*100,1);
% initialize states
x1    = 0;
x2    = 0;
xddot  =[0; 0];

x1_rm = x1;             %x1 reference model
x2_rm = x2;             %x2 reference model
x = [x1;x2];
xn=x;
edot1=[0;0];
epsilon=0;
e1=edot1;
x_rm = x;
V = zeros( (n1+1), n2 );%NN input matrix
W = zeros( (n2+1), n3 );%NN hidden matrix

%%
% initialize control
v_crm    = 0;           %v reference model
v_pd     = 0;           %v proportional-derivative controller
v_ad     = 0;           %v adaptive (NN)
v_h      = 0;           %V PCH
delta    = 0;           %command initilization
deltaCmd = 0;           %PCh cmd
deltaLim = 0.1;         % control limit
delta_i=0;
deltaErrHat=0;
lastControlUpdate = 0;  %control update flag






%noise parameters
an=0.01;
rps=800;
wn=0.01;




%% Data recording intialization


t=t0:dt:tf;
T_REC        = zeros(length(t),1);
X_REC        = zeros(length(t),1);
X_HAT_REC    = zeros(length(t),1);


XDOT_REC     = zeros(length(t),1);
XDOTDOT_REC  = zeros(length(t),1);
XRM_REC      = zeros(length(t),1);
XDOTRM_REC   = zeros(length(t),1);
XERR_REC     = zeros(length(t),1);
XDOTERR_REC  = zeros(length(t),1);
DELTACMD_REC = zeros(length(t),1);
DELTA_REC    = zeros(length(t),1);
DELTAERR_REC = zeros(length(t),1);
DELTAERR_HAT_REC= zeros(length(t),1);
VCRM_REC     = zeros(length(t),1);
VAD_REC      = zeros(length(t),1);
VH_REC       = zeros(length(t),1);
VPD_REC      = zeros(length(t),1);
SIGMA_REC    = zeros(length(t),1);
V_REC        = zeros((n1+1), n2,length(t));
W_REC        = zeros((n2+1),length(t));
vv           = zeros(length(t),3);
STORE        = zeros(length(t),1);
RESIDUAL_REC = zeros(10,length(t));
EDOT_REC     = zeros(n1,length(t));


COUNTER_REC  = zeros(length(t),1);
%% Main simulation loop (convert to a function with sub functions)
index=1;
for t=t0:dt:tf
    
%% compute reference model error

e=x_rm-x;%compute reference model error
%% If control? controller runs at a different update rate: controlDT   
 xref =xcmd(index); 
     controlDT=t-lastControlUpdate;

        %update reference model
        v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);
        v_h=deltaCmd-delta;
 
         %update NN
    xbar=[bv;x(1);x(2)];
    %% Sigma Update
    z=V'*xbar;%([n2*(n1+1)]*(n1+1)*1=(n2)*1
    [sigma,sigmadash]=sdash(z,a,bw,n2);
    
    v_ad = W'*sigma;             %% adaptation input
%%    
    %%W, V update
   [W,V,r]=W_V_int_simple(gammaW,gammaV,sigma,sigmadash,V,W,xbar,e,P,B,kappa,controlDT);
  



%% Update control  
          v_pd=[Kp Kd]*e;
          deltaCmd = v_crm + v_pd - v_ad;%Nu
          lastControlUpdate = t;          





 % end for
    
        delta = deltaCmd;
    
%% state and reference model Propogation
    [x,x_rm,xddot,deltaErr]= stateprop(x,x_rm,v_crm,v_h,delta,dt,controlDT);     
    
    %temp=xddot(2)*dt+temp;
    xn=x+an*sin(rps*index*dt)+randn(2,1)*wn;
   
    

   
    
%% data recording   

    
    T_REC(index)        = t;
    XCMD_REC(index)     =xref;
    
    X_REC(index)        = x(1);
    XDOT_REC(index)     = x(2);  
    XDOTDOT_REC(index)  = xddot(2);            %x2dot=xdotdot 
  

    
    
    
    XRM_REC(index)      = x_rm(1);
    XDOTRM_REC(index)   = x_rm(2);
    
    XERR_REC(index)     = e(1);
    XDOTERR_REC(index)  = e(2);
    
    
    DELTACMD_REC(index) = deltaCmd;
    DELTA_REC(index)    = delta;
    DELTAERR_REC(index) = deltaErr;
    VCRM_REC(index)     = v_crm;
    VAD_REC(index)      = v_ad;
    VH_REC(index)       = v_h;
    VPD_REC(index)      = v_pd;
    SIGMA_REC(index)    = norm(sigma);
    V_REC(:,:,index)    = V;
    W_REC(:,index)      = W;
    EDOT_REC(:,index)   =edot1;
  

    index = index+1;
 end

%% Data post processing


%create array sequences of elements of V for plotting
for ii=1:tf/dt;
vv(ii+1,1)=V_REC(1,1,ii);
vv(ii+1,2)=V_REC(2,1,ii);
vv(ii+1,3)=V_REC(3,1,ii);
vv(ii+1,4)=V_REC(1,2,ii);
vv(ii+1,5)=V_REC(2,2,ii);
vv(ii+1,6)=V_REC(3,2,ii);
vv(ii+1,7)=V_REC(1,3,ii);
vv(ii+1,8)=V_REC(2,3,ii);
vv(ii+1,9)=V_REC(3,3,ii);
vv(ii+1,10)=V_REC(1,4,ii);
vv(ii+1,11)=V_REC(2,4,ii);
vv(ii+1,12)=V_REC(3,4,ii);
vv(ii+1,13)=V_REC(1,5,ii);
vv(ii+1,14)=V_REC(2,5,ii);
vv(ii+1,15)=V_REC(3,5,ii);
vv(ii+1,16)=V_REC(1,6,ii);
vv(ii+1,17)=V_REC(2,6,ii);
vv(ii+1,18)=V_REC(3,6,ii);
end
%% Plotting



figure(1);
subplot(2,1,1)
plot(T_REC,X_REC, T_REC,XRM_REC,'--', T_REC,XCMD_REC,':');
xlabel('time (sec)');
ylabel('x (rad)');
title('Position');
legend('actual','ref model','command',0);
grid on;

subplot(2,1,2)
plot(T_REC,XDOT_REC, T_REC,XDOTRM_REC,'--');
xlabel('time (sec)');
ylabel('xDot (rad/s)');
title('Angular Velocity');
legend('actual','ref model',0);
grid on;

figure(2);
subplot(2,1,1);
plot(T_REC, XERR_REC);
xlabel('time (sec)');
ylabel('xErr (rad)');
title('Position Error');
grid on;
subplot(2,1,2);
plot(T_REC, XDOTERR_REC);
xlabel('time (sec)');
ylabel('xDotErr (rad/s)');
title('Angular Rate Error');
grid on;

figure(3);
plot(T_REC,DELTA_REC, T_REC,DELTACMD_REC,'--');
xlabel('time (sec)');
ylabel('torque');
title('Control Torque');
grid on;
legend('actual','command',0);

figure(4);
plot(T_REC,DELTAERR_REC, T_REC,VAD_REC,'--');
xlabel('time (sec)');
ylabel('torque');
title('Model Error and NN Output');
legend('\Delta','\nu_{ad}',0);
grid on;



figure(5)
subplot(211)
plot(T_REC,W_REC(1,:),T_REC,W_REC(2,:),T_REC,W_REC(3,:),T_REC,W_REC(4,:),T_REC,W_REC(5,:),T_REC,W_REC(6,:),T_REC,W_REC(7,:),T_REC,W_REC(8,:),T_REC,W_REC(9,:))
xlabel('time (sec)');
ylabel('W')
subplot(212)
plot(T_REC,vv(:,1),T_REC,vv(:,2),T_REC,vv(:,3),T_REC,vv(:,4),T_REC,vv(:,5),T_REC,vv(:,6),T_REC,vv(:,7),T_REC,vv(:,8),T_REC,vv(:,9),T_REC,vv(:,10),T_REC,vv(:,11),T_REC,vv(:,12),T_REC,vv(:,13),T_REC,vv(:,14),T_REC,vv(:,15),T_REC,vv(:,16),T_REC,vv(:,17),T_REC,vv(:,18))
xlabel('time (sec)');
ylabel('V')



