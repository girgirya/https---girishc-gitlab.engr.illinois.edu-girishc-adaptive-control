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
global RESIDUAL_REC a bw n1 n2 n3 gammaV gammaW_Full B P back_flag Wc Vc index WdotStore VdotStore mom VAD_REC X_REC XDOT_REC qmod_win bv qmod_coeff delta_LM Kalr old_res gammaV_Full
global deltaErrS beta 
global e_mjwt omega RFT_rr X deltaErr theta


WdotStore=0;
VdotStore=0;

%% %% main sim control
CL_on=1; %set to 0 if no concurrent learning, to 1 to have concurrent learning

%% integration parameters
%integration parameters
t0=0;
tf=150
dt=0.05;
%% reference model parameters
omegan_rm = 1;        % reference model natural freq (rad/s)
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

gammaW=3;              %W learning rate %3
gammaV=1;               %V learning rate%1
kappa=0.0;                %NN damping
mom=0.00;               %momentum term
Kalr=0.000;              %Adaptive Loop Recovery Term%off time default
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
%% for separate integration
Wt = zeros( (n2+1), n3 );%NN hidden matrix 
Wb = zeros( (n2+1), n3 );%NN hidden matrix
Vt=V*0;Vb=V*0;
WW=W;
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


%% Background learning parameter initialization
max_back_points=5;
xback=[0; 0; 0;];                %background learning vector
eback=[0;0];
residual=zeros(1,max_back_points);
rresidual=0;
deltaErrS=0;
pp=1;                    %background learning data point index (pointer to last update)
no_back_points=0;        %no. of background learning points
back_point_selector=1;
zz=pp;
SIGMADASH=0;
temp=0;
back_flag=1;
beta=0;

gammaWb=ones(1,max_back_points)*gammaW*1/3; %learning rate for background learning W matrix adaptation
gammaVb=ones(1,max_back_points)*gammaV*3; %learning rate for background learning V matrix adaptation

if CL_on==1
    gammaW_Full=[gammaW, gammaWb];   %these will be sent to W V int
    gammaV_Full=[gammaV,gammaVb];
else
    gammaW_Full=[gammaW, gammaWb*0];   %these will be sent to W V int
    gammaV_Full=[gammaV,gammaVb*0];
end
delta_LM=gammaVb;
old_res=0;


%% Qmod parameters
qmod_win=10;
qmod_coeff=0.0;
%% Least squares parameter initialization
freqs = [0:.02:5]; % analysis frequencies (low,increment,high) [Hz]
omega = (freqs.*(2*pi))';
e_mjwt = ones(size(freqs))';
RFT_rr=0;
X=zeros(max(size(freqs)),n2+1);
theta=W*0;

%% Smoothing parameter initialization
ns=3;
wn=0.02;%strength of gaussian random noise
an=0.18;%frequency of sinusoidal noise
rps=860*2*pi/60;%860 rpm in radians per second
x2_hat=0;
xdotdot_hat=0;
x_hat=[0,0,0]';%[x,xdot,xdotdot]'
F=expm(dt*[0,1,0;0,0,1;0,0,0]);
H=[1 0 0;0 1 0];%y=[x,xdot]'
Qk=diag([0 0 wn^2]);
Rk=diag([1,1])*wn^2;
PD=diag([wn^2 wn^2 wn^2]);
pdm=F*PD*F'+Qk;
x_hat_min=x_hat;
xhat_hat=x_hat;
xhat_hatS=x_hat;
PknS=PD;
NN=30;
%One step storing variables
xhat_hatSS=x_hat;
PDS=PD;
x_hatS=x_hat;
x_hat_minS=x_hat;
%backward filter parameters
%% stability metric parameters
Twin=100;
lambda=0;

%% Data recording intialization


t=t0:dt:tf;
T_REC        = zeros(length(t),1);
X_REC        = zeros(length(t),1);
X_HAT_REC    = zeros(length(t),1);

XHAT_minus_STORE=zeros(length(t),ns)';
XHAT_STORE   = XHAT_minus_STORE;

XDOT_REC     = zeros(length(t),1);
XDOT_HAT_REC = zeros(length(t),1);
XDOTDOT_REC  = zeros(length(t),1);
XDOTDOT_HAT_REC  = zeros(length(t),1);
XHAT_HAT_REC    =  XHAT_minus_STORE;
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
COVARIANCE_STORE=zeros(ns,ns,length(t));
COVARIANCE_minus_STORE=zeros(ns,ns,length(t));
Kalman_Gain_Store=zeros(3,2,length(t));
BETA_store = zeros(length(t),1);
Lambda_store=zeros(length(t),1);

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
    xbar=[bv;x_hat(1);x_hat(2)];
    %% Sigma Update
    z=V'*xbar;%([n2*(n1+1)]*(n1+1)*1=(n2)*1
    [sigma,sigmadash]=sdash(z,a,bw,n2);
    
    v_ad = W'*sigma;             %% adaptation input
%%    
    %%W, V update
   [W,V,r]=W_V_int(gammaW_Full,gammaV_Full,sigma,sigmadash,V,W,xbar,e,P,B,kappa,residual,xback,no_back_points,controlDT);
  

%qMOD
    q=sigma*0;
    c=0;
       
    BETA_store(index)=beta;
%% Update control  
          v_pd=[Kp Kd]*e;
          deltaCmd = v_crm + v_pd - v_ad;%Nu
          lastControlUpdate = t;          
%% store data points for background learning
   
if abs(xbar(1))+abs(xbar(2))+abs(xbar(3))>bv; %store only data points with
                                             %nonzero xbar(2) and xbar(3)
                                             %values since xbar(1) = bv=1
        Ssize=size(xback);
        if (xbar-xback(:,pp))'*(xbar-xback(:,pp))/(xbar'*xbar)>=0.3;
            
            
            
            no_back_points=no_back_points+1;
            if no_back_points<=max_back_points
             pp=no_back_points;
            elseif no_back_points > max_back_points
               %  pp=1;
               if back_point_selector >max_back_points
                   back_point_selector=1;
               end
%              
               pp=back_point_selector;
               back_point_selector=back_point_selector+1;
               no_back_points=no_back_points-1;
             end
            xback(:,pp)=xbar
            
            %gammaW_Full(1+pp)=gammaWb(1);
            
            eback(:,pp)=e;%Don't need
            %residual(:,pp)=r;%Don't need
            rresidual(:,pp)=r; %Don't need in default method%v_ad; if method 2
            %deltaErrS(:,pp)=(v_ad-B'*(edot1-A*e));%v_ad-r;%deltaErr;%v_ad-r';%deltaErr;%%x(2)-deltaCmd;%v_ad-r';%%deltaErr;%stroe model error
            delta_i(:,pp)=deltaErr;%Not used in Algorithm, used for plotting
            back_flag(pp)=0;       %0 indicates a new point has been chosen, but not yet inducted 
            %disp('Selected Data point:')
            %disp(pp)
            counter(pp)=index;
            xhat_hatSS(:,pp)=F*x_hat;%XHAT_minus_STORE(:,smoother_address(s_counter));
            s_counter(pp)=0;
            BnS(:,:,pp)=eye(3);
            
            %some stuff for Gram Scmidt
            
                
            
        end %if (xbar-xback(:,pp))'*(xbar-xback(:,pp))/(xbar'*xbar)>=0.4*xbar;       
        
end % if abs(xbar(1))+abs(xbar(2))+abs(xbar(3))>1;


 old_res=residual;
 for ii=1:no_back_points
    %______________________________________________________________________
    %This FOR loop inducts and passes the background learning
    %information to the NN adaptation function,
    %______________________________________________________________________
    if back_flag(ii)==1
    z=V'*xback(:,ii);
    [sigma,sigmadash]=sdash(z,a,bw,n2);
    
    
    residual(:,ii)=W'*sigma-deltaErrS(:,ii);%delta_i(:,ii);%
    
        if norm(old_res(ii))>=norm(residual(ii))  %if improvement
            delta_LM(ii) = max(delta_LM(ii)/5,0.001); %bring closer to Newton
         else %
            delta_LM(ii) =min(delta_LM(ii)*5,100000); %bring closer to Grad
         end
          
     
    end
    %);

 end% end for
    
        delta = deltaCmd;
    
%% state and reference model Propogation
    [x,x_rm,xddot,deltaErr]= stateprop(x,x_rm,v_crm,v_h,delta,dt,controlDT);     
    
    %temp=xddot(2)*dt+temp;
    xn=x+an*sin(rps*index*dt)+randn(2,1)*wn;
   
    
%% Smoothing/Filtering
   
%Forward Kalman filter 
    Yk=xn;%[X_REC(index),XDOT_REC(index)]';
    x_hat_min=F*x_hat;%propogate state estimate
   pdm=F*PD*F'+Qk; %extrapolate covariance
   rino=H*pdm*H'+Rk;
        %kalman gain calculation
        Kk  =pdm*H'*inv(rino);
   
   x_hat=x_hat_min+Kk*(Yk-H*x_hat_min);%update state
   PD=(eye(3)-Kk*H)*pdm;%%update covariance
   %if backflag(pp)==0

   %______________Backwards Fixed point smoother_______________________
   %This algorithm runs parallel branches of smoothing instances for each
   %selected data points for a predefined smoothing horizon NN
   %A point is designated for smoothing if back_flag(.)=0
   for kk=1:no_back_points
   if back_flag(kk) == 1 
    
   elseif back_flag(kk)==0;%>=100+1 && index<=101+NN
        adress=counter(kk);
        Ai=PD*F'*inv(F*PD*F'+Qk);%COVARIANCE_STORE(:,:,adress+s_counter(kk))*F'*inv(F*COVARIANCE_STORE(:,:,adress+s_counter(kk))*F'+Qk);%COVARIANCE_minus_STORE(:,:,adress+counter(kk)+1));
        BnS(:,:,kk)=BnS(:,:,kk)*Ai;
        xhat_hatSS(:,kk)=xhat_hatSS(:,kk)+BnS(:,:,kk)*(x_hat-x_hat_min);%%(XHAT_STORE(:,adress+s_counter(kk))-XHAT_minus_STORE(:,adress+s_counter(kk)));%
        s_counter(kk)=s_counter(kk)+1;
       if s_counter(kk) == NN
              if back_flag(kk)==0                  %Data point has been marked but not inducted
              deltaErrS(:,kk)=xhat_hatSS(3,kk)-DELTACMD_REC(counter(kk));%XDOTDOT_REC(counter(ii))-DELTACMD_REC(counter(ii));%%DELTAERR_REC(counter(ii))%% model error assigned for data point
              back_flag(kk)=1;            %Data point is ready for induction
              %disp('Added Data point number:'),disp(kk)
              end
       end
   end
   end

%____________Backwards iterated fixed point smoother (fixed lag)___________    
   if index<=NN %wait untill fixed lag reached
       xhat_hatS=x_hat;
       PknS=PD;
   elseif index>NN
       xhat_hatS=XHAT_minus_STORE(:,index-NN);
       Bn=1;
   for jj=index-NN:index-2 %initiate fixed point smoother, smoother 
                           %estimates xhat(index-NN) using measurements
                           %up to index
    
       Ai=COVARIANCE_STORE(:,:,jj)*F'*inv(F*COVARIANCE_STORE(:,:,jj)*F'+Qk);%COVARIANCE_minus_STORE(:,:,jj+1));
       Bn=Bn*Ai;
    
    xhat_hatS=xhat_hatS+Bn*(XHAT_STORE(:,jj)-XHAT_minus_STORE(:,jj));%F*XHAT_STORE(:,jj-1));%
    PknS=PknS+Bn*(COVARIANCE_STORE(:,:,jj)-COVARIANCE_minus_STORE(:,:,jj))*Bn';
   end
   XHAT_HAT_REC(:,index-NN)=xhat_hatS;
   end
   
%% Stability metric analysis

if index>Twin+300
    for ii=1:Twin
        edot=A*[XERR_REC(index-ii);XDOTERR_REC(index-ii)]+B*(VAD_REC(index-ii)-DELTAERR_REC(index-ii));
        lambda=lambda+dt*log(norm(edot)/norm([XERR_REC(index-ii);XDOTERR_REC(index-ii)]+0.00000001));
    end
    Lambda_store(index)=lambda;
    lambda=0;
else
    Lambda_store(index)=0;

end
    
%% data recording   
    global store
    
    T_REC(index)        = t;
    XCMD_REC(index)     =xref;
    
    X_REC(index)        = x(1);
    X_HAT_REC(index)    =x_hat(1);
    XDOT_REC(index)     = x(2);  
    XDOT_HAT_REC(index)  = x_hat(2);%x2
    XDOTDOT_REC(index)  = xddot(2);            %x2dot=xdotdot
    XDOTDOT_HAT_REC(index)=x_hat(3);%xdotdot_hat; 
    COVARIANCE_STORE(:,:,index)=PD;
    COVARIANCE_minus_STORE(:,:,index)=pdm;
    XHAT_STORE(:,index) =x_hat;
    XHAT_minus_STORE(:,index)=x_hat_min;
    
    
    Kalman_Gain_Store(:,:,index)=Kk;
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
    RESIDUAL_REC(1:max(size(residual)),index)=residual(1,:)';
    EDOT_REC(:,index)   =edot1;
    STORE(index)        =temp;
    
    x_hat_minS=x_hat_minS;
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

figure(6)

plot(T_REC,RESIDUAL_REC(1,:),T_REC,RESIDUAL_REC(2,:),T_REC,RESIDUAL_REC(3,:),T_REC,RESIDUAL_REC(4,:),T_REC,RESIDUAL_REC(5,:),T_REC,RESIDUAL_REC(6,:),T_REC,RESIDUAL_REC(7,:))
ylabel('\nu_{ad}-estimated model error')
title('Difference betweeen stored estimate of model error and current estimate of model error')
xlabel('time (sec)')


figure(7)
subplot(3,1,1)
plot(T_REC,X_REC,'-',T_REC,X_HAT_REC,':',T_REC,XHAT_HAT_REC(1,:),'-.')
subplot(312)
plot(T_REC,XDOT_REC,'-',T_REC,XDOT_HAT_REC,':',T_REC,XHAT_HAT_REC(2,:),'-.')
subplot(313)
plot(T_REC,XDOTDOT_REC,'-',T_REC,XDOTDOT_HAT_REC,':',T_REC,XHAT_HAT_REC(3,:),'-.')
legend('Xdotdot','estimate of xdotdot','smoothed estimate of xdotdot',0)



figure(8);
subplot(311);
plot(T_REC, XERR_REC);
xlabel('time (sec)');
ylabel('xErr (rad)');
title('Position Error');
grid on;
subplot(312);
plot(T_REC, XDOTERR_REC);
xlabel('time (sec)');
ylabel('xDotErr (rad/s)');
title('Angular Rate Error');
grid on;
subplot(313);
plot(T_REC, Lambda_store);
xlabel('time (sec)');
ylabel('lambda');
title('stability metric');
grid on;

save workspace