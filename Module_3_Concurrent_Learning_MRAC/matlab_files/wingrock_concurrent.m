%wingrock simulation with matched uncertainty
%Author: Girish Chowdhary
%distributed with adaptive control workshop in the 50th CDC

close all;
clear all;


%% sim params
t0=0;
tf=40;%170
dt=0.01;
t=t0:dt:tf;


%% initialization

x=[1.2 ;1];
n2=6;
Wstar=[0.8 0.2314 0.6918 -0.6245 0.0095 0.0214]';
Wstar_orig=Wstar;
W=Wstar*0+randn(1,6)'*0.0;

Kp = 1.5;                 % proportional gain
Kd = 1.3;                 % derivative gain

A=[0 1;-Kp -Kd];
B = [0; 1];
Q = eye(2);
P = lyap(A',Q);%A' instead of A because need to solve A'P+PA+Q=0
p=0;
[Lp,PP,E] = lqr(A,B,eye(2)*100,1);

%% adaptive control initialization and control parameters

gammaW=1;
gammaWbkk=1;%1.5;
gamWc=0; %set gamWc = 1 to use projection
lam_emod=0;%emod is off by default

nbg=0;
Xbg=Wstar*0;%the same size as sigma
add_new_points=1
pmax=30; %max dimension of history stack, once this is reached the singula value maximization algorithm triggers

%  load Xbg; %if starting with a previously recorded history stack, just
%  record into a mat file an existing Xbg matrix after a simulation run



%% commands
xref=0;
XREF=zeros(length(t),1);
 XREF(5/dt:7/dt)=0;
 XREF(15/dt:17/dt)=1;
 XREF(25/dt:27/dt)=-1;
%% reference model parameters
omegan_rm = 1;        % reference model natural freq (deg/s)
zeta_rm   = 0.5;        % reference model damping ratio
x_rm = x;
v_h=0;
v_ad=0;

v_crm=omegan_rm^2*(XREF(1)-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);




%% plotting array initialization 
index=1;

T_REC        = zeros(length(t),1);
X_REC        = zeros(length(t),1);
XDOT_REC     = zeros(length(t),1);
XRM_REC      = zeros(length(t),1);
XDOTRM_REC   = zeros(length(t),1);
XERR_REC     = zeros(length(t),1);
XDOTERR_REC  = zeros(length(t),1);
DELTACMD_REC = zeros(length(t),1);%control input
DELTAERR_REC = zeros(length(t),1);
W_REC = zeros(n2,length(t));
Wbdot_REC=W_REC*0;
VAD_REC = zeros(length(t),1);

%%
for t=t0:dt:tf
   


    e=x_rm-x;%compute reference model error
     xref=XREF(index);
      
    v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);
    v_pd=[Kp Kd]*e;
    deltaCmd = v_crm + v_pd - v_ad;%Nu
    delta = deltaCmd;
    v_h=deltaCmd-delta;
          
    [x,x_rm,xddot,deltaErr,v_crm]=wingrock_correct(x,x_rm,v_h,delta,dt,dt,Wstar,xref,omegan_rm,zeta_rm);
    x=x+randn(2,1)*0.000;
    sigma= [1;x(1);x(2);abs(x(1))*x(2);abs(x(2))*x(2);x(1)^3];
  
    
    Wd=-gammaW*e'*P*B*sigma-0.0*norm(e)*W;
    Wc=eye(n2)-Wd*pinv(Wd'*Wd)*Wd'*gamWc;%
    flag_point=0;
    %store data points
    oldrank=rank(Xbg*Xbg');
    old_eig=min(sqrt(eig(Xbg(:,1:nbg)*Xbg(:,1:nbg)')));
    old_det=det(Xbg*Xbg');
    if add_new_points==1
        if nbg==0
            nbg = nbg +1;flag_point=1;
            Xbg(:,nbg)=sigma;
            Timebg(nbg)=t;
            Delta_bg(:,nbg)=deltaErr;
            last_point=nbg;
        elseif norm(Xbg(:,last_point)-sigma)/norm(sigma)>=0.081 || rank([Xbg sigma],1e-5)>nbg% 
            nbg = nbg +1;flag_point=1;
            last_point=nbg;
            XbgOld=Xbg;
            Xbg(:,nbg)=sigma;
            Delta_bg(:,nbg)=deltaErr;
            Timebg(nbg)=index;
            if nbg>pmax
                [Xbg,Delta_bg,nbg,I]=data_point_remover(XbgOld,Delta_bg,sigma,old_eig,deltaErr,nbg);
                if norm(Xbg(:,1:nbg)-XbgOld)>0
                    Timebg(I)=index;
                    flag_point=1;
                else
                    flag_point=0;
                end
                
                last_point=nbg;
                
            end
         end
    end


%do concurrent learning
Wb=Wd*0;
summer=0;
  for kk=1:nbg
      if norm(Wc*Xbg(:,kk)-(Xbg(:,kk)))<1
         Wb=-gammaWbkk*(W'*Xbg(:,kk)-Delta_bg(kk))*Xbg(:,kk)+Wb;
      summer=Xbg(:,kk)*Xbg(:,kk)'+summer;
      end
  end
        wtilde=W-Wstar;
    
      tester(index)=wtilde'*Wc*summer*(eye(6)-Wc)*wtilde;
       vdot(index)=-e'*Q*e-wtilde'*(Wc*summer)*wtilde;
       lyap_v(index)=e'*P*e+wtilde'*inv(gammaW)*wtilde;


    Wdot=Wd+Wc*Wb-(norm(e))*W*lam_emod;

     W=W+dt*Wdot;
     v_ad=W'*sigma;
    
    
    %record for plotting
    T_REC(index)        = t;
    X_REC(index)        = x(1);
    XDOT_REC(index)     = x(2);
    XRM_REC(index)     = x_rm(1);
    XDOTRM_REC(index)   = x_rm(2);
    XERR_REC(index)     = e(1);
    XDOTERR_REC(index)  = e(2);
    DELTACMD_REC(index) = delta;%control input
    DELTAERR_REC(index) = deltaErr;
    W_REC(:,index)         =W;
    Wbdot_REC(:,index)        = Wb;
    WSTAR_REC(:,index)     =Wstar;
    VAD_REC(index)       =v_ad;
    FLAG_POINT(index)    = flag_point;
    MSV(index)          =old_eig;
    index = index+1;
end

%% plotting
figure(1);
subplot(2,1,1)
plot(T_REC,X_REC, T_REC,XRM_REC,'--');
xlabel('time (seconds)');
ylabel('deg');
title('roll angle');
legend('actual','ref model',0);
grid on;

subplot(2,1,2)
plot(T_REC,XDOT_REC, T_REC,XDOTRM_REC,'--');
xlabel('time (seconds)');
ylabel('xDot (deg/s)');
title('roll rate');
legend('actual','ref model',0);
grid on;

figure(2);
subplot(3,1,1);
plot(T_REC, XERR_REC);
xlabel('time (seconds)');
ylabel('xErr (deg)');
title('Position Error');
grid on;
subplot(3,1,2);
plot(T_REC, XDOTERR_REC);
xlabel('time (seconds)');
ylabel('xDotErr (deg/s)');
title('Angular Rate Error');
grid on;
subplot(3,1,3);
plot(T_REC, DELTACMD_REC);
xlabel('time (seconds)');
ylabel('\delta (deg)');
title('Angular Rate Error');
grid on;



figure(5)
plot(T_REC,W_REC(1,:),T_REC,WSTAR_REC(1,:),':',T_REC,W_REC(2,:),T_REC,W_REC(3,:),T_REC,W_REC(4,:),T_REC,W_REC(5,:),T_REC,W_REC(6,:))
hold on
plot(T_REC,WSTAR_REC(2,:),':',T_REC,WSTAR_REC(3,:),':',T_REC,WSTAR_REC(4,:),':',T_REC,WSTAR_REC(5,:),':',T_REC,WSTAR_REC(6,:),':')
legend('W(i)','W^*(i)',0)
xlabel('time (seconds)');
ylabel('W')


figure(12)
plot(X_REC,XDOT_REC)
grid on
xlabel('roll angle deg')
ylabel('roll rate deg/seconds')

figure(13)
plot(T_REC,DELTAERR_REC,T_REC,VAD_REC)
grid on
xlabel('Time')
ylabel('\Delta -\nu_{ad}')

figure(14)
plot(T_REC,MSV)
grid on
xlabel('Time seconds')
ylabel('\sigma_{min}(Z_k)')


figure(15)
%subplot(211)
plot(T_REC,lyap_v)
hold on
%subplot(212)
%plot(T_REC,FLAG_POINT)
grid on
xlabel('Time seconds')
ylabel('Lyapunov Function V')
%this code plots wher the data points were recorded: Warning: This is not
%optimized for speed, so it may take quite a while.
%for ii=1:max(size(T_REC))
%     if FLAG_POINT(ii)==1
%         figure(15)
%         plot(T_REC(ii),lyap_v(ii),'o')
%     end
%end
% nonzero=find(FLAG_POINT);

figure(16)
%subplot(211)
plot(T_REC,vdot)
hold on
%subplot(212)
%plot(T_REC,FLAG_POINT)
grid on
xlabel('Time seconds')
ylabel('Lie Derivative of Lyapunov Function')
nonzero=find(FLAG_POINT);




%% post processing
%save workspace
%save XBG Delta_bg Xbg last_point nbg