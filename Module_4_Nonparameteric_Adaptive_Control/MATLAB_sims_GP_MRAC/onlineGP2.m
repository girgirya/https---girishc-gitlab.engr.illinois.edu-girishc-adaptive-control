%============================== onlineGP ==================================
%
%  This code implements a sparse online GP algorithm by deleting the 
%  oldest element. 
%
%  This code is currently designed strictly for Gaussin kernels: if
%  you wish to extend it for non-Gaussian kernels, you MUST change
%  the values of k* away from unity!
%
%  Reference(s): 
%    None
% 
%  Inputs:
%    sigma  	 - bandwidth for the Gaussian kernel; either
%                  1 x 1 scalar or
%                  1 x d vector
%    noise      -  parameter for noise in GP; assumed given by oracle
%    m          -  the size of your budget
%    tol        -  tolerance for projection residual
%
%  Outputs:
%    see functions
%
%============================== onlineGP ==================================
%
%  Name:		onlineGP.m
%
%  Author: 		Hassan A. Kingravi
%
%  Created:  	2011/02/27
%  Modified: 	2012/02/29
%
%============================== onlineGP ==================================
function oGP = onlineGP2(sigma,noise,m,tol)
BV           = [];            % Basis vector set
K            = [];            % Kernel matrix
alpha        = [];            % mean parameter 
C            = [];            % inverted covariance matrix
Q            = [];            % inverted Gram matrix
current_size = [];
obs          = [];


%==[2] Setup the interface functions.
%
oGP.process = @process;
oGP.predict = @predict;
oGP.update  = @update;
oGP.get = @oGP_get;

  %------------------------------- process -------------------------------
  %
  %  Takes in a collection of data and generates an initial Gaussian
  %  process model with the associated kernel matrix, its inversion and
  %  the alpha vector. Currently, it's assumed that only a single point is
  %  passed in to initialize.
  %
  %  Inputs:
  %    data  	 - d x 1 data matrix passed in columnwise
  %    y         - 1 x 1 column vector of observations
  %
  %(
  function process(data,y)
    %create initial GP model 
    BV = data;
    obs = y';
    current_size = size(data,2);
    
    noise_x = noise + 1; % compute noise param
    Q = y/noise_x;
    C = -1/noise_x;    
    K = kernel(data,data,sigma);
    K = K + noise*eye(current_size);
    alpha = y/noise_x;
  end

  %------------------------------- predict -------------------------------
  %
  %  Given a new datapoint, predict a new value based on current model
  %
  %(
  function [f,var_x] = predict(x)
    k = kernel(x,BV,sigma)';
    f = k'*alpha;
    var_x = kernel(x,x,sigma) + k'*C*k;    
  end

  %------------------------------- update --------------------------------
  %
  %  Given a new data pair, update the model; remember, this is passed
  %  columnwise
  %(
  function update(x,y)      
    % first compute simple upate quantities
    k_t1 = kernel(x,BV,sigma)';   % pg 9, eqs 30-31   
    noise_x = noise + k_t1'*C*k_t1 + 1;
    q_t1 = (y - k_t1'*alpha)/(noise_x + noise);
    r_t1 = -1/(noise_x + noise);
    
    % compute residual projection update quantities 
    e_t1 = Q*k_t1; %residual vector pg 6, eq 16
    gamma_t1 = double(1-k_t1'*e_t1); %novelty of new point w.r.t RKHS: pg 7, eq 23
    eta_t1 = 1/(1+gamma_t1*r_t1); %numerical stability parameter    
    
    if gamma_t1 < tol
      % in this case, addition of point to basis doesn't help much, so 
      % don't add it, and compute update quantities in terms of old vectors
      % note that data, obs and gram matrix inverse not updated      
      s_t1 = C*k_t1 + e_t1;                  %pg 5, eqs 9, but modified
      alpha = alpha + q_t1*eta_t1*s_t1;
      C = C + r_t1*eta_t1*(s_t1*s_t1');                  
    else
      % in this case, you need to add the points
      current_size = current_size + 1;      
      
      %in this case, you can simply add the points    
      s_t1 = [C*k_t1; 1];    
      alpha = [alpha; 0] + q_t1*s_t1;
      C = [C zeros(current_size-1,1); zeros(1,current_size)] + r_t1*(s_t1*s_t1'); 
    
      % update basis vectors and observations
      BV = [BV x];
      obs = [obs; y];
    
      % update Gram matrix and inverse      
      K = [K k_t1; k_t1' 1]; 
      Q = inv(K);        
      
      if current_size <= m
        %do nothing   
      else
        % now you must delete one of the basis vectors; delete oldest
        index = 1;
%        current_size = current_size - 1;
%         BV(:,index) = [];
%         obs(index) = [];        
%         
%         %destroy matrices
%         alpha = [];
%         C     = [];
%         Q     = [];
%         K     = [];
%         
%         %recompute estimate
%         K = kernel(BV,BV,sigma) + 0.0001*eye(current_size);
%         Q = inv(K); % inverted gram matrix
%         K = K + noise*eye(current_size);
%         C = inv(K); % less numerically stable than cholesky    
%         alpha = C*obs;    
        
        
        %first compute scalar parameters
        a_s = alpha(index);
        c_s = C(index,index);
        q_s = Q(index,index);
        
        %compute vector parameters
        C_s = C(:,index);
        C_s(index) = [];
        Q_s = Q(:,index);
        Q_s(index) = [];

        %shrink matrices
        alpha(index) = [];
        C(:,index)   = [];
        C(index,:)   = [];
        Q(:,index)   = [];
        Q(index,:)   = [];
        K(:,index)   = [];
        K(index,:)   = [];
        
        %finally, compute updates
        alpha = alpha - (a_s/q_s)*(Q_s);
        C = C + (c_s/(q_s^2))*(Q_s*Q_s') - (1/q_s)*(Q_s*C_s' + C_s*Q_s');
        Q = Q - (1/q_s)*(Q_s*Q_s');
        
        current_size = current_size - 1;
        BV(:,index) = [];
        obs(index) = [];
      end          
     
    end       
  end

  %)
  %-------------------------------- get --------------------------------
  %ah
  %  Get a requested member variable.
  %
  %(
  function mval = oGP_get(mfield)

  switch(mfield)
    case {'basis','BV'}
	  mval = BV;
    case {'obs'}
	  mval = obs;      
	case {'K','kernel'}
	  mval = K;
	case {'Q'}
	  mval = Q;
	case {'current_size','size','current size'}
	  mval = current_size;      
  end

  end
  %)
  
end

%============================ Helper Functions ===========================

  %------------------------------- kernel ------------------------------
  % 
  %
  %
  function v =  kernel(x,y,sigma)

%  v = x'*y ;
  if(length(sigma) == 1) %same sigma
            d=x'*y;
            dx = sum(x.^2,1);
            dy = sum(y.^2,1);
            val = repmat(dx',1,length(dy)) + repmat(dy,length(dx),1) - 2*d;
            v = exp(-val./(2*sigma^2));
        else
            isigma = inv(diag(sigma.^2));
            d =  (x'*isigma)*y;
            dx = sum((x'*isigma)'.*x,1);
            dy = sum((y'*isigma)'.*y,1);
            val = repmat(dx',1,length(dy)) + repmat(dy,length(dx),1) - 2*d;
            v = exp(-val./2);
  end
  end

 %================================== kpca =================================

