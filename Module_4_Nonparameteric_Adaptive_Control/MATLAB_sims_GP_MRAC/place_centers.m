%=========================== place_centers ================================
%
%  This code implements a way to place centers in \R^2 around a given
%  domain. 
%
%  Reference(s): 
%    None
% 
%  Inputs:
%    amin       -  minimum value of domain 
%    amax       -  maximum value of domain 
%    num        -  number of points to place 
%
%  Outputs:
%    rbf_c      -  2 x num centers
%
%=========================== place_centers ================================
%
%  Name:		place_centers.m
%
%  Author:      Girish Chowdhary
%  Modifier:    Hassan A. Kingravi
%
%  Created:  	2012/06/03
%  Modified: 	2012/06/03
%
%=========================== place_centers ================================
function rbf_c = place_centers(amin,amax,num)

domain_length = sqrt(num);
rbf_c = zeros(2,num);
current_y = amin;
increment = (amax-amin)/(domain_length-1);
current_point = 1;

for i=1:domain_length
  current_x = amin;  
  for j=1:domain_length
    rbf_c(:,current_point) = [current_x; current_y];      
    current_x = current_x + increment;
    current_point = current_point + 1;
  end  
  current_y = current_y + increment;
end

end
