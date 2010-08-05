%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% REGRESSION MODULE FOR RESPONSE SURFACE METHODS
%%
%% Implemented from the article:
%% "A CONSERVATIVE MESH-FREE APPROACH 
%% FOR FLUID-STRUCTURE INTERFACE PROBLEMS" 
%% by Quaranta, Masarati and Mantegazza.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Choosing the weight function
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phiO = f_weight_function_givendist(dist,dim,r_max,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Inputs:
%%%%%  - nvars: the number of variables (i.e. the space dimension).
%%  - dist: distance from the desired point to the interpolant known point.
%%  - dim: space dimension.
%%  - r_max: maximum radius outside of which the weight must be nul.
%%  - type: the order of weight desired:
%%          !! See Wenland !! for 3D
%%          0 --> (1 - r)^2+
%%          2 --> (1 - r)^4+ (4r + 1)
%%          4 --> (1 - r)^6+ (35/5 r^2 + 18/3 r + 1)
%%          6 --> (1 - r)^8+ (32 r^3 + 25 r^2 + 8 r + 1)
%%          !! See Wenland !!
%%          -1 --> a weight = 1 
%%          -2 --> Levin Weights for interpolation 
%%
%% Outputs:
%%  - p: the polynomial basis evaluated at the given point.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Computing radii:


r = dist / r_max ;
r = r + 1e-8 ;  % to avoid singularity


d = dim;                  % space dimension
k = type;                 % weight order
l = floor(d/2) + k + 1 ;  % length of the weight vector coefficients

%if (r <= 1) %|| (max(abs(point-vars))/r_max <= 1)

  if type < 0
  
     switch abs(type)
         case 1
            
             phiO = ones(size(r)) ;
         
         case 2
            d_ave = mean(dist);
            r_ave = (dist + 1.e-8*r_max) / d_ave;
            phiO = 1 ./ (exp(r_ave.^2) - 1) .* exp(-1 ./ (1 - r).^2);
            if any(isinf(phiO))
                phi = zeros(size(phiO));
                phi(isinf(phiO)) = 1;
                phiO = phi;
            end
         otherwise 
             error('Unknown weight function type');
     end
  else
    
    clear phi;
    for j = 0 : l
        phi(j+1) = (-1)^j * nchoosek(l,j) ;
    end
            
    if k > 0
        
            for s = 0 : k - 1
            
            phi_prev = phi ;
            clear phi ;
            
            phi(0+1) = sum( phi_prev ./ ([0:length(phi_prev)-1] + 2) ) ;
            % phi_prev ./ ([0:length(phi_prev)-1]+2)
            phi(1+1) = 0 ;
            
            for j = 2 : l + 2*s + 2
                phi(j+1) = - phi_prev(j-2+1) / j ;
            end
            
        end
            
    end
    
    j = 0:l+2*k ;
    phi = phi / phi(1) ; % u
    
    if length(r) > 1
        for i = 1 : length(r) 
            phiO(i) = phi * (r(i).^j)' ;
        end
    else
        phiO = phi * (r.^j)' ;
    end
%     phi4 = (1-r) 
%     phi4 = (1-r).^3 .* (3*r+1) 
%     phi4 = (1-r).^5 .* (8*r.^2 + 5*r + 1) 
    
  end
  
%else
    
%    phi = 0.0 ;
    
%end