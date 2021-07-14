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
%%  Computing the derivatives of the weight function
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phi_i = f_weight_function_derivatives(x,x_known,n_der,dist,dim,r_max,d_ave,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Inputs:
%%%%%  - nvars: the number of variables (i.e. the space dimension).
%%  - x: the components of the point where the derivatives is to be computed, arranged as a row.
%%  - x_knwon: the components of the known points that influence the new point x.
%%  - n_der: the component of the derivative.
%%  - dist: distance from the desired point to the interpolant known point.
%%  - dim: space dimension.
%%  - r_max: maximum radius outside of which the weight must be nul.
%%  - type: the order of weight desired:
%%          !! See Henland !!
%%          0 --> (1 - r)^2+
%%          2 --> (1 - r)^4+ (4r + 1)
%%          4 --> (1 - r)^6+ (35/5 r^2 + 18/3 r + 1)
%%          6 --> (1 - r)^8+ (32 r^3 + 25 r^2 + 8 r + 1)
%%          !! See Henland !!
%%          -1 --> a weight = 1 
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

if (r <= 1) %|| (max(abs(point-vars))/r_max <= 1)

    if type < 0
  
        switch abs(type)
            case 1
             phi = 0.0 ;
         
         case 2
            r_ave = (dist + 1.e-8*r_max) / d_ave;
            phi = (-1 ./ (exp(r_ave.^2) - 1).^2) .* ( 2 * dist ./ d_ave^2  .* exp(r_ave.^2)) .* exp(-1 ./ (1 - r).^2) + ...
                1 ./ (exp(r_ave.^2) - 1) * (- 2/r_max ./ (1 - r).^3 *  exp(-1 ./ (1 - r).^2))
            if any(isnan(phi))
                phix = zeros(size(phi));
                phix(isinf(phi)) = 1;
                phi = phix;
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
    
    j = [0:l+2*k] ;
    phi = phi / phi(1) ;
    r = r + 1e-16 ;  % to avoid singularity
    
    phi = phi * ( j .* r.^(j-2) .* ( x_known(:,n_der)-x(n_der) )'/r_max )' ;
%     phi2 = ( -3 * (1-r).^2 .* (3*r+1) + 3 * (1-r).^3 ) .* sign( ( x_known(:,n_der)-x(n_der) )'/r_max ) 
%     phi2 = ( -5 * (1-r).^4 .* (8*r.^2 + 5*r + 1) + (1-r).^5 .* (16*r + 5) ) .* sign( ( x_known(:,n_der)-x(n_der) )'/r_max ) 
    
  end
  
else
    
    phi = 0.0 ;
    
end

phi_i = phi ;