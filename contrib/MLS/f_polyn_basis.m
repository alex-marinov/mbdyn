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
%%  Choosing the polynomial basis p
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = f_polyn_basis(vars,type, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% the polinamials may be evaluated in shifted form (x - x_p) ecc. in order to 
%% have smaller numbers
%%
%% Inputs:
%%%%%  - nvars: the number of variables (i.e. the space dimension).
%%  - vars: the values of the variables for that point arranged as a row.
%%  - type: the type of basis desired:
%%          1 --> first-order polynomial
%%          2 --> full second-order polynomial
%%          3 --> full third-order polynomial
%%          4 --> first-order polynomial + quadratic terms
%% Outputs:
%%  - p: the polynomial basis evaluated at the given point.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(vars,1) >  1
    error('MLS interpolation error: Wrong point size in polinomial basis evaluation');
end

if nargin > 2
    shift = varargin{1};
else
    shift = zeros(1,size(vars,2));
end

 vars = vars - shift;
 
 switch type
     
     case 0
         
         p = [ 1 ] ;
         
     case 1
         
         p = [ 1 , vars ] ;
         
     case 2
         
         % first-order terms:
         p = [ 1 , vars ] ;
         
         % second-order terms:
         k = 1 ;
         for i = 1 : length(vars) 
             
             for j = i : length(vars) 
                 
                 s(1,k) = vars(i) * vars(j) ;
                 k = k + 1 ;
                 
             end
             
         end   
         
         p = [ p , s ] ;
         
     case 3
         
         % first-order terms:
         p = [ 1 , vars ] ;
         
         k = 1 ;
         m = 1 ;
         for i = 1 : length(vars) 
             
             for j = i : length(vars) 

                 % second-order terms:
                 s(1,k) = vars(i) * vars(j) ;
                 k = k + 1 ;
                 
                 % third-order terms:
                 for n = i : length(vars) 
                     t(1,m) = vars(i) * vars(j) * vars(n) ;
                     m = m + 1 ;
                 end
                 
             end
             
         end   
         
         p = [ p , s , t ] ;
         
     case 4
         
         % first-order terms:
         p = [ 1 , vars ] ;
         
         % quadratic terms:
         k = 1 ;
         for i = 1 : length(vars) 
                              
                 s(1,k) = vars(i) ^ 2 ;
                 k = k + 1 ;
                              
         end   
         
         p = [ p , s ] ;
         
 end
