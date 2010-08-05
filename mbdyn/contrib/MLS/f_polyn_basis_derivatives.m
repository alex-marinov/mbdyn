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
%%  Computing the derivatives the polynomial basis p
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p_i = f_polyn_basis_derivatives(x,n_der,type, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Inputs:
%%%%%  - nvars: the number of variables (i.e. the space dimension).
%%  - x: the components of the point where the derivatives is to be computed, arranged as a row.
%%  - n_der: the component of the derivative.
%%  - type: the type of basis desired:
%%          1 --> first-order polynomial
%%          2 --> full second-order polynomial
%%          3 --> full third-order polynomial
%%          4 --> first-order polynomial + quadratic terms
%%          5 --> sine and cosine
%%          6 --> first-order polynomial + sine and cosine
%%          7 --> full second-order polynomial + sine and cosine
%%          8 --> first-order polynomial + sine^2 and cosine^2
%%          9 --> first-order polynomial + sine and cosine, for which data
%%          is normalized to max = pi.
%%          10 --> first-order polynomial + sine^2 and cosine^2, for which data
%%          is normalized to max = pi.
%%
%% Outputs:
%%  - p: the polynomial basis evaluated at the given point.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(x,1) >  1
    error('MLS interpolation error: Wrong point size in derivatives of polinomial basis evaluation');
end


if nargin > 3
    shift = varargin{1};
else
    shift = zeros(1,size(x,2));
end

 x = x - shift;

 switch type
     
     case 0
         
         p = [ 0 ] ;
         
     case 1
         
         p = [ 0 , zeros(1,size(x,2)) ] ;
         p(n_der+1) = 1 ;
         
     case 2
         
         % first-order terms:
         p = [ 0 , zeros(1,size(x,2)) ] ;
         p(n_der+1) = 1 ;
         
         % second-order terms:
         k = 1 ;
         for i = 1 : length(x) 
             
             for j = i : length(x) 
                 
                 if ( i==n_der & j==n_der )
                     s(1,k) = x(i) + x(j) ;
                 elseif i==n_der
                     s(1,k) = x(j) ;
                 elseif j==n_der
                     s(1,k) = x(i) ;
                 else
                     s(1,k) = 0 ;
                 end                     
                     
                 k = k + 1 ;
                 
             end
             
         end   
         
         p = [ p , s ] ;
         
     case 3
         
         % first-order terms:
         p = [ 0 , zeros(1,size(x,2)) ] ;
         p(n_der+1) = 1 ;
         
         k = 1 ;
         m = 1 ;
         for i = 1 : length(x) 
             
             for j = i : length(x) 

                 % second-order terms:
                 if ( i==n_der & j==n_der )
                     s(1,k) = x(i) + x(j) ;
                 elseif i==n_der
                     s(1,k) = x(j) ;
                 elseif j==n_der
                     s(1,k) = x(i) ;
                 else
                     s(1,k) = 0 ;
                 end                     
                 k = k + 1 ;
                 
                 % third-order terms:
                 for n = i : length(x) 
                     
                     if ( i==n_der & j==n_der & n==n_der )
                         t(1,m) = x(j) * x(n) + x(i) * x(n) + x(i) * x(j) ;
                     elseif ( i==n_der & j==n_der )
                         t(1,m) = x(n) * ( x(j) + x(i) ) ;
                     elseif ( i==n_der & n==n_der )
                         t(1,m) = x(j) * ( x(n) + x(i) ) ;
                     elseif ( j==n_der & n==n_der )
                         t(1,m) = x(i) * ( x(n) + x(j) ) ;
                     elseif ( i==n_der )
                         t(1,m) = x(j) * x(n) ;
                     elseif ( j==n_der )
                         t(1,m) = x(i) * x(n) ;
                     elseif ( n==n_der )
                         t(1,m) = x(i) * x(j) ;
                     else 
                         t(1,m) = 0 ;
                     end
                         
                     m = m + 1 ;
                 end
                 
             end
             
         end   
         
         p = [ p , s , t ] ;
         
     case 4
         
         % first-order terms:
         p = [ 0 , zeros(1,size(x,2)) ] ;
         p(n_der+1) = 1 ;
         
         % quadratic terms:
         k = 1 ;
         for i = 1 : length(x) 
                 
             if i==n_der
                 s(1,k) = 2 * x(i) ;
             else
                 s(1,k) = 0 ;
             end
             
             k = k + 1 ;
                              
         end   
         
         p = [ p , s ] ;
 end      
 p_i = p ;