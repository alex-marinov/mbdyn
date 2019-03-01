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
%%  Building the regression coefficients I 
%%  and their derivatives
%%  given the number of interpolant points
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [I,varargout] = f_regression_derivatives_kdtree(opoints, inpoints, tree, basis, n_points, phi_ord, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Inputs:
%%  - opoints:  the desired output point of the regression.
%%  - inpoints: the vectors of known point positions arranged as rows.
%%  - tree: pointer to kdtree structure built using kdtree function if
%%                  empty a new tree is built 
%%  - basis: the type of basis desired:
%%          1 --> first-order polynomial
%%          2 --> full second-order polynomial
%%          3 --> full third-order polynomial
%%          4 --> first-order polynomial + quadratic terms
%%          5 --> sine and cosine
%%  - n_points: number of points to be used in the local support for the interpolation.
%%  - phi_ord: the order of weight desired. For example, in 3 dimensions:
%%          0 --> (1 - r)^2+
%%          2 --> (1 - r)^4+ (4r + 1)
%%          4 --> (1 - r)^6+ (35/5 r^2 + 18/3 r + 1)
%%          6 --> (1 - r)^8+ (32 r^3 + 25 r^2 + 8 r + 1)
%%  - Weight: the weight of the trim points (Z)  optional 
%%    
%% Outputs:
%%  - I: the regression matrix
%%  - I_i: cell of sparse matrices for the derivatives of the regression coefficients on the i-th direction.
%%  - Index: the known points that are taken for the regression.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Total number of known points:
n_pt = size(inpoints,1) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total number of query points
n_qpt = size(opoints,1) ; 
% space dimension
dim = size(inpoints,2);
if (dim ~= size(opoints,2))
    error('Incopatible dimension between input points and output points');
end
    
if nargout > 1
    comp_derivatives = 1;
else
    comp_derivatives = 0;
end
if (nargin > 6)
    useW = 1;
    Weight = varargin{1};
else
    useW = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FINDING THE n_points NEAREST POINTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Index = zeros(n_qpt, n_points+1);
dist  = zeros(n_qpt, n_points+1);
if isempty(tree)
    tree = kd_buildtree(inpoints,0);
end
if (length(tree) > n_points)
    for i = 1 : n_qpt
        idx  = (kd_knn(tree,opoints(i,:),n_points+1,0))';
        Index(i,:) = idx;
        for j = 1 : n_points + 1
            dist(i,j) = sqrt(sum((opoints(i,:) - inpoints(Index(i,j),:)).^2));
        end
    end
elseif (length(tree) > n_points - 1) 
    for i = 1 : n_qpt
        idx  = (kd_knn(tree,opoints(i,:),n_points,0))';
        Index(i,1:length(idx)) = idx;
        Index(i,length(idx)+1:end) = idx(end);
        for j = 1 : n_points + 1
            dist(i,j) = sqrt(sum((opoints(i,:) - inpoints(Index(i,j),:)).^2));
        end
    end
else
    error('MLS Interpolation: There are not enought points reduce the number of points required');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structure for the sparse I matrix
r = zeros(n_qpt*n_points,1);
c = zeros(n_qpt*n_points,1);
h = zeros(n_qpt*n_points,1);
if comp_derivatives
    rd = cell(dim,1);
    cd = cell(dim,1);
    hd = cell(dim,1);
    for i = 1 : dim
        rd{i} = zeros(n_qpt*n_points,1);
        cd{i} = zeros(n_qpt*n_points,1);
        hd{i} = zeros(n_qpt*n_points,1);
    end
end
coef = 0;

p_dim = 0;
for i = 0 : basis
    p_dim = p_dim + nchoosek(dim + i - 1, i);
end
P = zeros(size(Index,2) - 1, p_dim);

for stp = 1 : n_qpt
    
%% BUILDING MATRIX P:
%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : size(Index,2) - 1
    
        P(i,:) = f_polyn_basis(inpoints(Index(stp,i),:),basis, opoints(stp,:)) ;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P = P';
%% BUILDING THE MATRIX FOR WEIGHTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  (abs(dist(stp,end) - dist(stp,end-1)) > 100*eps)
        r_max = 0.5 * (dist(stp,end) + dist(stp,end-1)) ;
    else
        r_max = 1.1 * dist(stp,end);
    end
    
    if useW
        phi = Weight(i)*f_weight_function_givendist(dist(stp,1:end-1),dim,r_max,phi_ord) ;
    else
        phi = f_weight_function_givendist(dist(stp,1:end-1),dim,r_max,phi_ord) ;
    end

    % scaling proposto da Levin
    %Scale = diag(sqrt( phi ./ (phi + 1) ));
    %Scale = eye(n_points);
    %P = Scale * P ;
    %Phi =   diag(phi + 1) ;
    Phi =   diag(phi) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BUILDING MATRICES A and b:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = P' * Phi * P ;
    b = P' * Phi ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FINDING THE REGRESSION VALUE FOR THE CHOSEN POINT:  f(point) = p'(point) * a(point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Building p:
    p_point = zeros(1, size(P, 2));
    p_point(1) = 1;

    tol = max(size(A)) * norm(A) * eps*10^0 ;
    rank_A=rank(A,tol) ;
    size_A=size(A,1) ;
    %svd(A) 
    if rank_A < size_A
       a = pinv(A,tol) * b ;  % it's the interpolation coefficients vector
    else
       a = A \ b ;            % it's the interpolation coefficients vector
    end
    h(coef+1:coef+n_points) = p_point * a ;
    r(coef+1:coef+n_points) = stp*ones(n_points,1);
    c(coef+1:coef+n_points) = Index(stp,1:n_points)';
    %if (sum(h(coef+1:coef+n_points)) < 0.999) || (sum(h(coef+1:coef+n_points)) > 1.001)
    %    disp('INTERPOLATION Warning:');        
    %    disp(['Interpolation point ', num2str(stp)]); 
    %    disp('the problem may be ill conditioned please try to change the interpolation parameters');
    %end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if comp_derivatives 
%% BUILDING REGRESSION DERIVATIVES: %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_i = zeros(dim, size(P, 2));
        phi_i = zeros(size(phi));
        Phi_i = zeros(size(Phi,1), size(Phi,2), dim);
        b_i = zeros(size(b,1), size(b,2), dim);
        gamma_i = zeros(size(a,1),size(a,2), dim);


%% Loop for each direction:
        for i = 1 : dim
    
%% BUILDING DERIVATIBES FOR THE BASE p_i:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p_i(i,:) = f_polyn_basis_derivatives(opoints(stp,:),i,basis, opoints(stp,:)) ;

    
%% BUILDING DERIVATIVES FOR THE WEIGHTS phi_i:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            adist = mean(dist(stp,1:n_points));
            for j = 1 : size(Index,2) - 1
                phi_i(j) = f_weight_function_derivatives(opoints(stp,:),inpoints(Index(stp,j),:),i,dist(stp,j),dim,r_max,adist,phi_ord) ;
            end

% figure;
% plot(vars(Index),phi,'*')
% grid on
% figure;
% plot(vars(Index),phi_i,'*')
% grid on

            Phi_i(:,:,i) = diag(phi_i);
    
%% BUILDING A_i and B_i:
%%%%%%%%%%%%%%%%%%%%%%%%
            b_i(:,:,i) = P' * Phi_i(:,:,i);
    
%% BUILDING DERIVATIVES FOR GAMMA gamma_i:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rank_A < size_A
                gamma_i(:,:,i) = pinv(A,tol) * b_i(:,:,i) ;  
            else
                gamma_i(:,:,i) = A \ b_i(:,:,i) ;      
            end
    
%% BUILDING DERIVATIVES FOR THE REGRESSION COEFFICIENTS I_i:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rd{i}(coef+1:coef+n_points) = stp*ones(n_points,1);
            cd{i}(coef+1:coef+n_points) = Index(stp,1:n_points)';
            hd{i}(coef+1:coef+n_points) = p_i(i,:) * a + p_point * ( gamma_i(:,:,i) - gamma_i(:,:,i)* P * a) ;
    
        end
    end
    coef = coef + n_points;
end
I = sparse(r, c, h, n_qpt, n_pt);
if comp_derivatives
    I_i = cell(dim,1);
    for i = 1 : dim
        I_i{i} = sparse(rd{i}, cd{i}, hd{i}, n_qpt, n_pt);
    end
    varargout{1} = I_i;
end
