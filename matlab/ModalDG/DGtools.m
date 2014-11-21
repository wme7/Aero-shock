classdef DGtools
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DGTOOLS class
    %   Build all one-dimensional Modal DG operators
    %
    %   by Manuel Diaz, NTU, 2013.10.12
    %   Copyright: all rights reserverd.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        kDeg
        nNodes
        xi
    end
    
    properties (Dependent = true, SetAccess = private)
        Vandermonde
        MassMatrix
        invMassMatrix
        CoefDiffMatrix
        legLeftEnd
        legRightEnd
    end
    
    methods (Static)
        function legP = legendreP(xi,l)
            % Construct Legendre Polynomials
            %**************************************************************
            % Compute the value of the legendre polinomial of degree 'kDeg'
            % for every column value 'xi'
            %**************************************************************
            x = sym('x'); 
            temp = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
            legP = subs(temp,xi);
        end
        
        function dlegP = dlegendreP(xi,l)
            % Construct derivatives of Legendre Polynomials
            %**************************************************************
            % Compute the derivative value of the legendre polinomial of
            % degree 'kDeg' for every column value 'xi'
            %**************************************************************
            x = sym('x'); 
            legP = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
            dlegP = subs(diff(legP,x),xi);
        end
               
        function h = DGnflux(f,df,u,strategy)
            %**************************************************************
            % General flux subroutine
            % We use 1 of 4 flux strategies/algorithms optimized for matlab.
            % coded by Manuel Diaz, 2012.12.21 (the end of the world)
            %**************************************************************
            
            % Parameters
            flx = f(u);     % flux value at every point of the domain
            dflx = df(u);   % flux slope at every point of the domain
            
            switch strategy
                case{1} % Roe Flux
                    
                    u_ave = (u(1,:) + u(2,:))/2;    % u @ cells boundaries
                    bool_p = df(u_ave) > 0; % if positive
                    bool_n = df(u_ave) <= 0;  % if negative
                    h = bool_p.*flx(1,:) + bool_n.*flx(2,:);
                    
                case{2} % LF
                    
                    alpha = max(max(abs(dflx)));
                    h = 0.5*(flx(1,:) + flx(2,:) - alpha*(u(2,:)- u(1,:)));
                    
                case{3} % LLF
                    
                    %fluxmat = [dflx(1:nx-1);dflx(2:nx)];
                    beta = max(abs(dflx));
                    h = 0.5*(flx(1,:) + flx(2,:) - beta.*(u(2,:)- u(1,:)));
                    
                case{4} % Upwind Flux
                    
                    % For dflux is constant along the domain!
                    a_p = max(max(dflx - abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a>0
                    a_n = max(max(dflx + abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a<0
                    h = a_p*flx(1,:) + a_n*flx(2,:);
                    
                otherwise
                    error('strategy not suported')
            end
        end
        
    end %Static Methods
     
    methods
        function obj = DGtools(localNodes) % Constructor
            obj.nNodes = length(localNodes);
            obj.kDeg = length(localNodes)-1;
            obj.xi = localNodes;
        end

        %%%%%%%%%%%%%%%%
        %  MAIN TOOLS: %
        %%%%%%%%%%%%%%%%
        
        function V = get.Vandermonde(obj) % Legendre Vandermonde Matrix
            % Construct Legendre Vandermonde Matrix, V,
            %**************************************************************
            % V is a matrix array in which each element is a Legendre
            % Polynomial of Degree n, n=0:kDeg evaluated at a column array
            % of points, xi.
            %
            % Coded by Manuel Diaz 2012.12.05
            %**************************************************************
            V = zeros(obj.kDeg);% Allocate
            for l = 0:obj.kDeg 	% All Polynomial Degrees up to kDeg
                j = l + 1;      % Dummy index
                for i = 1:obj.nNodes  % All Element Points in xi
                    V(i,j) = obj.legendreP(obj.xi(i),l);
                end
            end
        end
               
        %%%%%%%%%%%%%%%
        % MODAL TOOLS %
        %%%%%%%%%%%%%%%
        
        function legL = get.legLeftEnd(obj)
            % Scaled Legendre polynomials of degree 'l', l=0:kDeg,
            % evaluated at x = -1 
            %**************************************************************
            % same as evaluating legendreP(-1,l)
            %**************************************************************
            l = 0:obj.kDeg;
            legL = (-1).^(l); % LegP @ x_{i-1/2}^(+)
        end
    
        function legR = get.legRightEnd(obj)
            % Scaled Legendre polynomials of degree 'l', l=0:kDeg,
            % evaluated at x = +1 
            %**************************************************************
            % same as evaluating legendreP(+1,l)
            %**************************************************************
            l = 0:obj.kDeg;
            legR = (1).^(l);  % LegP @ x_{i+1/2}^(-)
        end
        
        function M = get.MassMatrix(obj)
            %**************************************************************
            % Mass Matrix
            %**************************************************************
            l = 0:obj.kDeg;
            M = diag(2./(2*l+1));
        end
        
        function invM = get.invMassMatrix(obj)
            %**************************************************************
            % Inverse of Mass Matrix
            %**************************************************************
            invM = inv(obj.MassMatrix);
        end
        
        function Dhat = get.CoefDiffMatrix(obj)
            %**************************************************************
            % Coefficients Differenciation Matrix
            %**************************************************************
            k = obj.kDeg; Dhat = zeros(k+1);
            for m = 0:k
                for j = m+1:2:k
                    Dhat(m+1,j+1) = (2*m+1);
                end
            end
        end
        
    end % Methods
end % Class

