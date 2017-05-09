function[eigvals, eigvecs] = quickSVD(C, varargin)
%% Run singular value decomposition on a matrix. Returns eigenvalues, eigenvectors with majority positive loadings.
%
% [eigvals, eigvecs] = quickSVD(C)
%   performs a full SVD on a matrix.
%
% [...] = quickSVD(C, 'svd', 'econ')
% Uses an economy sized svd decomposition rather than the full svd.
%
% [...] = qucikSVD(C, 'econ')
% Also performs the economy sized svd decomposition.
%
% [...] = quickSVD(C, 'svds', nEigs)
% Uses the svds decomposition and determines the first nEigs eigenvalues.
%
% [...] = quickSVD(C, 'svd')
% Performs the normal svd.
%
% ----- Inputs -----
%
% C: A matrix. May not contain NaN.
%
%
% ----- Outputs -----
%
% eigvals: A vector containing the eigenvalues of C.
%
% eigvecs: A matrix containing the eigenvectors of C. Each column is one
%       eigenvector.

% Parse Inputs
[svdFunc, svdArg] = parseInputs( varargin{:});

% Error check
errCheck(C);

% Run the SVD on each matrix
switch svdFunc
    case 'svd'
        if ~isnan(svdArg)
            [~,S,V] = svd(C, svdArg);  % Economy sized svd
        else
            [~,S,V] = svd(C);          % Full svd
        end

    case 'svds'
        [~,S,V] = svds(C, svdArg);     % Svds with nEigs
end

% Pull the eigenvalues off of the diagonal
eigvals = diag(S);

% Set the eigenvectors so that the majority of elements are positive
eigvecs = posColSign(V);
        
end
    
    
%%%%% Helper Functions %%%%%
function[svdFunc, svdArg] = parseInputs( varargin)
inArgs = varargin;

% Set the defaults
svdFunc = 'svd';
svdArg = NaN;

if ~isempty(inArgs)
    isSvdsArg = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
       
        if isSvdsArg
            if isscalar(arg)
                svdFunc = 'svds';
                svdArg = arg;
            else
                error('The svds flag must be followed by nEigs');
            end
        elseif strcmpi(arg, 'econ')
            svdFunc = 'svd';
            svdArg = 'econ';            
        elseif strcmpi(arg,'svd')
            % Do nothing
        elseif strcmpi(arg, 'svds') 
            if length(inArgs)>=k+1
                isSvdsArg = true;
            else
                error('The svds flag must be followed by nEigs');
            end            
        else
            error('Unrecognized Input');
       end
   end
end
end
          
function[] = errCheck(C)
% C is a matrix
if ~ismatrix(C)
    error('C must be a matrix');
% No NaN
elseif hasNaN(C)
    error('C may not contain NaN');
end
end

