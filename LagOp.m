classdef (InferiorClasses = {?internal.econ.LagIndexedArray}, Sealed) LagOp
%LAGOP Create a lag operator polynomial (LagOp) object
%
% Syntax:
%
%   A = LagOp(coefficients)
%   A = LagOp(coefficients,param1,val1,param2,val2,...)
%
% Description:
%
%   Create a lag operator polynomial, A(L), by specifying the coefficients
%   and, optionally, the corresponding lags.
%
% Input Argument:
%
%   coefficients - The coefficients of the lag operator polynomial.
%     Generally, coefficients is a cell array of square matrices. For
%     convenience, coefficients may also be specified in other ways:
%
%     o As a vector, representing a univariate time series polynomial with
%       multiple lags.
%
%     o As a matrix, representing a multivariate time series polynomial
%       with a single lag.
%
%     o As an existing LagOp object, to be updated according to the
%       optional inputs.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Lags'       Vector of integer lags associated with the polynomial 
%                coefficients. If specified, the number of lags must be the
%                same as the number of coefficients. If unspecified,
%                coefficients are associated with lags 0, 1, ...,
%                numCoefficients-1
%
%   'Tolerance'  Nonnegative scalar tolerance used to determine which lags 
%                are included in the object. The default is 1e-12. Specifying
%                a tolerance greater than the default excludes lags with
%                near-zero coefficients. A lag is excluded if the magnitudes 
%                of all elements of the coefficient matrix are less than or 
%                equal to the specified tolerance.
%
% Output Argument:
%
%   A - Lag operator polynomial (LagOp) object.
 
% Copyright 2010 The MathWorks, Inc.
 
properties (Access = public)
  Coefficients % Lag indexed cell array of nonzero polynomial coefficients
end
 
properties (Dependent)
  Lags      % Polynomial lags associated with nonzero coefficient
  Degree    % Polynomial degree (the highest lag associated with a nonzero coefficient)
  Dimension % Polynomial dimension (the number of time series to which it may be applied)
end
 
properties (Constant,Hidden)
  ZeroTolerance = 1e-12;     % Element-wise tolerance for excluding zero coefficient matrices
  SingularTolerance = 1e+10; % Condition number tolerance for warning about near-singular A0 coefficient
end
 
methods % Dependent property GET methods
  function L = get.Lags(OBJ)
    L = getLags(OBJ.Coefficients);
  end
  function P = get.Degree(OBJ)
    L = OBJ.Lags;
    if isempty(L)
       P = 0;
    else
       P = max(L);
    end
  end
  function N = get.Dimension(OBJ)
    N = getDimension(OBJ.Coefficients);
  end
end
 
methods (Access = public)
 
%-------------------------------------------------------------------------
function OBJ = LagOp(Coefficients,varargin)
 
if nargin == 0
 
   % MATLAB classes should construct a scalar object in a default state in
   % the absence of input arguments. In this case, the default constructor
   % syntax is equivalent to an identity operator: If there are no inputs,
   % then a zero-lag, 1-D lag operator polynomial is created.
 
   Coefficients = 1;
end
 
% Check input arguments:
 
if nargin > 5
    
   error(message('econ:LagOp:LagOp:TooManyInputs'))
      
end
 
parser = inputParser;
parser.addRequired('Coefficients',@OBJ.checkCoefficients);
parser.addParamValue('Lags',0:numel(Coefficients)-1,@OBJ.checkLags);
parser.addParamValue('Tolerance',OBJ.ZeroTolerance,@OBJ.checkTolerance);
parser.parse(Coefficients,varargin{:});
 
Coefficients = parser.Results.Coefficients;
Lags         = parser.Results.Lags;
tolerance    = parser.Results.Tolerance;
 
if isa(Coefficients,'LagOp')
   if isempty(varargin(1:2:end)) || ~any(strcmpi('Lags',varargin(1:2:end)))
      Lags = union(Coefficients.Lags, 0);
   end
   C = cell(1, numel(Lags));
   for L = 1:numel(Lags)
       C{L} = Coefficients.Coefficients{Lags(L)};
   end
   Coefficients = C;
end
 
if any(strcmpi('Lags',varargin(1:2:end))) % Handle convenience forms
 
   if iscell(Coefficients) || isvector(Coefficients)
 
      % Coefficients specified as either:
      %
      % (1) a cell array of square matrices, or
      % (2) a vector (representing a 1-dimensional, multi-lag polynomial)
      %
      % In either case, the number of elements of Lags must be the same as
      % the number of elements of Coefficients.
 
      if numel(Lags) ~= numel(Coefficients)
          
         error(message('econ:LagOp:LagOp:InconsistentInputLengths'))
           
      end
   else
 
      % Coefficients specified as a square matrix (representing a 
      % multidimensional, single-lag polynomial). In this case, Lags must
      % be a scalar.
 
      if numel(Lags) > 1
          
         error(message('econ:LagOp:LagOp:MatrixCoefficientsNonScalarLags'))
           
      end
   end
 
else
 
   if isnumeric(Coefficients) && ~isvector(Coefficients)
      Lags = 0;
   end
 
end
 
if isnumeric(Coefficients) 
   if isvector(Coefficients)
      Coefficients = num2cell(Coefficients);
   else
      Coefficients = {Coefficients};
   end
end
 
% Retain only unique lags. Calling the UNIQUE function will sort the lags,
% and in the event of duplicate lags will also retain and assign the
% coefficient associated with the LAST instance. This is designed to mimic
% the behavior of core MATLAB cell arrays.
 
[Lags,indices] = unique(Lags);
Coefficients = Coefficients(indices);
Dimension = size(Coefficients{1},1);
 
% Identify and retain only sufficiently nonzero coefficient matrices:
 
indices = true(numel(Lags),1);
 
for L = 1:numel(Lags)
    indices(L) = any( (abs(Coefficients{L}(:)) > tolerance) | (isnan(Coefficients{L}(:))) ); 
end
 
Lags = Lags(indices);
Coefficients = Coefficients(indices);
OBJ.Coefficients = internal.econ.LagIndexedArray(Coefficients,Lags,Dimension);
 
end
 
end % Methods (Access = public)
 
methods (Hidden)
 
%-------------------------------------------------------------------------
function varargout = subsref(OBJ, args)
%SUBSREF Subscripted referencing of Lag Operator Polynomial (LagOp) objects
%
% Syntax:
%
% For lag operator polynomial object A, valid subscripted references are:
% 
%  A(Y)                   - Filter time series data Y (same as A.filter(Y))
%
%  A.Coefficients         - Get the contained LagIndexedArray object
%  A.Coefficients(...)    - Get a subset of the contained LagIndexedArray
%                           object
%  A.Coefficients{...}    - Get a subset of n-by-n coefficients matrices
%  A.Coefficients{L}(i,j) - Get elements (i,j) of a coefficient matrix at
%                           lag L
%
%  A.Lags      - Get polynomial lags associated with nonzero coefficients
%  A.Degree    - Get polynomial degree, defined as the highest lag 
%                associated with a specified nonzero coefficient
%  A.Dimension - Get dimension of each coefficient matrix (the number of
%                time series to which the LagOp object may be applied)
%
%  A.MethodName(...) - Invoke the specified method (same as
%                      MethodName(A,...))
 
nSubscripts = LagOp.checkSubscripts(args); % Error checking
 
try
 
   switch args(1).type
 
     case '()' % Filter time series data, A(Y(t))
 
       if nSubscripts > 1
 
          error(message('econ:LagOp:subsref:InvalidAccessType'));
          
       end
       
       if isempty(args(1).subs)
           
          error(message('econ:LagOp:subsref:MissingTimeSeries'))
            
       else
           
          [varargout{1:nargout}] = OBJ.filter(args(1).subs{:}); % A(Y(t)) == A.filter(Y(t))
          
       end
 
     case '.' % Property or method referencing
 
        switch args(1).subs
 
          case {'Degree','Dimension','Lags'} % Dependent properties
 
            % For dependent properties, the following are OK:
            %
            % (1) OBJ.Property
            % (2) OBJ.Property( )
            % (3) OBJ.Property(:)
 
            if (nSubscripts > 1) && ~(isempty(args(2).subs) || strcmp(args(2).subs,':'))
 
               error(message('econ:LagOp:subsref:InvalidAccessType'));
               
            else
                
               varargout = {OBJ.(args(1).subs)};
               
            end
 
          case {'Coefficients'} % Contained LagIndexedArray object
 
            if nSubscripts > 1
 
               if (nargout == 0) && strcmp(args(2).type,'{}') && ...
                  ( (numel(args(2).subs{:}) > 1) || (strcmp(args(2).subs,':') && (numel(OBJ.Lags) > 1)) )
 
                   % "{}" syntaxes with no LHS output arguments and more
                   % than one RHS index, such as
                   %
                   % (1) A.Coefficients{[0 2 3]}
                   % (2) A.Coefficients{[0:3]}
                   % (3) A.Coefficients{:} (but only when the object has
                   %                        more than one nonzero lag)
                   %
                   % are disallowed because VARARGOUT cannot effect a
                   % comma-delimited list in this situation.
 
                   error(message('econ:LagOp:subsref:InsufficientLHSArgumentList'))
                     
               else
 
                   if strcmp(args(2).type, '()') && (isempty(args(2).subs) || strcmp(args(2).subs, ':'))
 
                   %  Handle the following:
                   %
                   %  (1) A.Coefficients()
                   %  (2) A.Coefficients(:)
 
                      [varargout{1:nargout}] = OBJ.Coefficients;
 
                   else
 
                   %  Handle the following:
                   %
                   %  (1) A.Coefficients{L} (single RHS lag index)
                   %  (2) [...] = A.Coefficients{[...]}
                   %  (3) C = A.Coefficients{L}(i,j)
                   %  (4) A.Coefficients([...])
                   %  (5) C = A.Coefficients([...])
 
                      [varargout{1:nargout}] = subsref(OBJ.Coefficients,args(2:end));
 
                   end
 
               end
 
            else % nSubscripts = 1
 
                [varargout{1:nargout}] = OBJ.Coefficients;
 
            end
 
          otherwise % Invoke a method
 
            [varargout{1:nargout}] = builtin('subsref',OBJ,args);
 
        end
 
     otherwise
 
        error(message('econ:LagOp:subsref:InvalidAccessType'));
 
   end
 
catch exception
 
   % Trap common errors to provide a more meaningful message:
   
   if strcmpi(exception.identifier,'MATLAB:unassignedOutputs')
       
      error(message('econ:LagOp:subsref:IncorrectLHSArgumentList'))
        
   else
       
      exception.throwAsCaller();
      
   end
 
end
 
end % SUBSREF
 
%-------------------------------------------------------------------------
function OBJ = subsasgn(OBJ,args,value)
%SUBSASGN Subscripted assignment of Lag Operator Polynomial (LagOp) objects
%
% Syntax:
%
% For LagOp object A, LagIndexedArray object B, and lag L, valid
% assignments are:
%
%  A.Coefficients         =  {...} - Assign all coefficients of a LagOp
%  A.Coefficients(...)    =  {...} - Assign subset of coefficients of a
%                                    LagOp
%  A.Coefficients         = B(...) - Assign all coefficients of a LagOp 
%                                    from a LagIndexedArray
%  A.Coefficients(...)    = B(...) - Assign subset of coefficients of a 
%                                    LagOp from a LagIndexedArray
%  A.Coefficients{L}      =  [...] - Assign subset of n-by-n coefficient
%                                    matrices
%  A.Coefficients{L}(i,j) =  [...] - Assign elements of a coefficient
%  matrix
 
nSubscripts = LagOp.checkSubscripts(args);
 
try
 
   switch args(1).subs
 
      case {'Degree','Dimension','Lags'} % Dependent properties
 
         error(message('econ:LagOp:subsasgn:DependentProperty'))
          
 
      case {'Coefficients'}
 
         if isa(value,'internal.econ.LagIndexedArray') % Assignment from a contained LagIndexedArray object
 
            if nSubscripts == 1
 
               % Assignments of the form: A.Coefficients = B(...):
 
               RHS_Lags = sort(getLags(value)); % RHS "assigned from" lags
               OBJ = LagOp({value{RHS_Lags(:)}},'Lags',RHS_Lags,'Tolerance',0); %#ok
 
            else
 
               OBJ = LagOp(subsasgn(OBJ.Coefficients,args(2:end),value));
 
            end
 
         elseif iscell(value) || isnumeric(value)      % Value is either a cell array or a matrix
 
            if nSubscripts == 1
 
               if iscell(value) && ~isempty(value)
                  [nR,nC] = size(value{1}); % Check dimensions
               else
                  [nR,nC] = size(value);    % Check dimensions
               end
 
               if (nR ~= OBJ.Dimension) || (nC ~= OBJ.Dimension)
                   
                  error(message('econ:LagOp:subsasgn:InconsistentDimension'))
                    
               end
 
               OBJ = LagOp(value);
 
            else
 
               OBJ = LagOp(subsasgn(OBJ.Coefficients,args(2:end),value));
 
            end
 
         else
 
            error(message('econ:LagOp:subsasgn:InvalidAccessType'));
 
         end
 
       otherwise
 
         error(message('econ:LagOp:subsasgn:InvalidAccessType'));
 
   end
 
catch exception
 
   exception.throwAsCaller();
 
end
 
end  % SUBSASGN
 
%-------------------------------------------------------------------------
function disp(OBJ)
% DISP Display method
 
D = OBJ.Dimension;   % Get the dimension
N = numel(OBJ.Lags); % Get the number of lags
 
spaces = '                    ';
fprintf( '    %d-D Lag Operator Polynomial:\n', D)
fprintf( '    -----------------------------\n')
string = [spaces(1:8) 'Coefficients:'];
 
% Print coefficient information:
 
if (N > 0) && (N <= 12) && (D == 1)
    C = OBJ.Coefficients;
    format = repmat(' %g',1,N);
    fprintf([string  ' [' format(2:end) ']\n'], C{OBJ.Lags})
else
    fprintf([string  ' [Lag-Indexed Cell Array with %d Non-Zero Coefficients]\n'], N)
end
 
% Print lags information:
 
if N <= 12
 
   % If short enough, then print them all:
 
   format = repmat(' %d',1,N);
   fprintf([spaces(1:16)      'Lags: [' format(2:end) ']\n'],OBJ.Lags)
elseif sum(diff(OBJ.Lags)) == (N-1)
 
   % Pretty-print consecutive lags:
 
   fprintf([spaces(1:16) 'Lags: [%d ... %d]\n'], min(OBJ.Lags), max(OBJ.Lags))
else
 
   % Too many non-consecutive lags, so just print how many there are:
 
   fprintf([spaces(1:16) 'Lags: [Lags Associated with %d Non-Zero Coefficients]\n'],N)
end
 
fprintf([spaces(1:14)    'Degree: %d\n'],OBJ.Degree)
fprintf([spaces(1:11) 'Dimension: %d\n'],OBJ.Dimension)
 
end
 
end % Methods (Hidden)
 
methods (Static, Hidden) % Error Checking Utilities
 
%-------------------------------------------------------------------------
function [nSubscripts, errorStruct] = checkSubscripts(args)
% CHECKSUBSCRIPTS Common subscript/index error checking
 
errorStruct.message = 'Invalid subscript or access syntax.';
errorStruct.identifier = 'econ:LagOp:checkSubscripts:InvalidAccessType';
 
nSubscripts = length(args);
 
% Ensure multiple references are of acceptable forms:
 
if nSubscripts > 3
 
   error(message('econ:LagOp:checkSubscripts:InvalidAccessType'));
   
end
 
if nSubscripts >= 2
   if ~strcmp('.', args(1).type)
      error(message('econ:LagOp:checkSubscripts:InvalidAccessType'));
   else
      if ~any(strcmp(args(2).type,{'()' '{}'}))
         error(message('econ:LagOp:checkSubscripts:InvalidAccessType'));
      end
   end
end
 
end
 
%-------------------------------------------------------------------------
function OK = checkCoefficients(coefficients)
% Check coefficients
 
if isa(coefficients,'LagOp')
   OK = true;
   return;
end
 
if isempty(coefficients)
    
   error(message('econ:LagOp:checkCoefficients:EmptyCoefficients'))
end
 
if iscell(coefficients)
    
   [nR1, nC1] = size(coefficients{1});
   
   if (nR1 ~= nC1) || (ndims(coefficients{1}) > 2)
       
      error(message('econ:LagOp:checkCoefficients:NonSquareCoefficientMatrix'))
        
   end
   
   for i = 2:numel(coefficients)
       [nR, nC] = size(coefficients{i});
       if (nR ~= nC) || (nR ~= nR1) || (nC ~= nC1) || (ndims(coefficients{i}) > 2)
           
          error(message('econ:LagOp:checkCoefficients:NonSquareCoefficientCells'))
            
       end
   end
   
elseif isnumeric(coefficients)
    
   if ~isvector(coefficients)
      [nR, nC] = size(coefficients);
      if (nR ~= nC) || (ndims(coefficients) > 2)
          
         error(message('econ:LagOp:checkCoefficients:NonSquareCoefficientMatrix'))
           
      end
   end
   
else
    
   error(message('econ:LagOp:checkCoefficients:NonSquareCoefficientMatrix'))
     
end
 
OK = true;
 
end
 
%-------------------------------------------------------------------------
function OK = checkLags(lags)
% Check lags
 
if ~isnumeric(lags)
    
   error(message('econ:LagOp:checkLags:LagsNonNumeric'))
     
elseif ~isvector(lags)
    
   error(message('econ:LagOp:checkLags:LagsIsArray'))
     
elseif any(mod(lags,1) ~= 0) || any(lags < 0)
    
   error(message('econ:LagOp:checkLags:NegativeLags'))
     
end
 
OK = true;
 
end
 
%-------------------------------------------------------------------------
function OK = checkDegree(degree)
% Check degree
 
if ~isnumeric(degree) || ~isscalar(degree) || (degree < 0) || (mod(degree,1) ~= 0)
    
   error(message('econ:LagOp:checkDegree:InvalidTolerance'))
end
 
OK = true;
 
end
 
%-------------------------------------------------------------------------
function OK = checkWindow(window)
% Check the MRDIVIDE/MLDIVIDE termination window
 
if ~isnumeric(window) || ~isscalar(window) || (window <= 0) || (mod(window,1) ~= 0)
    
   error(message('econ:LagOp:checkWindow:InvalidWindow'))
     
end
 
OK = true;
 
end
 
%-------------------------------------------------------------------------
function OK = checkTolerance(tolerance)
% Check tolerance
 
if ~isnumeric(tolerance) || ~isscalar(tolerance) || (tolerance < 0)
    
   error(message('econ:LagOp:checkTolerance:InvalidTolerance'))
     
end
 
OK = true;
 
end
 
end % Methods (Static, Hidden)
 
end % Class definition
