function v = multiindex(n,p)
% Recursive determination of all possible multi-indices
% 0 <= j_1 , j_2 , ... , j_n satisfying j_1 + j_2 + ... + j_n <= p
%
% Input parameters:
%  n   number of components in multi-indices
%  p   maximum degree
%
% Output parameter:
%  v   matrix containing multi-indices
%       first column:  linear numbering (row index)
%       second column: degree of multi-index in this row
%       other columns: multi-indices j_1 , ... , j_n
%

I = multiindex_2(n,p);

[nindex,dummy] = size(I);

v = zeros(nindex,n+2);

v(:,3:end) = I;

v(:,1) = 1:nindex;

v(:,2) = sum(I,2);

% end of main part of function multiindex.m


function I_mp=multiindex_2(m,p,combine,varargin)
% MULTIINDEX Generate a table of multiindices.
%   I_MP=MULTIINDEX(M,P,COMBINE,OPTIONS) generate a table of multiindices using
%   the standard block scheme i.e. generating all multi-indices up to
%   degree P in all M (random) variables. (Limitation to certain
%   limiters/norms will be added later). If combine is not specified or
%   evaluates to true then the homogeneous multiindices will be combined
%   into one large (sparse) array I_MP. Otherwise I_MP is a cell array
%   where I_MP{q+1} represents the multiindices with degree q.
%
% Options:
%   use_sparse: true, {false}
%     Return the result as a sparse array.
%   lex_ordering: true, {false}
%     Returns the result in lexicographical ordering (like e.g. A. Keese) instead
%     of ordering by degree first (this option is obviously ignored if COMBINE is
%     false)
%   full: true, {false}
%     Return the full tensor product multiindex set.
%
% Example (<a href="matlab:run_example multiindex">run</a>)
%   % To generate the polynomial chaos for 2 random variables up to
%   % polynomial order 4
%   I=multiindex(2,4);
%   disp(I);
%   % Get output as sparse array
%   I=multiindex(2,4,[],'use_sparse',true);
%   disp(I); % convert from sparse
%
%   % To generate the polynomial chaos for 5 random variables up to
%   % polynomial order 3, using only the homogeneous chaos of order 3
%   I=multiindex(5,3,false);
%   I3=full(I{3+1});
%   disp(I3)
%
% See also MULTIINDEX_ORDER, MULTIINDEX_COMBINE, MULTIINDEX_FACTORIAL

%   Elmar Zander
%   Copyright 2006, Institute of Scientific Computing, TU Braunschweig.
%   $Id$
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

error( nargchk( 2, inf, nargin ) );

if nargin<3 || isempty(combine)
    combine=true;
end

options=varargin2options( varargin );
[use_sparse,options]=get_option( options, 'use_sparse', false );
[lex_ordering,options]=get_option( options, 'lex_ordering', false );
[full,options]=get_option( options, 'full', false );
check_unsupported_options( options, mfilename );

if full
    I_mp=multiindex_full(m,p,use_sparse);
else
    I_mp=multiindex_complete(m,p,use_sparse);
end

if combine
    I_mp=cell2mat( I_mp(:) );
    if lex_ordering
        I_mp=sortrows(I_mp,m:-1:1);
    end
end


function I_kp=multiindex_full(m,p,use_sparse)
I_kp=cell(1,m*p+1);
for q=0:m*p
    if use_sparse
        I_kp{q+1}=sparse(q==0,0);
    else
        I_kp{q+1}=zeros(q==0,0);
    end
end

% Now iterate over the number of random variables.
for k=1:m
    % Backup the old multiindex set for later use.
    I_k1p=I_kp;

    for q=0:k*p
        if use_sparse
            I_kp{q+1}=sparse(0,k);
        else
            I_kp{q+1}=zeros(0,k);
        end
    end
    for q=0:(k-1)*p
        for r=0:p
            I_kp{q+r+1}=[I_kp{q+r+1}; [I_k1p{q+1}, r*ones(size(I_k1p{q+1},1),1)]];
        end
    end
end


function I_kp=multiindex_complete(m,p,use_sparse)
% The (old) idea of the algorithm is the following:
% We do a recursion on the number of random variables, not on the order (in
% my opinion its easier and faster that way). For just one random variable
% the multiindices are then trivial (0..p). For m+1 random variables we
% take the result from m random variables, which are stored by order, and
% for each order of the new set take all sets from m with lower order and
% add the remaining difference as last column. At the end we combine all
% bins which contain multiindices of homogeneous order into one large set.
% An advantage of this approach, besides its simplicity, is that large
% parts can be vectorized, and it runs pretty fast even for big sets.
%
% This has now been changed (more or less trivially) to a non-recursive
% algorithm. We just start with one random variable, create the
% multiindex set, and then iteratively build the multiindex sets for one
% more random variable each time. (BTW: the reason for this change was that
% for large number of random variable, i.e. about 100 or more, the matlab stack
% was exhausted by the recursive algorithm).

% Note: algorithm has been changed to start from m=0
% Start with one random variable. Create a cell array like this
%   {[0],[1],[2],...,[p]}
% Since we have only one random variable (m=1), in each cell i there is
% just one monomial of homogeneous order q=i-1
I_kp=cell(1,p+1);
for q=0:p
    if use_sparse
        %I_kp{q+1}=sparse(1,1,q);
        I_kp{q+1}=sparse(q==0,0);
    else
        %I_kp{q+1}=q;
        I_kp{q+1}=zeros(q==0,0);
    end
end

% Now iterate over the number of random variables.
for k=1:m
    % Backup the old multiindex set for later use.
    I_k1p=I_kp;

    % Get number of nonzero elements and number for multiindex set I_mp
    % nonzero and count are arrays that contain the respective values
    % indexed by order of the homogeneous indices (or polynomials). Then
    % use this information to allocate (sparse) arrays of the right size.
    [count,nonzero]=multiindex_stats(k,p);
    I_kp=cell(1,p+1);
    for q=0:p
        if use_sparse
            I_kp{q+1}=spalloc(count(q+1),k,nonzero(q+1));
        else
            I_kp{q+1}=zeros(count(q+1),k);
        end
    end

    for q=0:p
        % Copy indices from m-1 random vars in the new multiindex field to
        % the right position
        I_kp{q+1}( :, 1:end-1 ) = catmat(I_k1p{(q+1):-1:1});
        % Now fill the right most column such that I_kp{q+1} still has
        % homogeneous order q
        I_kp{q+1}(:,end)=q-sum(I_kp{q+1},2);
    end
end

function [count,nonzero]=multiindex_stats(m,p)
% MULTIINDEX_STATS Compute number of multiindices and of non-zero exponents.
count=ones(1,p+1);
nonzero=[0 ones(1,p)];

for k=2:m
    for q=p:-1:0
        count(q+1)=sum(count(1:(q+1)));
        nonzero(q+1)=sum(nonzero(1:(q+1))) + sum(count(1:q));
    end
end


function A=catmat( varargin )
% CATMAT Concatenate multiindices from a cell array into one.

% Due to some stupidity one side of the Mathworks CAT does not work if only
% one array if passed (for CAT could just return the array as it is...). So
% we have to do that here ourselves. Given, if you explicitly call CAT with
% just one argument that would not make sense, but if your code shall be
% oblivious as to the number of arrays you're processing it does indeed.
if length(varargin)==1
    A=varargin{1};
else
    A=cat(1,varargin{:});
end

% workaround for an octave bug
if issparse(varargin{1}) && ~issparse(A)
    A=sparse(A);
end


function ok=check_unsupported_options( options, mfilename )
% CHECK_UNSUPPORTED_OPTIONS Check whether unsupported options were passed.
%   OK=CHECK_UNSUPPORTED_OPTIONS( options, mfilename ) checks whether there
%   are entries left over in the OPTIONS struct, indicating that the user
%   probably called the function with a wrong or unsupported option. This
%   function provided a nice way to inform the user of such mistakes (most
%   probably spelling mistakes). CHECK_UNSUPPORTED_OPTIONS returns true, if
%   there were no fiels left over in OPTIONS, and false otherwise. In the
%   latter case a warning about unsupported options for function MFILENAME
%   is issued.
%
%   Note 1: for this method to work you have to call GET_OPTION with two
%   output arguments, so that used options are eliminated from the options
%   structure.
%
%   Note 2: pass mfilename literally for the second argument (i.e. pass the
%   return value of the buildin function 'mfilename' which tells you the
%   name of the current script, and is thus exactly what you want.)
%
% Example (<a href="matlab:run_example check_unsupported_options">run</a>)
%   % declare your own function to take varargs as options
%   %   function my_function( arg1, arg2, varargin )
%   % for demonstration we use args and set it arbitratily
%     args={'opt1', 'val1', 'usopt1', 'foo', 'usopt2', 'bar' };
%     options=varargin2options( args );
%     mfile='my_function';
%     [opt1,options]=get_option( options, 'opt1', 'default1' );
%     [opt2,options]=get_option( options, 'opt2', 'default2' );
%     [opt3,options]=get_option( options, 'opt3', 'default3' );
%     check_unsupported_options( options, mfile );
%
% See also VARARGIN2OPTIONS, GET_OPTION

%   Elmar Zander
%   Copyright 2007, Institute of Scientific Computing, TU Braunschweig.
%   $Id$
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

fields=fieldnames( options );
if isempty(fields) || (length(fields)==1 && strcmp( fields{1}, 'supported_fields__') )
    if nargout>0; ok=true; end
    return;
end

ok=false;
if isfield( options, 'supported_fields__' )
    supported_fields=options.supported_fields__;
    options=rmfield( options, 'supported_fields__' );
else
    supported_fields=[];
end

fields=fieldnames(options);
for i=1:length(fields)
    mode='warning';
    message=sprintf( 'unsupported option detected: %s', fields{i} );
    check_boolean( ok, message, mfilename, 'depth', 2, 'mode', mode );
end

if ~isempty( supported_fields )
    fields=sprintf( ', %s', supported_fields{:} );
    fprintf( 'Valid options for "%s" are: %s \n', mfilename, fields(3:end) );
end

if nargout==0
    clear ok;
end


function [val,options]=get_option( options, field, default )
% GET_OPTION Get a user option or return the default.
%   VAL=GET_OPTION( OPTIONS, FIELD, DEFAULT ) return the value of
%   OPTIONS.FIELD where FIELD is a string containing the structure field
%   containing the option, or DEFAULT is the field is not present. Useful
%   inside functions that can have a bunch of optional arguments.
%   If OPTIONS is also specified as output argument the field (if present)
%   is removed from the struct. This feature can be used to make sure only
%   valid options were passed to the function (struct should be empty after
%   all options have been queried.)
%
% Example (<a href="matlab:run_example get_option">run</a>)
%   function retval=my_function( arg1, arg2, arg3, varargin );
%     options = varargin2options( varargin );
%     [option1,options] = get_option( options1, 'option1', 1234 );
%     check_unsupported_options( options, mfilename );
%
% See also VARARGIN2OPTIONS, CHECK_UNSUPPORTED_OPTIONS

%   Elmar Zander
%   Copyright 2006, Institute of Scientific Computing, TU Braunschweig.
%   $Id$
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.



if ~isstruct(options)
    error( 'util:get_option:wrong_param', 'First argument to get_option must be a struct (maybe you interchanged options and field?)' );
elseif ~ischar(field)
    error( 'util:get_option:wrong_param', 'Second argument to get_option must be a string (maybe you interchanged options and field?)' );
end;

% set field or use default
if isfield( options, field )
    val=options.(field);
else
    val=default;
end

% if second output argument present, store field in supported field list
% and remove the current field if it was present
if nargout>1
    if isfield( options, field )
        options=rmfield(options,field);
    end
    if ~isfield( options, 'supported_fields__' )
        options.supported_fields__={};
    end
    options.supported_fields__{end+1}=field;
end


function log_flush(varargin)
% LOG_FLUSH Short description of log_flush.
%   LOG_FLUSH Long description of log_flush.
%
% Example (<a href="matlab:run_example log_flush">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2011, Inst. of Scientific Computing, TU Braunschweig
%   $Id$ 
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

status=get(0,'Diary');
if strcmp( status, 'on' )
    diary( 'off' );
    diary( 'on' );
end


function stop_check

global last_check
if isempty(last_check)
    last_check=tic;
end

if toc(last_check)>=1
    log_flush();
    filename=fullfile( getenv('HOME'), 'sglib_stop' );
    if exist( filename, 'file' )
        delete(filename);
        disp('Stopping matlab, entering debugger, press F5 to continue' );
        keyboard;
    end
    last_check=tic;
end


function options=varargin2options( args )
% VARARGIN2OPTIONS Convert variable argument list to options structure.
%   OPTIONS=VARARGIN2OPTIONS( ARGS ) returns the variable arguments as
%   an options structure. This allows the user to pass the arguments in
%   different forms; see the following examples.
%   OPTIONS=VARARGIN2OPTIONS() returns an empty options structure.
%   OPTIONS=VARARGIN2OPTIONS( {STRARG1, VAL1, STRARG2, VAL2, ...}) returns an
%   a structure with the pairs STRARGN and VALN converted to fields in the
%   returned options structure (i.e. options.(STRARGN)=VALN).
%   OPTIONS=VARARGIN2OPTIONS( {OPTS} ) returns the options structure as it
%   was passed to this function.
%
% Example (<a href="matlab:run_example varargin2options show">run</a>)
%   % declare your own function to take varargs as options
%   function my_function( arg1, arg2, varargin )
%     options=varargin2options( varargin );
%   % now suppose my_function has a debug option which can take on boolean
%   % values true or false, then you can call as either ...
%   my_function( arg1, arg2, 'debug', true );
%   % ... or ...
%   options.debug = true;
%   my_function( arg1, arg2, options );
%
% See also GET_OPTION

%   Elmar Zander
%   Copyright 2006, Institute of Scientific Computing, TU Braunschweig.
%   $Id$
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

stop_check;

% we need a cell array (don't pass varargin{:} any more, always pass
% varargin itself)
if ~iscell(args)
    error( 'util:varargin2options:no_cell', 'Argument must be a cell array' );
end

% if no args present return empty structure
if isempty(args)
    options=struct();
    return;
end

% one arg as structure means we can return the struct as is
if isstruct( args{1} )
    if length(args)>2
        error( 'util:varargin2options:struct', 'Wrong option specification: when struct is given only one argument is allowed' );
    end
    options=args{1};
    return
end

% options wrapped in a cell
if length(args)==1 && iscell(args{1})
    args=args{1};
end

names=args(1:2:end);
values=args(2:2:end);
if ~iscellstr(names)
    error( 'util:varargin2options:invalid_names', 'Wrong option specification: not all option names are strings: %s', evalc( 'disp(names);' ) );
end
if length(names)~=length(values)
    error( 'util:varargin2options:missing_value', 'Wrong option specification: not all option names have a corresponding value' );
end
[unames,ind]=unique(names);
unames;%#ok (we only need ind)
options=cell2struct( values(ind), names(ind), 2 );



