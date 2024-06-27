function fout = tilefigure(varargin)
%TILEFIGURE Creates tiled figures
%   fh = tilefigure(tileSize,key,value,...) starts creating tiled figures
%   from top-left corner. The size of figure is determined by tileSize =
%   [number of columns, number of rows] (default value is tileSize = [4,
%   3]). The tile size can be followed by the key-value pairs supported by
%   figure command. The tiling is restarted when different tileSize is
%   given

persistent tileSize numberOfDisplayedFigures



% Get number of optional arguments
numberOfVarArgs = length(varargin);

if numberOfVarArgs>0
    % If some input arguments are given thre are posibilities 1) first
    % argument is vector of tile size possibly folowed by key-value pairs
    % 2) only key-value pairs are given
    if isnumeric(varargin{1}) && mod(numberOfVarArgs-1,2)==0 && iscellstr(varargin(2:2:end))
        
        % Reset the number of ploted tiles if a nan value if given in tile dimension and
        % terminate
        if any(isnan(varargin{1}(:)))
            numberOfDisplayedFigures = 0;
            return
        % Reset tiles if different tile size is given and continue
        elseif ~isequal(tileSize(:),varargin{1}(:))
            tileSize = varargin{1};
            numberOfDisplayedFigures = 0;
        end
        
        % Remove the first argument from varargin
        varargin(1) = [];
    elseif mod(numberOfVarArgs,2)~=0 || ~iscellstr(varargin(1:2:end))
        error('Incorrect input arguments')
    else
        
    end
        
end

% If the function is called for the first time and tile size input argument
% is not given,set default tile size and set the number of displayed
% figures to zero
if isempty(tileSize)
    tileSize = [4 3];
    numberOfDisplayedFigures = 0;
end


% Compute maximum number of figures
maxFigures = tileSize(1)*tileSize(2);


if numberOfDisplayedFigures == maxFigures
    error('Cannot display more than %i figures.',maxFigures)
end

% Compute normalized width and height of figure
normWH = 1./tileSize;

% Increase the number of diplayed figures
numberOfDisplayedFigures = numberOfDisplayedFigures + 1;

% Compute indices of tile for current figure
[col,row] = ind2sub(tileSize,numberOfDisplayedFigures);

% Note the screen is indexed bottom-up, therefore tileSize(2)-row is used
% to get top-down ordering

% Display figure
f = figure('Units','normalized',...
    'OuterPosition',[normWH(1)*(col-1) normWH(2)*(tileSize(2)-row) normWH],...
    varargin{:});

% Set Units of figure to default 'pixels'
set(f,'Units','pixels');


% Assign output argument if required
if nargout==1
    fout = f;
end


end