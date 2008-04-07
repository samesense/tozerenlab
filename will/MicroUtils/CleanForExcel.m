function ExcelOutput=CleanForExcel(ExcelInput,varargin)
%   CleanForExcel
%       A utility function which formats cells in preperation for xlswrite.
%       It fills empty cells with a space characted (' ') and concatenates
%       nested cells into a single cell.
%
%
%   ExcelOutput=CleanForExcel(ExcelInput)
%
%
%   Optional Properties
%
%   EMPTY_CHAR          Defines the character(s) which are used to fill
%                       empty cells. DEFAULT = ' '
%
%   SPACE_CHAR          Defines the chararcter(s) which are used for
%                       concatenating nested cells. DEFAULT = ', '
%
%


EMPTY_CHAR=' ';
SPACE_CHAR=', ';

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i+1})
            case 'empty_char'
                if ischar(varargin{i+1})
                    EMPTY_CHAR=varargin{i+1};
                else
                    error('CleanForExcel:BAD_EMPTYCHAR','Arguement for EMPTY_CHAR must be a char-array')
                end

            case 'space_char'
                if ischar(varargin{i+1})||isempty(varargin{i+1})
                    SPACE_CHAR=varargin{i+1};
                else
                    error('CleanForExcel:BAD_SPACECHAR','Arguement for SPACE_CHAR must be a char-array or []')
                end
                
            otherwise
                error('CleanForExcel:UNKNOWN_ARG','An unknown arguement was provided: %s',varargin{i+1})
        end
    end
end




ExcelOutput=ExcelInput;

for i=1:numel(ExcelInput)
    if iscell(ExcelInput{i})&&~isempty(ExcelInput{i})
        temp=ExcelInput{i};
        ExcelOutput{i}=temp{1};
        for k=2:length(temp)
            ExcelOutput{i}=[ExcelOutput{i} SPACE_CHAR temp{k}];
        end
    elseif isempty(ExcelInput{i})
        ExcelOutput{i}=EMPTY_CHAR;
    end

end