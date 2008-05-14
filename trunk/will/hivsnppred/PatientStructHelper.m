function varargout=PatientStructHelper(INPUT_STRUCT,varargin)
%   PatientStructHelper
%       A helper function for adding and retrieving fields from cell arrays
%       of structures.
%
%
%
%   OUTPUT_STRUCT=PatientStructHelper(INPUT_STRUCT,...
%               {'field_name1',INDEX_ARRAY1,DATA1},...
%               {'field_name2',INDEX_ARRAY2,DATA2},
%               {'field_nameN',INDEX_ARRAYN,DATAN})
%
%       In this syntax the a field is added to the INPUT_STRUCT from the
%       data in the order according to the INDEX_ARRAY.  Any data that
%       already exists in the INPUT_STRUCT is carried over to the
%       OUTPUT_STRUCT.  If there is already a field with the new field_name
%       the old data is REPLACED with the new_data.
%
%
%
%   [DATA1 IND_ARRAY1 DATA2 IND_ARRAY2 ... DATAn IND_ARRAYn ] =
%       PatientStructHelper(INPUT_STRUCT, {'field_name1', 'FLAG1'}, ...
%                           {'field_name2', 'FLAG2'},...
%                           {'field_nameN', 'FLAGn'})
%
%
%       In this syntax the function will retrieve the data from the input
%       struct and deal it into the output maxtrices.  Where and IND_ARRAY
%       is the same size as DATA and indexes which row in INPUT_STRUCT
%       produced the corresponding row in DATA.
%
%       FLAG properties control the way in which the INPUT_STRUCT is parsed
%       into the DATA matrix.
%
%           LeaveCells      Leaves the data inside cells or encapusulates
%                           multiple rows into a single cell.  This ensures
%                           that the DATA matrix is the same size as
%                           INPUT_STRUCT.  Any INPUT_STRUCTs missing data
%                           are returned as empty-cells in DATA.  IND_ARRAY
%                           will be 1:length(INPUT_STRUCT) by definition.
%
%           ExplodeCells    Any INPUT_STRUCT which has a multi-row field
%                           will be 'EXPLODED' so that each row of DATA
%                           corresponds to only 1 row of 1 INPUT_STRUCT.
%
%           LeaveNumeric    The contents of the field in INPUT_STRUCT will
%                           be wrapped into a cell in DATA.  This ensures
%                           that the DATA matrix is the same size as
%                           INPUT_STRUCT.  Any INPUT_STRUCTs missing data
%                           are returned as empty-cells in DATA.  IND_ARRAY
%                           will be 1:length(INPUT_STRUCT) by definition.
%
%           ExplodeNumeric  Any INPUT_STRUCT which has a multi-row field
%                           will be 'EXPLODED' so that each row of DATA
%                           corresponds to only 1 row of 1 INPUT_STRUCT.
%

if nargout==1
    %%%%%Want to add fields to the structure
    FINAL_STRUCT=INPUT_STRUCT;
    for i=1:length(varargin)
        if iscell(varargin{i})&&length(varargin{i})==3
            temp_arg=varargin{i};

            data=temp_arg{3};

            if ischar(temp_arg{1})
                field_name=temp_arg{1};
            else
                warning('PatientStructHelper:BAD_FIELD_NAME','Arguement group-%d was not added because FIELD_NAME was not a char-array.',i);
                continue
            end

            if isvector(temp_arg{2})
                ind_array=temp_arg{2};
                if max(ind_array)>length(INPUT_STRUCT)
                    warning('PatientStructHelper:OUT_OF_BOUNDS','Arguement group-%d was not added because IND_ARRAY contained a reference outside of INPUT_STRUCT.',i);
                    continue
                elseif any(ind_array)<0
                    warning('PatientStructHelper:NEG_INDEXES','Arguement group-%d was not added because IND_ARRAY contained negative values at locations [%d].',i,find(temp_arg<0));
                    continue
                end

            else
                warning('PatientStructHelper:BAD_IND_ARRAY','Arguement group-%d was not added because IND_ARRAY was not a vector.',i);
                continue
            end


            %%%%%%IFF everything conforms then add the data
            FINAL_STRUCT=AddField(FINAL_STRUCT,ind_array,field_name,data);

        else
            warning('PatientStructHelper:BAD_FIELD','An instruction must be a 1x3 cell array containing {field_name,INDEX_ARRAY,DATA}.  Arguement group-%d was not added.',i);
            continue
        end


    end

    %%%%%%%%%%%After all are parsed return the resulting structure

    varargout{1}=FINAL_STRUCT;

else

    if nargout/2==length(varargin)
        varargout=cell(nargout,1);

        for i=1:length(varargin)
            temp_arg=varargin{i};

            if ~ischar(temp_arg{1})
                error('PatientStructHelper:BAD_FIELD_NAME','FIELD_NAME must be a char-array');
            end

            switch lower(temp_arg{2})
                case {'leavecells','leavecell'}
                    RETURN_CELL=true;
                    [varargout{2*i} varargout{2*(i-1)+1}]=ReturnCellFieldGrouped(INPUT_STRUCT,temp_arg{1});
                case {'explodecells','explodecell'}
                    RETURN_CELL=true;
                    [varargout{2*i} varargout{2*(i-1)+1}]=ReturnCellFieldNotGrouped(INPUT_STRUCT,temp_arg{1});
                case 'leavenumeric'
                    RETURN_CELL=false;
                    [varargout{2*i} varargout{2*(i-1)+1}]=ReturnNumericFieldGrouped(INPUT_STRUCT,temp_arg{1});
                case 'explodenumeric'
                    RETURN_CELL=false;
                    [varargout{2*i} varargout{2*(i-1)+1}]=ReturnNumericFieldNotGrouped(INPUT_STRUCT,temp_arg{1});

                otherwise
                    error('PatientStructHelper:BAD_ARG','An unknown arguement was provided: %s.',temp_arg{2});
            end
        end
        
    else
        error('PatientStructHelper:INCORRECT_SYNTAX','The number of output arguement PAIRS [%d] does not equal the number of FIELD requests [%d].',nargout/2,length(varargin));
    end

end

    function new_struct=AddField(old_struct,index_array,field_name,data)
        new_struct=cell(size(old_struct));
        for k=1:length(new_struct)
            temp_struct=old_struct{k};
            if isfield(temp_struct,field_name)
                rmfield(temp_struct,field_name);
            end

            spots=find(index_array==k);
            if ~isempty(spots)
                temp_struct=setfield(temp_struct,field_name,data(spots,:));
            end
            new_struct{k}=temp_struct;
        end
    end

    function [index_array data]=ReturnCellFieldGrouped(input_struct,field_name)

        data=cell(size(input_struct));
        index_array=(1:length(data))';
        for k=1:length(input_struct)
            if isfield(input_struct{k},field_name)
                data{k}=getfield(input_struct{k},field_name);
            end
        end

    end

    function [index_array data]=ReturnCellFieldNotGrouped(input_struct,field_name)
        T_data=cellfun(@(x)GetData(x,field_name),input_struct(:),'uniformoutput',false);

        data=cell(10*length(T_data),1);
        index_array=zeros(size(data));
        counter=1;

        for k=1:length(T_data)
            temp=T_data{k};
            if ~isempty(temp)
                num_elements=length(temp);
                data(counter:counter+num_elements-1)=temp(:);
                index_array((counter:counter+num_elements-1))=k;
                counter=counter+num_elements;
            end
        end
        data=data(1:counter-1,:);
        index_array=index_array(1:counter-1);
    end

    function [index_array data]=ReturnNumericFieldNotGrouped(input_struct,field_name)
        T_data=cellfun(@(x)GetData(x,field_name),input_struct(:),'uniformoutput',false);

        data=zeros(10*length(T_data),1);
        width=1;
        index_array=zeros(size(data));
        counter=1;

        for k=1:length(T_data)
            temp=T_data{k};
            if ~isempty(temp)
                num_elements=size(temp,1);
                if size(temp,2)>width
                    data=[data NaN(10*length(T_data),size(temp,2)-width)];
                    width=size(data,2);
                end

                data(counter:counter+num_elements-1,:)=temp;
                index_array((counter:counter+num_elements-1))=k;
                counter=counter+num_elements;
            end
        end
        data=data(1:counter-1,:);
        index_array=index_array(1:counter-1);
    end

    function data1=GetData(input,field_name)
        if isstruct(input)&&isfield(input,field_name)
            data1=getfield(input,field_name);
        else
            data1={};
        end
    end

    function [index_array data]=ReturnNumericFieldGrouped(input_struct,field_name)

        data=zeros(size(input_struct));
        index_array=(1:length(data))';
        for k=1:length(input_struct)
            data(k)=getfield(input_struct{k},field_name);
        end

    end

end
