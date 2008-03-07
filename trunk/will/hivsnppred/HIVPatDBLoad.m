function [PAT_DB]=HIVPatDBLoad(DIRECTORY,FLAG)
%   HIVPatDBLoad
%       Loads the Patient information in the format provided by the
%       http://hivdb.stanford.edu/ database.
%
%
%   PAT_DATA = HIVPatDBLoad(DIRECTORY)
%       This loads all of the data in the specified directory
%
%
%
%   PAT_DATA                A cell of structs containing the RX-history,
%                           CD4 history, and Viral-Load history.
%
%
%
%   This importer requires that the Databases have the following format:
%
%       "BASEFILENAME"_PR_fasta.txt     In FASTA format
%       "BASEFILENAME"_RT_fasta.txt     In FASTA format
%       "BASEFILENAME"_RNA.txt          Tab-delimited text file
%               PtID    RNAdate     Vload   Vmatch
%       "BASEFILENAME"_CD4.txt          Tab-delimited text file
%               PtID    CD4date     CD4count



if nargin==0
    file_list=dir;
else
    file_list=dir(DIRECTORY);
end





base_filename=[];

for j=1:length(file_list)
    if length(file_list(j).name)>14
        if strcmp(file_list(j).name(end-7:end),'_CD4.txt')
            base_filename=[base_filename; {file_list(j).name(1:end-8)}];
        elseif strcmp(file_list(j).name(end-7:end),'_RNA.txt')
            base_filename=[base_filename; {file_list(j).name(1:end-8)}];
        elseif strcmp(file_list(j).name(end-12:end),'_PR_fasta.txt')
            base_filename=[base_filename; {file_list(j).name(1:end-13)}];
        elseif strcmp(file_list(j).name(end-12:end),'_RT_fasta.txt')
            base_filename=[base_filename; {file_list(j).name(1:end-13)}];
        end
    end
end

base_filenames=unique(base_filename);

display('The following databases will be loaded:');
display(base_filenames);

display('Now loading')

PAT_DB=[];
for k=1:length(base_filenames)
    display(base_filenames(k))
    PAT_DB=[PAT_DB; GetOneDatabase(base_filenames(k))];
end







    function [THIS_DB]=GetOneDatabase(BASE)
        BASE_FILENAME=BASE{1};
        %%%%%%%%%%%%%%%%%%%%%%FILE CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        try
            [CD4file m1]=fopen([DIRECTORY '\' BASE_FILENAME '_CD4.txt']);
        catch
            fclose('all');
            warning('FindError:CD4File',['Cannot find ' BASE_FILENAME '_CD4.txt'])
            rethrow(lasterror);
        end

        try
            [RNAfile m2]=fopen([DIRECTORY '\' BASE_FILENAME '_RNA.txt']);
        catch
            fclose('all');
            warning('FindError:RNAFile',['Cannot find ' BASE_FILENAME '_RNA.txt'])
            rethrow(lasterror);
        end

        try
            [RXfile m3]=fopen([DIRECTORY '\' BASE_FILENAME '_RX.txt']);
        catch
            fclose('all');
            warning('FindError:RXFile',['Cannot find ' BASE_FILENAME '_RX.txt'])
            rethrow(lasterror);
        end


        %%%%%%%%%%%%%%%%%%%%%%%LOAD RX DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        try
            RX_header_line=fgetl(RXfile);
            RX_headers=textscan(RX_header_line,'%s');
            RX_headers=RX_headers{1};
            DRUGS=RX_headers(7:end);
            format=[];
            for i=1:length(RX_headers)
                if ismember(RX_headers(i),DRUGS)
                    format=[format '%d'];
                else
                    format=[format '%s'];
                end
            end
            RX_data=textscan(RXfile,format);
        catch
            fclose('all');
            rethrow(lasterror)
            error('FileError:RXFile','%s','Could not Process RX Data properly')
        end
        fclose(RXfile);

        RX_pat_ids=unique(RX_data{1});
        drug_spots=ismember(RX_headers,DRUGS);
        drug_data=[RX_data{drug_spots}];

        RX_pat_vectors=cell(length(RX_pat_ids),1);
        for i=1:length(RX_pat_ids)
            spots=strmatch(RX_pat_ids(i),RX_data{1});
            dates=str2double(RX_data{5}(spots));
            RX_vals=drug_data(spots,:);
            RX_pat_vectors{i}=[dates RX_vals];
        end

        %%%%%%%%%%%%%%%%%%%%%%LOAD CD4 DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        try
            CD4_header_line=fgetl(CD4file);

            CD4_headers=textscan(CD4_header_line,'%s');
            CD4_headers=CD4_headers{1};

            CD4_data=textscan(CD4file,'%s%n%n');
            %[PtID CD4Date CD4Count]
        catch
            fclose('all');
            rethrow(lasterror)
            error('FileError:CD4File','%s','Could not Process CD4 clinical data');
        end
        fclose(CD4file);

        CD4_pat_ids=unique(CD4_data{1});
        CD4_pat_vectors=cell(length(CD4_pat_ids),1);
        for i=1:length(CD4_pat_ids)
            spots=strmatch(CD4_pat_ids(i),CD4_data{1});
            dates=CD4_data{2}(spots);
            CD4_counts=CD4_data{3}(spots);
            CD4_pat_vectors{i}=[dates CD4_counts];
        end

        %%%%%%%%%%%%%%%%%%%%%LOAD RNA DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        try
            RNA_header_line=fgetl(RNAfile);
            RNA_headers=textscan(RNA_header_line,'%s');
            RNA_headers=RNA_headers{1};
            RNA_data=textscan(RNAfile,'%s%n%n%s');
            %[PtID RNADate VLoad VLoadmatch]
        catch
            fclose('all');
            rethrow(lasterror)
            error('FileError:RNAFile','%s','Could not process RNA file properly')
        end

        RNA_pat_ids=unique(RNA_data{1});
        RNA_pat_vectors=cell(length(RNA_pat_ids),1);
        for i=1:length(RNA_pat_ids)
            spots=strmatch(RNA_pat_ids(i),RNA_data{1});
            dates=RNA_data{2}(spots);
            RNA_counts=RNA_data{3}(spots);
            RNA_matches=zeros(length(spots),1)+(cell2mat(RNA_data{4}(spots))=='>')-(cell2mat(RNA_data{4}(spots))=='<');
            RNA_pat_vectors{i}=[dates RNA_counts RNA_matches];
        end




%        if ~isempty(FLAG)&&FLAG
            %%%%%%%%%%%%%%%%%%%%%%LOAD SEQUENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            empty_seq.Header='PtID 0 | Alias 302-150155 | Week 0';
            empty_seq.Sequence='';
            try
                PR_DATA=fastaread([DIRECTORY '\' BASE_FILENAME '_PR_fasta.txt']);
            catch
                PR_DATA=empty_seq;
            end
            try
                RT_DATA=fastaread([DIRECTORY '\' BASE_FILENAME '_RT_fasta.txt']);
            catch
                RT_DATA=empty_seq;
            end

            RT_info=arrayfun(@(x)(textscan(x.Header,'%*s %s %*s %*s %*s %*s %*s %n',1)),RT_DATA,'uniformoutput',false);
            PR_info=arrayfun(@(x)(textscan(x.Header,'%*s %s %*s %*s %*s %*s %*s %n',1)),PR_DATA,'uniformoutput',false);

            RT_info=[RT_info{:}];
            PR_info=[PR_info{:}];
            
            RT_dates=[RT_info{2:2:end} 0];
            PR_dates=[PR_info{2:2:end} 0];
            
            RT_pat_ids=[RT_info{1:2:end}];
            PR_pat_ids=[PR_info{1:2:end}];
            
           
            PR_DATA=[PR_DATA; empty_seq];
            RT_DATA=[RT_DATA; empty_seq];
            
%        end

        struct_list=genvarname([RX_headers;{'RNALoad';'CD4Count'}]);
        shared_pat_ids=intersect(intersect(RX_pat_ids,RNA_pat_ids),CD4_pat_ids);

        THIS_DB=cell(length(shared_pat_ids),1);
        for i=1:length(shared_pat_ids)
            RX_spot=strmatch(shared_pat_ids(i),RX_pat_ids);
            CD4_spot=strmatch(shared_pat_ids(i),CD4_pat_ids);
            RNA_spot=strmatch(shared_pat_ids(i),RNA_pat_ids);

%             if ~isempty(FLAG)&&FLAG
                PR_spot=strmatch(shared_pat_ids(i),PR_pat_ids);
                RT_spot=strmatch(shared_pat_ids(i),RT_pat_ids);
                
                %point things to a dummy sequence if its missing
                if isempty(PR_spot)
                    PR_spot=length(PR_DATA);
                end
                if isempty(RT_spot)
                    PR_spot=length(PR_DATA);
                end
                
                this_struct=struct('Study',BASE_FILENAME,'PtID',shared_pat_ids(i),'RX_vals',RX_pat_vectors{RX_spot},'RX_names',{DRUGS},...
                    'CD4_vec',CD4_pat_vectors{CD4_spot},'RNA_vec',RNA_pat_vectors{RNA_spot},'PR_seqs',{{PR_DATA(PR_spot).Sequence}},...
                    'PR_dates',PR_dates(PR_spot),'RT_seqs',{{RT_DATA(RT_spot).Sequence}},'RT_dates',RT_dates(RT_spot) );
%             else
% 
%                 this_struct=struct('Study',BASE_FILENAME,'PtID',shared_pat_ids(i),'RX_vals',RX_pat_vectors{RX_spot},'RX_names',{DRUGS},...
%                     'CD4_vec',CD4_pat_vectors{CD4_spot},'RNA_vec',RNA_pat_vectors{RNA_spot});
% 
%             end
            THIS_DB{i}=this_struct;
        end

    end
end























