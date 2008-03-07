function [NORMALIZED_SPOTS SITE_NAME SITE_MATCH]=ParseDatabaseResults(FILENAME,HIV_ANNOT_CELL,DBTYPE)
%   ParseELMResults
%       Parses a results file from http://elm.eu.org/links.html and
%       determines the corresponding NT position of each motiff result.
%
%
%       [NORMALIZED_SPOTS SITE_NAME SITE_MATCH] = ...
%               ParseELMResults(FILENAME,HIV_ANNOT_CELL)
%
%       FILENAME        The filename of an ELM formated file
%       HIV_ANNOT_CELL  A cell array with the following format:
%                   {[gene_start gene_stop]} {GENE_NAME}
%
%       NORMALIZED_SPOTS    The location of each ELM entry (start stop)
%                           relative to the genome start.
%       SITE_NAME           The ELM entry name.
%       SITE_MATCH          The matching AA sequence found in by ELM.
%
%





switch lower(DBTYPE)
    case 'elm'
        elm_fid=fopen(FILENAME,'rt');
        if elm_fid==-1
            error('ParseDatabaseResults:BAD_FILE','The file %s could not be opened.',FILENAME);
        end

        try
            raw_data=textscan(elm_fid,'%s%n%n%s%s%*s','delimiter','\t');
        catch
            fclose(elm_fid);
            error('ParseDatabaseResults:BAD_FORMAT','A bad file format was provided.')
        end
        fclose(elm_fid);

        %adjust for AA-NT converstion
        starts=raw_data{2}*3;
        stops=raw_data{3}*3;

        gene_starts=cell2mat(HIV_ANNOT_CELL(:,1));
        [TF LOC]=ismember(lower(raw_data{1}),HIV_ANNOT_CELL(:,2));

        NORMALIZED_SPOTS=[starts(TF) stops(TF)]+repmat(gene_starts(LOC(TF),1),[1 2]);
        SITE_NAME=raw_data{4};
        SITE_MATCH=raw_data{5};

    case {'match','transfac'}
        match_fid=fopen(FILENAME,'rt');
        if match_fid==-1
            error('ParseDatabaseResults:BAD_FILE','The file %s could not be opened.',FILENAME);
        end

        raw_data=textscan(match_fid,'%s%n%s%*n%*n%s%q','delimiter','\t');

        NORMALIZED_SPOTS=[raw_data{2} raw_data{2}+cellfun('length',raw_data{4})];
        SITE_NAME=raw_data{5};
        SITE_MATCH=raw_data{4};

    otherwise
        error('Unknown Database Type provided')
end

end