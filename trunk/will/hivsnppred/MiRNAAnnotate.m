function OUT_STRUCT=MiRNAAnnotate(PAT_STRUCT,FILENAME)
%   MiRNAAnnotate
%       Loads the Patient Data from the pyHybrid refined-output file into
%       the PAT_STRUCT.
%
%
%   OUT_STRUCT=MiRNAAnnotate(PAT_STRUCT,FILENAME)
%
%   PAT_STRUCT          A patient structure that was created with
%                       PatStructHelper.
%
%   FILENAME            The path to the refined output file generate by
%                       pyHybrid.
%
%
%   OUT_STRUCT          The patient structure appeneded with the fields:
%                       SimpleHumanMiRNA
%                           A 1xN double array where the valueindicates the
%                           min p-value of that miRNA being found in the
%                           sequence.  NaN if the value was not found.
%                       HumanMiRNAnames
%                           A 1xN cell-array which has the corresponding
%                           miRNA labels.
%
%



fid = fopen(FILENAME,'rt');

BLOCK_SIZE = 10000;
seqNames = cell(BLOCK_SIZE,1);
data = cell(0,3);

while size(seqNames,1) == BLOCK_SIZE

    display(['Reading Block: ' num2str(size(data,1)/BLOCK_SIZE)])
    inputdata = textscan(fid,'%s%s%[^\n]',BLOCK_SIZE,'delimiter','\t');
    [seqNames miRNAnames pos_data]=deal(inputdata{:});

    pat_mask = strncmpi('p',seqNames,1);
    if any(pat_mask)
        for i = find(pat_mask)'
            pos_data{i} = unique(reshape(str2num(pos_data{i}),3,[])','rows');
        end
        data = [data; seqNames(pat_mask) miRNAnames(pat_mask) pos_data(pat_mask)];
    end
end
fclose(fid);

allmiRNAnames = unique(data(:,2));


patInds = regexpi(data(:,1),'Pat_ind_(\d*)_\d','tokens');
patInds = cat(1,patInds{:});
patInds = cat(1,patInds{:});

patInds = cellfun(@(x)(str2double(x)),patInds);

[junk miRNAinds] = ismember(data(:,2),allmiRNAnames);

pat_miRNA_data = NaN(length(PAT_STRUCT),length(allmiRNAnames));

for i = 1:size(pat_miRNA_data,1)
    for k = 1:size(pat_miRNA_data,2)
        spots = find(patInds==i&miRNAinds==k);
        if ~isempty(spots)
            temp = cat(1,data{spots,3});
            pat_miRNA_data(i,k) = min(temp(:,3),[],1);
        end
    end
end



[miRNAnamesmulti{1:length(PAT_STRUCT)}]=deal(allmiRNAnames');
miRNAnamesmulti=cat(1,miRNAnamesmulti{:});

OUT_STRUCT = PatientStructHelper(PAT_STRUCT,{'SimpleHumanMiRNA',1:length(PAT_STRUCT),pat_miRNA_data},{'HumanMiRNAnames',1:length(PAT_STRUCT),miRNAnamesmulti});















