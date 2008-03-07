function SIGNATURE_STRUCT=LoadGSEASignatures(DIRECTORY)
%   LoadGSEASignatures
%       Creates a struct-array containing the gene-lists compatible with
%       the matlab GSEA analysis program.
%
%
%   SIGNATURE_STRUCT=LoadGSEASignatures(DIRECTORY)
%       Loads all signatures in DIRECTORY.
%





SIGNATURE_STRUCT=[];

sig_files=dir([DIRECTORY '\*.GMT']);

for i=1:length(sig_files)
    fid=fopen([DIRECTORY '\' sig_files(1).name],'rt');

    data=textscan(fid,'%s na %[^\n]','bufsize',10000);

    gene_symbols=regexp(data{2},'\w+','match');

    temp_struct=cell2struct([data{1} gene_symbols],{'Name','GeneSymbols'},2);

    [temp_struct.GroupName]=deal(sig_files(i).name(1:end-3));

    SIGNATURE_STRUCT=[SIGNATURE_STRUCT; temp_struct];
end


% fid=fopen([DIRECTORY '\GO\gene_ontology.obo'],'rt');
% 
% %get rid of header info
% while strcmp(junk,'[Term]')
%     junk=fgetl(fid);
% end
% 
% GO_ids=struct('GO_ID',num2cell((1:65010)',2),'Name',[],'Type',[],'Inheritence',[]);
% 
% while 1
%     junk=fgetl(fid);
%     if junk==-1
%         break
%     end
%     switch [junk(1:find(junk==':',1)) 'a']
%         case {'a'}
%             pointer=[];
%         case {'id:a'}
%             pointer=str2num(junk(find(junk==':',1,'last')+1:end));
%         case {'name:a'}
%             GO_ids(pointer).Name=strtrim(junk(find(junk==':',1,'last')+1:end));
%         case {'name_space:a'}
%             GO_ids(pointer).Type=
%         otherwise
%             display('nope')
%         
%         
%         
%         
%     end
%     
%     
%     
% end
% 








%end