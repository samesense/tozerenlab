function HIV_REF=HIVrefLoader(FILENAME,SubInput,SkipBads)
%   HIVrefLoader
%       HIVrefLoader reads a FASTA formatted file which has the complete
%       HIV genomic region followed by the the NT sequence of each protien
%       region.  This produces and HIV_REF which is compatible with all of
%       the GenerateRead functions.
%
%   HIV_REF=HIVrefLoader(FILENAME)
%
%

if iscell(FILENAME)&&length(FILENAME)>1
    HIV_REF = [];
    
    for i = 1:length(FILENAME)
        HIV_REF = [HIV_REF; HIVrefLoader(FILENAME{i},SubInput{i},SkipBads)];
    end
    
    
else
    testData = genbankread(FILENAME);
    [gene_names{1:length(testData.CDS)}] = deal(testData.CDS(:).gene);
    [product_names{1:length(testData.CDS)}] = deal(testData.CDS(:).product);
    [notes{1:length(testData.CDS)}] = deal(testData.CDS(:).note);
    
    [standardizedNames cannotDefine] = nameMapping(gene_names,product_names,notes);
    
    if SkipBads&&any(cannotDefine)
        HIV_REF=[];
        return
    end
    
    [trans_aa{1:length(testData.CDS)}] = deal(testData.CDS(:).translation);
    [locs{1:length(testData.CDS)}] = deal(testData.CDS(:).indices);
    locsval = cell2mat(cellfun(@(x)(x(1:2)),locs','uniformoutput',false));
    HIV_REF=struct('Subtype',SubInput,'name',testData.LocusName,...
        'Sequence',upper(testData.Sequence),...
        'GeneNames',{standardizedNames(~cannotDefine)'},'AAseqs',{trans_aa(~cannotDefine)'},...
        'GenePos',locsval(~cannotDefine,:),'AppendedGenome',[trans_aa{~cannotDefine}],...
        'AAPos',cumsum([1 cellfun('length',trans_aa(~cannotDefine(1:end-1)))]));
end

    function [standardNames junked] = nameMapping(genes_annot,product_annot,notes_annot)
        standardNames = cell(size(genes_annot));
        junked = false(size(genes_annot));
        
        known_names{1} = 'gag';
        known_names{2} = 'pol[^y]';
        known_names{3} = 'vif';
        known_names{4} = 'vpr';
        known_names{5} = 'tat';
        known_names{6} = 'rev';
        known_names{7} = 'vpu';
        known_names{8} = 'env';
        known_names{9} = 'nef';
        
        strs = cat(2,genes_annot',product_annot');
        strs(cellfun('isempty',strs)) = {' '};
               
        possible_choices = regexpi(strs,'\w*','match');

        for i = 1:size(possible_choices,1)
            temp = cat(2,possible_choices{i,:});
            if isempty(temp)
                continue
            end
            tempStrs = cell(1,2*size(temp,2));
            tempStrs(:) = {' '};
            tempStrs(1:2:end) = temp(:);
           
            pos=regexpi(cat(2,tempStrs{:}),known_names,'start');
            
            if ~all(cellfun('isempty',pos))&&sum(~cellfun('isempty',pos))==1
                standardNames{i}=known_names{find(~cellfun('isempty',pos),1)}(1:3);
            end
        end
        
        %second pass
        
        strs = cat(2,genes_annot',product_annot',notes_annot');
        strs(cellfun('isempty',strs)) = {' '};
        strs(:,3) = cellfun(@(x)(x(1,:)),strs(:,3),'uniformoutput',false);
                
        possible_choices = regexpi(strs,'\w*','match');
        
        for i = 1:size(possible_choices,1)
            if isempty(standardNames{i})
                temp = cat(2,possible_choices{i,:});
                tempStrs = cell(1,2*size(temp,2));
                tempStrs(:) = {' '};
                tempStrs(1:2:end) = temp(:);

                pos=regexpi(cat(2,tempStrs{:}),known_names,'start');

                if all(cellfun('isempty',pos))
                    display('Could not find any in:')
                    display(cat(2,tempStrs{:}))
                    junked(i)=true;
                elseif sum(~cellfun('isempty',pos))>1
                    display('Found multiple in:')
                    display(cat(2,tempStrs{:}))
                    display(cat(2,known_names{~cellfun('isempty',pos)}))
                    junked(i)=true;
                else

                    standardNames{i}=known_names{find(~cellfun('isempty',pos),1)}(1:3);
                end
            end
        end

        
    end
end