function SIG_GENES=GetSigGenes(GE_DATA,ANOT_DATA,CHIP_DATA,ANAL_TYPE,varargin)
%   GetSigGenes
%       Get Significantly Altered Genes
%
%   SIG_GENES=GetSigGenes(GE_DATA,ANOT_DATA,ANAL_TYPE)
%
%       GE_DATA         Gene Expression Structure created from
%                       LoadAnnotatedGEData
%       ANOT_DATA       Probe Annotation Data created from
%                       LoadAffyNotation
%       CHIP_DATA       Affymatrix Chip Annotation Data created from
%                       LoadAffyNotation
%       ANAL_TYPE       The type of Analysis
%           'fold-change'   Determine Significantly Altered genes by a >1.5
%                           or <0.7 fold-change
%           't-test'        Determine Significantly Altered genes by a
%                           t-test value of a<0.001
%           'SAM'           Determine Significantly Altered genes by SAM
%                           analysis with 100-permutations and FDR of 1.0
%
%
%       SIG_GENES       Structure of Entrez-Gene IDS of significantly
%                       altered genes and information about the test.
%           SIG_GENES.Entrez `          Cell of Entrez Gene IDS for each
%                                       disease-state
%           SIG_GENES.GOMolFun          Significantly OverEnriched GO
%                                       Molecular Function Catagories
%           SIG_GENES.GOCellComp        Significantly OverEnriched GO
%                                       Cellular Compartment Catagories
%           SIG_GENES.GOBioProc         Significantly OverEnriched GO
%                                       Biological Process Catagories
%           SIG_GENES.KEGGPath          Significantly OverEnriched KEGG
%                                       Pathways
%           SIG_GENES.DiseaseNames      Cell of Disease Names
%           SIG_GENES.DataName          Dataset Name
%           SIG_GENES.AnalType          Type of Analysis
%           SIG_GENES.ChipType          Affymatrix Chip name
%
%
%   GetSigGenes(...,'Waitbar',WaitBarHandle,[start stop])
%
%       WaitBarHandle       The handle of a waitbar object
%       [start stop]        The proportion of the wait inwhich this
%                           instance should cover.
%
%

cutoff=0.01;   %hypergeometric significance cutoff
invalid_probes=cellfun(@isempty,ANOT_DATA(:,2));
valid_indexes=find(~invalid_probes);

if  isempty(varargin)
    h=waitbar(0,'');
    bar_prefix=[GE_DATA.DataName ': '];

    bar_start=0;
    bar_stop=1;
elseif strcmpi(varargin{1},'waitbar')
    h=varargin{2};
    bar_prefix=[get(h,'UserData') GE_DATA.DataName ': '] ;
    bar_start=varargin{3}(1);
    bar_stop=varargin{3}(2);
end


switch ANAL_TYPE
    case 'fold-change'
        normal=mean(GE_DATA.GEData(~invalid_probes,find(GE_DATA.DiseaseStates==0)),2);

        disease=zeros(sum(~invalid_probes),max(GE_DATA.DiseaseStates));
        if max(GE_DATA.DiseaseStates)>1
            for i=1:max(GE_DATA.DiseaseStates)
                disease(:,i)=mean(GE_DATA.GEData(~invalid_probes,GE_DATA.DiseaseStates==i),2);
            end
        else
            disease(:,1)=mean(GE_DATA.GEData(~invalid_probes,GE_DATA.DiseaseStates==1),2);
        end

        sig_gene_bool=(normal*ones(1,max(GE_DATA.DiseaseStates))*1.5)<disease...
            |(normal*ones(1,max(GE_DATA.DiseaseStates))*0.7)>disease;
    case 't-test'
        if max(GE_DATA.DiseaseStates)>1
            sig_gene_bool=false(sum(~invalid_probes),max(GE_DATA.DiseaseStates));
            for i=1:max(GE_DATA.DiseaseStates)
                sig_gene_bool(:,i)=ttest2(GE_DATA.GEData(~invalid_probes,GE_DATA.DiseaseStates==0),GE_DATA.GEData(~invalid_probes,GE_DATA.DiseaseStates==i),.001/sum(~invalid_probes),[],[],2)==1;
            end
        else
            sig_gene_bool=ttest2(GE_DATA.GEData(~invalid_probes,GE_DATA.DiseaseStates==0),GE_DATA.GEData(~invalid_probes,GE_DATA.DiseaseStates==1),.001/sum(~invalid_probes),[],[],2)==1;
        end
    case 'SAM'
        sig_gene_bool=false(sum(~invalid_probes),max(GE_DATA.DiseaseStates));
        if max(GE_DATA.DiseaseStates)>1
            for i=1:max(GE_DATA.DiseaseStates)

                %%%%%%%%%%%%%%%%Display Stuff%%%%%%%%%%%%
                set(h,'userdata',[GE_DATA.DiseaseNames{i} ': '])
                %END%%%%%%%%%%%%Display Stuff%%%%%%%%%%%%

                cols=GE_DATA.DiseaseStates==0|GE_DATA.DiseaseStates==i;
                genes=SAM_mod(GE_DATA.GEData(~invalid_probes,cols),GE_DATA.DiseaseStates(cols),100,1,...
                    'waitbar',h,[bar_start+.7*(bar_stop-bar_start)*(i/max(GE_DATA.DiseaseStates)) ...
                    bar_start+0.7*(bar_stop-bar_start)*(i+1)/max(GE_DATA.DiseaseStates)]);
                if ~isempty(genes)
                    sig_gene_bool(genes,i)=true(length(genes),1);
                end
            end
        else

            %%%%%%%%%%%%%%%%Display Stuff%%%%%%%%%%%%
            set(h,'userdata',[GE_DATA.DiseaseNames{1} ': '])
            %END%%%%%%%%%%%%Display Stuff%%%%%%%%%%%%


            genes=SAM_mod(GE_DATA.GEData(~invalid_probes,:),GE_DATA.DiseaseStates,100,1,...
                'waitbar',h,[bar_start bar_start+0.7*(bar_stop-bar_start)]);
            if ~isempty(genes)
                sig_gene_bool(genes)=true(length(genes),1);
            end
        end
end

%%%%%%%%%%%%CONVERT SIG-GENES TO ENTREZ IDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIG_GENES.Entrez=cell(1,max(GE_DATA.DiseaseStates));
SIG_GENES.KeggPath=cell(1,max(GE_DATA.DiseaseStates));
SIG_GENES.GONames=cell(1,max(GE_DATA.DiseaseStates));
SIG_GENES.GeneNames=cell(1,max(GE_DATA.DiseaseStates));

waitbar(bar_start+0.7*(bar_stop-bar_start),h,[bar_prefix 'Collating Genes'])
if max(GE_DATA.DiseaseStates)>1
    for i=1:max(GE_DATA.DiseaseStates)
        if sum(sig_gene_bool(:,i))~=0
            SIG_GENES.Entrez{i}=ANOT_DATA(valid_indexes(sig_gene_bool(:,i)),2);
            SIG_GENES.KeggPath{i}=CHIP_DATA.KeggNames(get_HyperTest(i,'KeggArray'));
            SIG_GENES.GONames{i}=CHIP_DATA.GONames(get_HyperTest(i,'GOArray'));
            SIG_GENES.GeneNames{i}=ANOT_DATA{valid_indexes(sig_gene_bool(:,i)),5};
        end
    end
else
    if sum(sig_gene_bool(:,1))~=0

        SIG_GENES.Entrez{1}=ANOT_DATA(valid_indexes(sig_gene_bool(:,1)),2);
        SIG_GENES.KeggPath{1}=CHIP_DATA.KeggNames(get_HyperTest(1,'KeggArray'));
        SIG_GENES.GONames{1}=CHIP_DATA.GONames(get_HyperTest(1,'GOArray'));
        SIG_GENES.GeneNames{1}=ANOT_DATA(valid_indexes(sig_gene_bool(:,1)),5);
    end
end

SIG_GENES.DiseaseNames=GE_DATA.DiseaseNames';
SIG_GENES.DataName=GE_DATA.DataName;
SIG_GENES.AnalType=ANAL_TYPE;
SIG_GENES.ChipType=GE_DATA.ChipType;

waitbar(bar_stop,h,[bar_prefix 'Finished'])

    function entrez_ids=get_ids(ID)
        entrez_ids=[ANOT_DATA{strcmp(ANOT_DATA(:,1),ID),2}]';
    end

    function sig_indexes=get_HyperTest(row,field_name)
        input=getfield(CHIP_DATA,field_name);
        defective_picked=sum(input(valid_indexes(sig_gene_bool(:,row),:)));
        total_lot_size=size(input(valid_indexes,:),1);
        num_defective=sum(input(valid_indexes,:),1);
        num_picked=sum(sig_gene_bool(:,row));
        data=hygecdf(defective_picked,total_lot_size,num_defective,num_picked);
        sig_indexes=(data<cutoff);
    end


end