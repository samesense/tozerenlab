function OverlapDisplay(INPUT_DATA,varargin)
%   OverlapDisplay
%       Displays the Pair-wise Overlap data from CheckGeneListOverlap
%
%
%   OverlapDisplay(INPUT_DATA)
%
%
%   OverlapDisplay(INPUT_DATA,'background:type',BKG_DATA)
%       Uses the Background Data from random shuffling of the annotation
%       produced by CheckRandom


if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'background:SAM'
                BKG_DATA_AVAIL=true;
                BKG_DATA=varargin{i+1};
                ANAL_TYPE='SAM';
            case 'background:fold'
                BKG_DATA_AVAIL=true;
                BKG_DATA=varargin{i+1};
                ANAL_TYPE='fold';
            case 'background:ttest'
                BKG_DATA_AVAIL=true;
                BKG_DATA=varargin{i+1};
                ANAL_TYPE='ttest';
            otherwise
                warning(['Unknown Parameter' varargin{1}])
        end
    end
else
    BKG_DATA_AVAIL=false;
end


DISPLAY_FIGURE=figure;

good_disease_list=diag(INPUT_DATA.P_Overlap)~=-1;
good_disease_index=find(good_disease_list);
plotting_matrix=zeros(1+sum(good_disease_list));

if BKG_DATA_AVAIL
    bkg_nums=zeros(length(BKG_DATA(1).Fold_L),length(BKG_DATA));
    these_nums=cellfun(@length,INPUT_DATA.GeneOverlap(:));
    
    switch ANAL_TYPE
        case 'fold'
            for i=1:length(BKG_DATA)
                bkg_nums(:,i)=BKG_DATA(i).Fold_L;
            end
        case 'ttest'
            for i=1:length(BKG_DATA)
                bkg_nums(:,i)=BKG_DATA(i).TTest_L;
            end

        case 'SAM'
            for i=1:length(BKG_DATA)
                bkg_nums(:,i)=BKG_DATA(i).SAM_L;
            end
    end
    
    these_stats=reshape(sum(bkg_nums>(these_nums*ones(1,length(BKG_DATA))),2)/length(BKG_DATA),size(INPUT_DATA.P_Overlap));
    
    P_OVERLAP_FEATURES=INPUT_DATA.P_Overlap;
    P_OVERLAP_FEATURES(P_OVERLAP_FEATURES(:)==-1)=2*ones(1,sum(P_OVERLAP_FEATURES(:)==-1));
    
    P_OVERLAP_FEATURES(P_OVERLAP_FEATURES>=0&P_OVERLAP_FEATURES<=1)=((these_stats(P_OVERLAP_FEATURES>=0&P_OVERLAP_FEATURES<=1))<0.05);
    
    
else
    P_OVERLAP_FEATURES=INPUT_DATA.P_Overlap;
    P_OVERLAP_FEATURES(P_OVERLAP_FEATURES>=0&P_OVERLAP_FEATURES<=1)=(P_OVERLAP_FEATURES(P_OVERLAP_FEATURES>=0&P_OVERLAP_FEATURES<=1))<0.05;
    P_OVERLAP_FEATURES(P_OVERLAP_FEATURES(:)==-1)=2*ones(1,sum(P_OVERLAP_FEATURES(:)==-1));
end

plotting_matrix(1:end-1,1:end-1)=P_OVERLAP_FEATURES(good_disease_list,good_disease_list);

go_overlap=cellfun(@isempty,INPUT_DATA.GOOverlap(:));
go_overlap=reshape(~go_overlap,size(INPUT_DATA.GOOverlap));
go_overlap_list=INPUT_DATA.GOOverlap(good_disease_index,good_disease_index);
[go_I go_J]=find(go_overlap(good_disease_index,good_disease_index));

kegg_overlap=cellfun(@isempty,INPUT_DATA.KeggOverlap(:));
kegg_overlap_list=INPUT_DATA.KeggOverlap(good_disease_index,good_disease_index);
kegg_overlap=reshape(~kegg_overlap,size(INPUT_DATA.KeggOverlap));
[kegg_I kegg_J]=find(kegg_overlap(good_disease_index,good_disease_index));


hold on

pcolor(plotting_matrix)
scatter(go_I+.5,go_J+.5,'k','filled')
set(gca,'YTick',[1.5:length(INPUT_DATA.DiseaseNameVector(good_disease_list))+1],...
    'YTickLabel',INPUT_DATA.DiseaseNameVector(good_disease_list),'FontSize',10);
xticklabel_rotate([1.5:length(INPUT_DATA.DiseaseNameVector(good_disease_list))+1],90,INPUT_DATA.DiseaseNameVector(good_disease_list));
set(gca,'position',[0.2 0.3 0.6 0.6])
scatter(kegg_I+.5,kegg_J+.5,'k','filled')
colorbar('YTick',[0 1 2],'YTickLabel',{'Insignificant Overlap','Significant Overlap','No Overlap'})

hold off

while(1)
    [x y button]=ginput(1);
    if button==3
        break
    end

    display(INPUT_DATA.DiseaseNameVector{good_disease_index(floor(x))})
    display(INPUT_DATA.DiseaseNameVector{good_disease_index(floor(y))})


    display('Overlapping Kegg Pathways')
    display(kegg_overlap_list{floor(x),floor(y)}')
    display('Overlapping GO Terms')
    display(go_overlap_list{floor(x),floor(y)})

    %go_spot=find([ceil(x) ceil(y)]==[go_I go_J]);
    %go_kegg=find([ceil(x) ceil(y)]==[kegg_I kegg_J]);



    Articles=GetPUBMED([INPUT_DATA.DiseaseNameVector{good_disease_index(floor(x))} '+AND+' INPUT_DATA.DiseaseNameVector{good_disease_index(floor(y))}]);
    if isempty(Articles)
        display('No Articles Found')
    else
        for i=1:length(Articles)
            display(Articles(i).Title)
        end
    end

end
