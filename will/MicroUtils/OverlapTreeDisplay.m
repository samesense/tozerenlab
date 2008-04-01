function OverlapTreeDisplay(INPUT_DATA,varargin)

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





good_disease_list=diag(INPUT_DATA.P_Overlap)~=-1;
good_disease_index=find(good_disease_list);

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
    sig_mask=these_stats>=0&these_stats<=0.05;

else
    P_OVERLAP_FEATURES=INPUT_DATA.P_Overlap;
    sig_mask=(P_OVERLAP_FEATURES>=0&P_OVERLAP_FEATURES<=0.05);
end

gene_dist_vec=cellfun(@length,INPUT_DATA.GeneOverlap(:));
gene_dist_mat=1./(reshape(gene_dist_vec,size(INPUT_DATA.GeneOverlap))+1);
gene_dist_mat(~sig_mask)=1;
gene_tree=seqlinkage(gene_dist_mat(good_disease_list,good_disease_list),'average',INPUT_DATA.DiseaseNameVector(good_disease_list));

GO_dist_vec=cellfun(@length,INPUT_DATA.GOOverlap(:));
GO_dist_mat=1./(reshape(GO_dist_vec,size(INPUT_DATA.GeneOverlap))+1);
GO_tree=seqlinkage(GO_dist_mat(good_disease_list,good_disease_list),'average',INPUT_DATA.DiseaseNameVector(good_disease_list));

Kegg_dist_vec=cellfun(@length,INPUT_DATA.KeggOverlap(:));
Kegg_dist_mat=1./(reshape(Kegg_dist_vec,size(INPUT_DATA.GeneOverlap))+1);
Kegg_tree=seqlinkage(Kegg_dist_mat(good_disease_list,good_disease_list),'average',INPUT_DATA.DiseaseNameVector(good_disease_list));

plot(gene_tree)
title('Linkage By Gene Overlap')
plot(GO_tree)
title('Linkage By Gene Ontology Overlap')
plot(Kegg_tree)
title('Linkage By Kegg Overlap')

level_distance=pdist(gene_tree,'criteria','levels','squareform',true);

for i=1:length(level_distance)
    genes=cell2mat(INPUT_DATA.GeneOverlap(good_disease_index(i),good_disease_index(i)));
    close_spots=find(level_distance(:,i)<=3);
    arrayfun(@MultiIntersect,close_spots)

end

    function MultiIntersect(input)
        genes=intersect(genes,cell2mat(INPUT_DATA.GeneOverlap(good_disease_index(input),good_disease_index(i))));
    end
end