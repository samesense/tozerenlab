require 'bio'

serv = Bio::KEGG::API.new
pathways = serv.list_pathways('hsa')
pathways.each do |path|
	      genes = serv.get_genes_by_pathway(path.entry_id)
	      genes.each do |gene|
	      		 puts gene + "\t" + path.entry_id
	      end
	      # puts path.entry_id + "\t" + path.definition
end

#genes = serv.get_genes_by_organism('hsa', 0, 40000)
## ls = serv.get_genes_by_pathway('path:eco00020')
#genes.each do |gene|
#	   paths = serv.get_pathways_by_genes([gene])
#	   puts gene, paths.length
#end