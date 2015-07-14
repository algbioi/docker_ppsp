#!/usr/bin/ruby -w
#calculate binning fragmentation 
#its assumed that there are scaffolds and contigs and prediction per contig
#measure how the binning is fragmented for each scaffold and then average it
#newick file minimizes calls to the biosql database
#the contigs prediction file is suppused to be in the following format (tab-separated)
#fasta_header prediction_1 prediciton_2 ... prediction_n
#where prediction_n is the most specific prediciton
#the fasta_header must have Scaffold_m_ at beginning where m is the scaffold identifier

#IF SOME CONTIG IS NOT PREDICTED THEN IT IS IGNORED IN CONSISTENCY CALCULATION
#BUT IT IS PRINTED IN THE GRAPH FILES


mydir = File.dirname(__FILE__)

require mydir+"/lib/rubygems/trollop/trollop.rb"

#parse arguments
opts = Trollop::options do
  
  version "Scaffold-contig consistency 1.0.0 (c) 2010 Kaustubh Patil\n* are necessary options"

  opt :pred, "Predictions file (tab delimited in PP format, only the first and last columns are considered)*", :default => nil, :type => String
  opt :newick, "Newick file*", :default => nil, :type => String
  opt :fasta, "FASTA file (headers should be contid ids)*", :default => nil, :type => String
  opt :sqlite, "SQLite database file contaning the taxonomy information*", :default => nil, :type => String
  opt :scafcont, "Scaffold-contig mapping file (tab delimited)*", :default => nil, :type => String
  opt :graph, "Generate graph files (default: false)", :default => false
  opt :graphf, "Graph file prefix (default: predictions file)", :default => nil, :type => String
  opt :pop, "Population (necessary for creating graph files)", :default => nil, :type => String
  opt :pop_restrict, "Restrict analysis to scaffolds with at-least one pop contig", :default => true 
  opt :min_n_cont, "Minimum number of contigs inside scaffold (0 means no limit)", :default => 0
  opt :min_l_scaf, "Minimum scaffold length (0 means no limit)", :default => 0
  opt :min_l_cont, "Minimum contig length (0 means no limit)", :default => 0
  opt :max_l_cont, "Maximum contig length (0 means no limit)", :default => 0
  opt :tree_extend, "Maximum contig length", :default => true
  opt :ncbi_taxon_id, "Whether the predictions are NCBI taxonomic ids (default false)", :default => false
  opt :n_levels, "Number of consistent taxonomic levels", :default => nil, :type => :int
  opt :missing_scaf_id, "Missing scaffold id (0: ignore, 1: add, 2: error)", :default => 0
  opt :sort_cont, "Sort contigs in a scaffold using identifiers", :default => true
  opt :ranks, "Taxonomic ranks (comma separated, no space)", :default => "species,genus,family,order,class,phylum,superkingdom,root"
  opt :unit, "Unit for lengths in smmary (bp/kb/mb)", :default => "kb"
  opt :warn, "Print warnings on STDERR", :default => true
end

#check for sanity of arguments
Trollop::die :pred, "must exist" if opts[:pred].nil? || !File.exist?(opts[:pred])
Trollop::die :newick, "must exist" if opts[:newick].nil? || !File.exist?(opts[:newick])
Trollop::die :fasta, "must exist" if opts[:fasta].nil? || !File.exist?(opts[:fasta])
Trollop::die :scafcont, "must exist" if opts[:scafcont].nil? || !File.exist?(opts[:scafcont])
Trollop::die :sqlite, "must exist" if opts[:sqlite].nil? || !File.exist?(opts[:sqlite])

Trollop::die :unit, "must be bp/kb/mb" if !["bp","kb","mb"].member? opts[:unit]

raise "Please provide population (--pop) if you want to generate graph files" if opts[:graph] && opts[:pop].nil?

require mydir+"/lib/sqlite_utils.rb"
t = SqliteUtils::Taxonomy.new(opts[:sqlite])

#should it create graph files? ONLY USE FOR SINGLE POPULATION OTHERWISE WONT WORK
#files names are automatically created using the "contigs_prediction_file"
#graph files contain two files
#1. Co-ordinate matrix file (row col value)
#2. color code file
if opts[:graph]
	opts[:graphf] = opts[:pred] if opts[:graphf].nil?
	graph_mat_file = File.open(opts[:graphf]+"_graph_mat.txt","w")
	graph_col_file = File.open(opts[:graphf]+"_graph_col.txt", "w")
	graph_rname_file = File.open(opts[:graphf]+"_graph_rname.txt", "w")
	graph_legend_file = File.open(opts[:graphf]+"_graph_legend.txt", "w")
end

#define colors
graph_col = {
	#if something is not predicted
	#shade of blue
	"na" => "109",

	#shades of green
	"correct" => "518",

	"same_genus"   => "518",
	"same_family"  => "517",
	"same_order"   => "516",
	"same_class"   => "515",
	"same_phylum"  => "514",
	"same_domain"  => "513", #tis is yellowish
	"same_superkingdom"  => "513", #ncbi calls domain superkingdom
        
	#shades of red
	"different_genus"   => "501",
        "different_family"  => "503",
        "different_order"   => "504",
        "different_class"   => "505",
        "different_phylum"  => "506",
        "different_domain"  => "507",
	"different_superkingdom"  => "507",
	
	#grey
	"root"    => "220",
	"same_root"    => "220",
	"different_root"    => "220"
	}

#ranks
ranks = opts[:ranks].split(",")

#if the fasta file is given then get header and length of each contig
contig_lengths = Hash.new(0)
if !opts[:fasta].nil?
	fasta = Bio::FastaFormat.open(opts[:fasta])
	fasta.each do |seq|
		head = seq.definition
		raise "Nil header" if head.nil?
		len = seq.length
		raise "No sequence for: #{head}" if len.nil? || len==0
		contig_lengths[head] = len
		#puts "#{head}:#{len}"
	end

	fasta = TRUE #used for printing information
end

#get scaffolds and contigs
scaffolds = Hash.new(nil)       #scaffold id for each contig, here the scaffold part is taken away
contigs = Hash.new{|hash, key| hash[key] = Array.new;} #array of headers for each scaffold id
n_bp = Hash.new(0) # number of bp of a scaffold
n_non_scaff = -1 #if some contigs havent been assigned to any scaffold then treat them as their own scaffold
fp =  File.open(opts[:scafcont], "r")
n_line=0
while(line=fp.gets)

        n_line += 1

        next if line.nil?
        next if line[0].chr=="#"

        line.chomp!
        next if line.empty?
        line = line.split(/\t/)

	raise "Not two entries at line #{n_line} of scafcont file!" if line.length!=2

        #get scaffold id
        if line[0].nil? ||  line[0].empty? ||  line[0]==" "
                if opts[:missing_scaf_id]==0
                        scaf = nil
                elsif opts[:missing_scaf_id]==1
                        scaf = n_non_scaff.to_s
                        n_non_scaff = n_non_scaff-1
                else
                        raise "\nNo scaffold id at line #{n_line}\n\n"
                end
        else
                scaf = line[0]
        end

	cont = line[1]

	scaffolds[cont] = scaf

	next if scaf.nil? #dont count things for missing scaffolds

        contigs[scaf] << cont
        n_bp[scaf] += contig_lengths[cont]
end
fp.close

#get the predictions for each contig
fp = File.open(opts[:pred], "r")

#get predictions for each contig in the preds hash

preds = Hash.new(nil) 		#prediction for each individual contig
myrank = Hash.new(nil)
mytaxid = Hash.new(nil)
myname = Hash.new(nil)
n_line = 0
while(line=fp.gets)

	n_line += 1

	next if line.nil?
	next if line[0].chr=="#"

	line.chomp!
	next if line.empty?
	line = line.split(/\t/)

	if line.length<=1
		STDERR.puts "#WARN: no predicitons for #{line.to_s} (at line #{n_line})" if opts[:warn]
		#next #commenting this so that no prediction contigs will be used for analysis
	end

	#get contig and scaffold id for this prediction
	cont = line[0]
	scaf = scaffolds[cont]

	next if scaf.nil? #this is to ignore missing scaffolds

	line[0] = ""

	#get prediction
	#which is the last field of a line
	pred = line[1]
	line.each do |x|
		next if x.empty? || x==" " || x=="	"
		pred=x
	end

	if opts[:ncbi_taxon_id] && pred.to_i<=0
		STDERR.puts "#WARN: no predicitons for #{line.to_s} (at line #{n_line})" if opts[:warn]
		pred = "1"
	end

	preds[cont] = pred if line.length>1

	if myrank[pred].nil? && !pred.nil?
		tax_id = pred
		tax_id = t.tax_id(pred) if !opts[:ncbi_taxon_id]
		myrank[pred] = t.rank(tax_id)
		mytaxid[pred] = tax_id
		myname[tax_id] = pred
	end

	raise "Line #{n_line}: nil rank for #{pred}" if myrank[pred].nil?
	raise "Line #{n_line}: nil taxid for #{pred}" if mytaxid[pred].nil?
	

end
fp.close

#check if population exists and quit if not
raise "The population #{opts[:pop]} doesn't exits in predictions!" if opts[:graph] && !preds.member?(opts[:pop])

#select population scaffolds to process
if !opts[:pop].nil? && opts[:pop_restrict]
	scaff_to_delete = []
	contigs.each do |scaf,conts|
		mypreds = []
		conts.each do |cont|
			mypreds << preds[cont]
		end
		scaff_to_delete << scaf if !mypreds.member? opts[:pop]
	end
	#now delete them
	scaff_to_delete.each {|scaf| contigs.delete(scaf)}	
end

#sort contigs if asked
if opts[:sort_contigs]
	contigs.each {|s,c| contigs[s] = c.sort}
end

#get the tree
fp = File.open(opts[:newick], "r")
line = fp.gets
line = line.chomp
fp.close

#by default the parse takes internal (numeric) node names as bootstrap values so disable it
opt=Hash.new
opt[:bootstrap_style] = :disabled 
tree = Bio::Newick.new(line, options=opt).tree
#tree.each_node {|n| puts n}

#collect names and ranks of the nodes in the tree (not everything was observed in the predictions)
tree.nodes.each do |nd|
	tax_id = nd.name
	myname[tax_id] = t.scientific_name(tax_id)	if  myname[tax_id].nil?
	name = myname[tax_id]
	myrank[name] = t.rank(tax_id) if myrank[name].nil?
	raise "\nNo rank found for #{name}\n\n" if myrank[name].nil?
	mytaxid[name] = tax_id
end

#put some information and interpretation
puts "#=======================================================================#"
puts "#Scaffold-contig assignment consistency"
puts "#Author (c)   : Kaustubh Patil (patil@mpi-inf.mpg.de)"
puts "#Generated    : " + Time.new.inspect
puts "#Input file   : #{opts[:pred]}"
puts "#Tree file    : #{opts[:newick]}"
puts "#Sequence file: #{opts[:fasta]}"
puts "#Each line summarizes the consistency for one scaffold with the associated contigs"
#puts "#scaff info   : Scaffold information in terms of number of contigs and their lengths"
puts "#Each contig is annotated with one of the following annotations."
puts "#               * = assignment to scaffold clade"
puts "#               + = taxonomically consistent higher rank assignment"
puts "#               - = taxonomically inconsistent assignment"
puts "#Clade statistics at the end of the file."
puts "" #empty line
puts "#Scaffold_id Scaffold_assignment Contig_assignments"

dist_total = 0
n_scaffolds = 0
n_contigs = 0
n_total_bp = 0
dist_most_spec_pred = Hash.new(0)    #average distance for this prediction
n_scaf_most_spec_pred = Hash.new(0) #how many scaffolds have this most specific prediction
n_cont_most_spec_pred = Hash.new(0) #how many contigs have this most specific prediction
n_bp_most_spec_pred = Hash.new(0)   #base-pairs
n_cont = Hash.new(0)	#how many total contigs for this most specific prediction
path = Hash.new{|hash, key| hash[key] = Array.new;} #path of a most specific prediction
n_in_path = Hash.new(0) #how many contigs have prediction in the path of the most specific prediction
n_bp_in_path = Hash.new(0)
graph_legend = Hash.new(nil)
contigs.each do |scaffold_id, contig_ids|

	n_scaffolds += 1

	#gather information for this scaffold
	n_contigs_this_scaffold = contig_ids.length
	#get lengths of contigs in this scaffold if fasta file was provided
	contig_lengths_this_scaffold = []
	scaffold_length = nil
	if fasta
		contig_ids.each { |c|  contig_lengths_this_scaffold << contig_lengths[c] }
		scaffold_length = n_bp[scaffold_id]
	end

	#skip scaffold if asked for
	next if opts[:min_l_scaf] > 0 && scaffold_length < opts[:min_l_scaf]
	next if opts[:min_n_cont] > 0 && n_contigs_this_scaffold < opts[:min_n_cont]

	#get most specific prediction
	most_spec_pred_index = ranks.length+1
	most_spec_pred = "root"

	if !opts[:pop].nil? && !opts[:pop].empty?
		most_spec_pred = opts[:pop]
	else
		#get it from data
		contig_ids.each do |ci|
			v = preds[ci] #this is the prediction for this contig
			if v.nil? || v.empty?
				STDERR.puts "#WARN: Nil prediction for contig #{ci}" if opts[:warn]
				next
			end

			#see if this prediction is in the tree
			node = tree.get_node_by_name(mytaxid[v])
			if node.nil? && !opts[:tree_extend]
				#STDERR.print "Finding in tree #{v} (#{mytaxid[v]})"
				#find in the tree
				ranks.each do |r|
					parent = t.parent_atrank(mytaxid[v], r)
					next if parent.nil?
					node = tree.get_node_by_name(parent)
					next if node.nil?
					myname[parent] = t.scientific_name(parent) if myname[parent].nil?
					v = myname[parent]
					preds[ci] = v
					break
				end
				raise "Couldn't map #{v} in the tree" if node.nil?
				#STDERR.print "...mapped to #{v} (#{mytaxid[v]})\n"
			end

			myrank[v] = t.rank(mytaxid[v]) if myrank[v].nil?
			s = ranks.index(myrank[v])
			raise "\nCouldn't find #{v} (#{myrank[v]}) in the rank array (#{s}) or is nil\n" if myrank[v].nil? || s.nil? 
			if s<most_spec_pred_index
				most_spec_pred_index = s
				most_spec_pred = v
			end

		end

		#if there are multiple most spec predictions and 
		#fasta file is given then assign the longest one
		most_spec_pred_lens = Hash.new(0)
		if fasta
			#get combined lengths for all most spec preds
			(0...contig_ids.length).each do |iv|
				c = contig_ids[iv]
				v = preds[c]
				s = ranks.index(myrank[v])
				most_spec_pred_lens[v] = most_spec_pred_lens[v]+contig_lengths[c] if s==most_spec_pred_index
			end
			#now assign the longest one
			most_spec_pred = nil
			most_spec_pred_len = 0
			most_spec_pred_lens.each do |msp, mspl|
				if most_spec_pred.nil?
					most_spec_pred = msp
					most_spec_pred_len = mspl
				elsif most_spec_pred_len<mspl
					most_spec_pred = msp
					most_spec_pred_len = mspl
				end
			end
		end
	end

	next if most_spec_pred.nil?

	#count how many have this most specific prediction
	n_scaf_most_spec_pred[most_spec_pred] += 1

	#the newcik tree is supposed to contain ncbi tax ids
	most_spec_node = tree.get_node_by_name(mytaxid[most_spec_pred])

	#
	#now we have the most specific prediction in this scaffold in most_spec_pred & most_spec_node
	#
	
	#add this path to the tree if it doesn't exist and user asks for
	if most_spec_node.nil? && opts[:tree_extend]
		child = most_spec_pred
		child_node = Bio::Tree::Node.new(mytaxid[child])
		tree.add_node(child_node)
		parent_node = nil
		ranks.each do |r|
			parent = t.parent_atrank(mytaxid[child], r)
			if !parent.nil?
				parent_node = tree.get_node_by_name(parent)
				parent_node = Bio::Tree::Node.new(parent) if parent_node.nil?
				tree.add_node(parent_node)
				tree.add_edge(child_node, parent_node)
				mytaxid[parent]=parent
				child = parent
				child_node = parent_node
			end
		end
		#add connection to root
		tree.add_edge(tree.root, child_node)
		#get more info about this node
		tax_id = mytaxid[most_spec_node]
		myname[tax_id] = most_spec_pred
		myrank[most_spec_pred] = t.rank(tax_id) if myrank[most_spec_pred].nil?		
	end

	most_spec_node = tree.get_node_by_name(mytaxid[most_spec_pred])
	raise "\nMost specific node not in the tree #{most_spec_pred} #{mytaxid[most_spec_pred]}\n" if most_spec_node.nil? && !most_spec_pred.nil?

	#get distance and statistics w.r.t to the most_spec_pred
	dist = 0
	n_non_spec_pred = 0

	#get the path to be considered as consistent
	path = []
	par_node = nil
	if opts[:n_levels].nil?	#get the consistency till root node
		par_node = tree.root
	elsif opts[:n_levels]==0	#get the consistency only for most_spec_pred
		par_node = most_spec_node
	else				#get the consistency untill asked ancestrors
		par_node = most_spec_node
		n_levels_traced = n_levels
		while n_levels_traced!=0 do
			par_node = tree.parent(par_node)
			n_levels_traced-=1
			break if par_node==tree.root
		end
	end
	
	raise "\nParent node is nil for prediction #{v}\n\n" if par_node.nil?
	
	p = tree.path(most_spec_node, par_node)
	p.each {|x| path << x.name}

	annotated_pred = [] #annotated consistent/inconcistent predictions
	(0...contig_ids.length).each do |i|
		
		cont_id_current = contig_ids[i] #get the corresponding contig_id
		cont_len_current = contig_lengths[cont_id_current]
		v = preds[cont_id_current] #v is the prediction now (this will be name for PP type file)
		
		if v.nil?
			tax_id_current = "0"
			rank_current = "na"
		else
			tax_id_current = mytaxid[v]

			raise "\nNo tax_id found for #{v}\n\n" if tax_id_current.nil?
			rank_current = myrank[v]
			raise "\nNo rank found for #{v}\n\n" if rank_current.nil?
		end

		#skip the contig if asked for
		next if opts[:min_l_cont] > 0 && cont_len_current < opts[:min_l_cont]
		next if opts[:max_l_cont] > 0 && cont_len_current > opts[:max_l_cont]
		
		#get overall statistics
		if !v.nil?
			n_bp_most_spec_pred[most_spec_pred] += contig_lengths[cont_id_current]
        	        n_total_bp += contig_lengths[cont_id_current]
			n_cont[most_spec_pred] += 1
		end
	
		n_contigs += 1 #increment this here as its used for graph matrix
		
		#get node for this prediction
		node = tree.get_node_by_name(tax_id_current)

		#raise "\nNode #{mytaxid[v]} not found\n" if node.nil?
		if node.nil? && tax_id_current != "0" #add the node to tree
			child = v
			child_node = Bio::Tree::Node.new(mytaxid[child])
			tree.add_node(child_node)
			parent_node = nil
			ranks.each do |r|
				 parent = t.parent_atrank(mytaxid[child], r)
				 if !parent.nil?
				 	parent_node = tree.get_node_by_name(parent)
					parent_node = Bio::Tree::Node.new(parent) if parent_node.nil?
					tree.add_node(parent_node)
					tree.add_edge(child_node, parent_node)
					mytaxid[parent]=parent
					child = parent
					child_node = parent_node
				 end
			end
			#add connection to root
			tree.add_edge(tree.root, child_node)
			tax_id = mytaxid[v]
			myname[tax_id] = v
			myrank[v] = t.rank(tax_id) if myrank[v].nil?		
		end
		node = tree.get_node_by_name(tax_id_current)
		raise "\nNode #{mytaxid[v]} not found\n" if node.nil? && tax_id_current != "0"
		
		this_prediction_is = ""
		if v==most_spec_pred
			n_cont_most_spec_pred[most_spec_pred] += 1
			n_bp_in_path[most_spec_pred] += contig_lengths[contigs[scaffold_id][i]]
			annotated_pred << v+"*"
			this_prediction_is = "correct"
		elsif path.member?(tax_id_current)
			n_non_spec_pred += 1
			n_in_path[most_spec_pred] +=1
			annotated_pred << v+"+"
			n_bp_in_path[most_spec_pred] += contig_lengths[contigs[scaffold_id][i]]
			this_prediction_is = "consistent"
		elsif v.nil?
			n_non_spec_pred += 1
			annotated_pred << "Not predicted-"
			this_prediction_is = ""
		else
			n_non_spec_pred += 1
			annotated_pred << v+"-"
			this_prediction_is = "inconsistent"
		end		
	
		if !v.nil?	
			dist += (tree.path(node, most_spec_node).length-1).to_f #path contains the nodes themselves
		end

		#print stuff in the graph files
		if opts[:graph]		

			#first the color
			lca = tree.root
			lca_rank = "same_root"
			col_to_use = nil
			legend = nil #used for the color
			if this_prediction_is.empty?
				col_to_use = graph_col["na"]
				legend = "Not predicted"
			elsif this_prediction_is=="correct"
				col_to_use = graph_col["correct"]
				lca = node
				legend = v
			elsif this_prediction_is=="consistent"
				col_to_use = graph_col["same_"+rank_current]
				legend = v
			elsif this_prediction_is=="inconsistent"
				lca = tree.lowest_common_ancestor(node, most_spec_node)
				raise "\nCouldn't find LCA\n\n" if lca.nil?
				#lca.name is tax_id from the tree and lca_name is the prediction name
				lca_name = myname[lca.name]
				lca_rank = myrank[lca_name]
				col_to_use = graph_col["different_"+lca_rank]
				legend = "Inconsistent below "+lca_rank
			end
			#if col_to_use.nil?
			#	col_to_use = "black"
			#end
			raise "No color found" if col_to_use.nil?
			graph_col_file.puts "#{col_to_use}"

			#now the bar_info
			bar_len = contig_lengths[cont_id_current]
			graph_mat_file.puts "#{n_scaffolds}\t#{n_contigs}\t#{bar_len}"

			#legend info
			#just store it now print it all together
			raise "\nMessed up legend for #{legend}: #{graph_legend[v]} (#{col_to_use}) [!check if a label is forced!]\n\n" if !graph_legend[legend].nil? && graph_legend[legend]!=col_to_use
			graph_legend[legend] = col_to_use	
		end
		
	end #end processing of constituent contigs

	graph_rname_file.puts scaffold_id if opts[:graph]

	dist /= contig_ids.length.to_f
	dist_total += dist
	dist_most_spec_pred[most_spec_pred] += dist
	dist = "%.3f" % dist.to_s	

	#print it
	print scaffold_id+"\t"+most_spec_pred+"\t"
	(0...n_contigs_this_scaffold).each do |i|
		print "#{contig_ids[i]}(#{annotated_pred[i]}"
		print ",#{contig_lengths_this_scaffold[i]}bp" if fasta
		print ")"
		print "," if i<(n_contigs_this_scaffold-1)	
	end
	#print "\t#{dist}"	
	print "\n"
end

#print statistics
#

unit = 1.0
unit = 1000.0 if opts[:unit]=="kb"
unit = 1000000.0 if opts[:unit]=="mb"


n_total_bp_unit = n_total_bp.to_f/unit
avg_contig_consistency = 0.0
avg_bp_consistency = 0.0
n_clades = 0
puts "@===============Clade statistics===============@"
puts "@Average distance:\t#{dist_total.to_f/n_contigs.to_f}"
puts "@Total scaffolds:\t#{n_scaffolds}"
puts "@Total contigs:\t#{n_contigs}"
puts "@Total base-pairs:\t#{n_total_bp}bp (#{n_total_bp_unit}#{opts[:unit]})"
puts "@==============================================@"
n_scaf_most_spec_pred.each do |key, val|

	#here key is most spec prediction and val is number
	
	percent_scaffods = 100*val.to_f/n_scaffolds.to_f
	percent_contigs =  100*n_cont[key].to_f/n_contigs.to_f
	percent_contigs_most_spec = 100*n_cont_most_spec_pred[key].to_f/n_cont[key].to_f

	percent_consistent_pred = 100*n_in_path[key].to_f/n_cont[key].to_f 
	percent_consistent_pred = 100 if percent_consistent_pred.nan?
	percent_consistent_pred = percent_consistent_pred+percent_contigs_most_spec

	average_distance = dist_most_spec_pred[key].to_f/n_cont[key].to_f

	#calculate percentage lengths	
	percent_bp = 100*n_bp_most_spec_pred[key].to_f/n_total_bp.to_f
	percent_contig_bp = 100*(n_bp_in_path[key].to_f) / n_bp_most_spec_pred[key].to_f

	#print it
	#puts "@==========================================@"

	#change to unit
	n_bp_most_spec_pred_unit = n_bp_most_spec_pred[key].to_f/unit
	n_bp_in_path_unit = n_bp_in_path[key]/unit
	
	mystring = sprintf("@%s: %s scaffolds, %.2f#{opts[:unit]} (%.2f%% of total #{opts[:unit]}), %.2f distance",key,val,n_bp_most_spec_pred_unit,percent_bp, average_distance)
	puts mystring
	
	mystring = sprintf("@%s scaffold-contig consistency: %.2f%% (%d contigs from %d), %.2f%% (%.2f#{opts[:unit]} from %.2f#{opts[:unit]})", key,percent_consistent_pred, n_in_path[key]+n_cont_most_spec_pred[key], n_cont[key], percent_contig_bp, n_bp_in_path_unit, n_bp_most_spec_pred_unit)
	puts mystring

	n_clades += 1
	avg_contig_consistency += percent_consistent_pred
	avg_bp_consistency += percent_contig_bp
	
	#old puts stuff
	#puts "@#{key}: #{val} scaffolds, #{n_bp_most_spec_pred_unit}#{opts[:unit]} (#{percent_bp}% of total bp)"
	#puts "@Contigs of these are #{key}: #{n_cont_most_spec_pred[key]} (#{percent_contigs_most_spec}%)"
	#puts "@Scaffold-Contig consistency: #{percent_consistent_pred}% (#{n_in_path[key]+n_cont_most_spec_pred[key]} contigs from #{n_cont[key]}), #{percent_contig_bp}% (#{n_bp_in_path[key]}bp from #{n_bp_most_spec_pred[key]}bp)"
end

mystring = sprintf("@Summary: Number of clades %d, Average contig consistency %.3f%%, Average bp consistency %.3f%%", n_clades, avg_contig_consistency/n_clades, avg_bp_consistency/n_clades)
puts mystring

puts "@=============== End Of clade statistics ===============@"
puts "#=============== End Of File ===============#"


graph_legend.sort{|a,b| a[1]<=>b[1]}.each {|p,c| graph_legend_file.puts "#{p}\t#{c}"}

#close graph files
if opts[:grapf]
	graph_mat_file.close
	graph_col_file.close
	graph_rname_file.close
	graph_legend_file.close
end

