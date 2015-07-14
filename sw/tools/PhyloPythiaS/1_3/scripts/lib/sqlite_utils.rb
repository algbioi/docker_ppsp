#Author: Kaustubh Patil (patil AT mpi-inf DOT mpg DOT de)
#Copyright (c) Kaustubh Patil 2010
#Disclaimer: This program is provided WITHOUT ANY WARRANTY.
#Author accepts no responsibility to anyone for the consequences of using it or for whether
#it serves any particular purpose or works at all. Permission is hereby granted to copy, modify and
#redistribute this software for non-commercial use under the restriction that all disclaimers, copyrights
#and acknowledgments remain intact.

#methods to query the custom built sqlite3 database
#containing ncbi taxonomy
#also parent_taxon_id is not the parent's ncbi_taxon_id but the taxon_id
#thats why methods are different in biosql and sqlite
#in sqlite parent_taxon_id is its ncbi_taxon_id

mydir=File.dirname(__FILE__)

#create paths
$LOAD_PATH.unshift("#{mydir}/rubygems/sqlite3-ruby-1.2.5/lib")
$LOAD_PATH.unshift("#{mydir}/rubygems/bioruby/lib")

#some requires
require "rubygems"
require "sqlite3"
require "bio"
require "bio/tree"

module SqliteUtils

class Taxonomy 

	attr_reader :db_file, :db

	# class initializer
	def initialize(db_file)
		raise "Could not open the database with provided file!" if !File.exists?(db_file)
		@db_file = db_file
		@db=SQLite3::Database.new(@db_file)
	end

	# close the sqlite database
	def close
		@db.close
	end

	#get if its a valid ncbi_taxon_id
	def ncbi_taxon_id_exists(tax_id)
		res = @db.execute("select taxon_id from taxon where ncbi_taxon_id=?",tax_id)
		return res[0][0] if res.length==1
		return nil if res.length==0
		raise "\nMultiple ncbi_taxon_id: #{tax_id}\n" if res.length>1
	end

	def scientific_name(tax_id)
	
		return nil if tax_id.nil? || tax_id.to_s.empty?  #changed I

		#first get the taxon_id from given ncbi_taxon_id
		taxon_id = get_taxon_id(tax_id)

		return nil if taxon_id.nil?
		
		#raise "Couldnot find taxon_id for #{tax_id}" if taxon_id.nil?

		res = @db.execute("select name from taxon_name where taxon_id=? and name_class='scientific name'",taxon_id)
		return nil if res.length==0
		raise "More than one parent..." if res.length>1
		return res[0][0]
	end

	def all_name(tax_id)
	
		raise "\nInvalid argument in all_name\n" if tax_id.nil? || tax_id.empty?

		#first get the taxon_id from given ncbi_taxon_id
		taxon_id = get_taxon_id(tax_id)

		raise "Couldnot find taxon_id for #{tax_id}" if taxon_id.nil?

		res = @db.execute("select name from taxon_name where taxon_id=?",taxon_id)
		return nil if res.length==0
		res = res.flatten
		return res
	end

	def rank(tax_id)
		# return "root" if tax_id is 1, its useful to
		# differentiate with other "no rank" nodes
		return nil if tax_id.nil? || tax_id.empty?
		return "root" if tax_id.to_i==1
		res = @db.execute("select node_rank from taxon where ncbi_taxon_id=?",tax_id)
		raise "More than one ranks..." if res.length > 1
		return nil if res.length == 0
		return res[0][0]
	end

	def tax_id(name)
		#returns taxonomic id given any name	
		return nil if name.empty? || name.nil?
		#escape the name
		name = SQLite3::Database.quote(name)
		res = @db.execute("select taxon_id from taxon_name where name=?",name)
		#raise "Could not find tax_id for: #{name}" if res.length == 0
		return nil if res.length == 0
                #raise "Found multiple tax_id for: #{name}" if res.length > 1
		taxon_id = res[0][0]

		#get corresponding ncbi_taxon_id
		t = get_ncbi_taxon_id(taxon_id)
		return t
		
	end

	def tax_id_rlike(name)
		return nil if name.empty? || name.nil?
		#escape the name
		res = @db.execute("select taxon_id from taxon_name where name r like ?",name)
		
		return nil if res.length == 0
		return nil if res.length > 1
		
		taxon_id = res[0][0]
		
		#get corresponding ncbi_taxon_id
		t = get_ncbi_taxon_id(taxon_id)
		return t
	end
	
	#get ncbi_taxon_id given taxon_id
	def get_taxon_id(ncbi_taxon_id)
		res = @db.execute("select taxon_id from taxon where ncbi_taxon_id=?",ncbi_taxon_id)
		return nil if res.length == 0
                raise "Found multiple taxon_id for ncbi_taxon_id: #{ncbi_taxon_id}" if res.length > 1
		return res[0][0]
	end

        #get ncbi_taxon_id given taxon_id
	def get_ncbi_taxon_id(taxon_id)
		res = @db.execute("select ncbi_taxon_id from taxon where taxon_id=?",taxon_id)
		raise "Could not find ncbi_taxon_id for taxon_id: #{taxon_id}" if res.length == 0
		raise "Found multiple ncbi_taxon_id for taxon_id: #{taxon_id}" if res.length > 1
		return res[0][0]
	end

	def parent(tax_id)
		#get immediate parent
		#first get the taxon_id
		taxon_id = get_taxon_id(tax_id)	
		res = @db.execute("select parent_taxon_id from taxon where taxon_id=?",taxon_id)
                raise "Could not find parent for: #{tax_id}" if res.length == 0
		#return nil if res.length == 0
		raise "Found multiple parents for: #{tax_id}" if res.length > 1
		parent_taxon_id = res[0][0]

		#no need to convert to ncbi_taxon_id as it is ncbi_taxon_id
		return parent_taxon_id
	end

	def children(tax_id)
		#get the taxon_id
		#no taxon_id necessary for sqlite
		#taxon_id = get_taxon_id(tax_id)

		res = @db.execute("select taxon_id from taxon where parent_taxon_id=?",tax_id)
		children = []
		res = res.flatten
		res.each {|child| children << get_ncbi_taxon_id(child)}
		return children
	end

	def parent_atrank(tax_id,rank)
		#get parent node of a given ncbi_taxon_id at given rank
		#will return first occurence of a parent at gievn rank (e.g. if rank is "no rank")
		#this method can be made more efficient and understandable

		return nil if tax_id.nil? || tax_id.to_i <= 0
		return "1" if rank == "root"

		#get taxon_id
		taxon_id = get_taxon_id(tax_id)
		return nil if taxon_id.nil?

		my_rank = ""
		parent_taxon_id = taxon_id # set self as parent
		niter=0
		while TRUE
			taxon_id = taxon_id.to_s
			res = @db.execute("select parent_taxon_id,node_rank from taxon where taxon_id=?",taxon_id)
			return nil if res.length == 0 # no parent at this level
			#raise "No parent found: #{tax_id}" if res.length == 0
			raise "Multiple parent found for: #{tax_id}" if res.length > 1
			parent_taxon_id = res[0][0]
			my_rank = res[0][1]
			if parent_taxon_id.to_s == "1"
				#reached highest level
				return nil
			elsif my_rank==rank
				#return ncbi_taxon_id of current taxon_id
				return get_ncbi_taxon_id(taxon_id)
			else #continue to search upwards
				#convert to taxon_id for sqlite
				taxon_id = get_taxon_id(parent_taxon_id)
			end
			niter+=1
			raise "parent_atrank: unusual looping, more than 100 iterations!" if niter>100
		end
		
		#get ncbi_taxon_id
		return get_ncbi_taxon_id(taxon_id)

	end

	def path_toroot(tax_id)
		#returns path to the root including the node itself
		tax_id=tax_id.to_s
		path = []
		path << tax_id
		while tax_id.to_i!=1	#1 is the root node
			p = parent(tax_id)
			next if p.nil?
			tax_id = p
			path << tax_id
		end
		return path
	end

	#following three defs are here just for convenience
	def tree_first_child(tree, node)
		children = tree.children(node)
		return nil if children.empty?
		return children.first
	end

        def tree_last_child(tree, node)
		children = tree.children(node)
		return nil if children.empty?
		return children.last
	end

	def tree_path_to_leaf(tree,node)
		path = [node]
		child = tree.children(node).first
		while !child.nil?
			path << child
			child = tree.children(child).first
		end
		return path
	end

	# get a tree for given leaves
	# ranks are allowed ranks in the taxonomy (should be bottom up, 
	# e.g. ["genus","family","order","class","phylum","superkingdom"])
	# allranks=T gives equal branch lengths
	# add_nodes adds given nodes to the path irrespective of their rank
	def tree_from_nodes(nodes, ranks=nil, allranks=FALSE, add_nodes=FALSE)

		phylotree = Bio::Tree.new
		#get all ranks if not provided
		if ranks.nil? || ranks.empty?
			#ordering of ranks cannot be derived from the database
			#so set it to a default
			ranks = ["genus","family","order","class","phylum","superkingdom","root"]
		end

		# for each leaf add corresponding parents
		rank_cache = Hash.new(nil)
		path_to_root = []
		dummy_node = 0
		nodes.each do |leaf|
			#some of the nodes could have been already added through path of some other node
			#node=phylotree.get_node_by_name(leaf.to_s)
			#node = Bio::Tree::Node.new(leaf.to_s) if node.nil?
			#phylotree.add_node(node) #if the node exists then this does nothing
			
			path_to_root = []
			path_to_root << leaf if add_nodes
			ranks.each do |r|
				parent = parent_atrank(leaf, r)
				if parent.nil? && allranks
					dummy_node = dummy_node-1
					parent = dummy_node.to_s
				end
				path_to_root << parent if !parent.nil?
			end

			#take care of nil paths
			#path has length of 1 if the node has no parents or is not present
			next if path_to_root.length <= 1
			child =  path_to_root[0]
			path_to_root.each do |p|
				node=phylotree.get_node_by_name(p.to_s)
				node = Bio::Tree::Node.new(p.to_s) if node.nil?
				phylotree.add_node(node) #if the node exists then this does nothing			
				child_node = phylotree.get_node_by_name(child.to_s)
				phylotree.add_edge(node,child_node) if node != child_node
				#print "\nAdding link between: " + p.to_s + " and " + child.to_s
				child = p

				#raise "Nil child" if child.nil?
			end

			#print phylotree.newick(options={:indent=>false})
	
		end
		#set the tree root
		node=phylotree.get_node_by_name(path_to_root[-1].to_s)
		phylotree.root=node
		phylotree
	end

	#remove nonsense nodes but not the root node as the original method does
	#donot remove nodes with name in except
	def remove_nonsense_nodes_but_not_root(tree,except=[])
		nodes = tree.nodes
		nonsense_nodes = []
		nodes.each {|node| nonsense_nodes << node if tree.out_degree(node)==2 && tree.root!=node}

		nonsense_nodes.each do |node|
			#first check if the node exists
			#is this redundant?
			name = tree.get_node_name(node)
			node = tree.get_node_by_name(name)
			next if node.nil?
			
			child = tree.children(node)[0]
			parent = tree.parent(node)

			if !except.member?(name)
				tree.clear_node(node)
				tree.add_edge(child, parent)
			end

		end

		tree

	end

end # class Taxonomy ends

end # end of SqliteUtils module

