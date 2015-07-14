#Changes in the code due to changes in the biosql_1_0_1 schema
#for new schema taxon_id != ncbi_taxon_id

# utilities
# some requires
$LOAD_PATH.unshift("rubygems/bioruby/lib")
$LOAD_PATH.unshift("rubygems/narray")
$LOAD_PATH.unshift("rubygems/narray/lib")
$LOAD_PATH.unshift("rubygems/ruby-mysql")

require 'mysql'
require 'bio'
require 'bio/tree'

module BioSqlUtils

class Taxonomy 

	attr_accessor :host, :user, :pass, :db

	# class initializer
	#def initialize(host="infno6740",user="patil",pass="z6Ly3rI3",db="biosql_1_0_0")
	def initialize(host="bioinfodb",user="patil",pass="z6Ly3rI3",db="patil_biosql_1_0_1")
		# some mysql settings
		@host = host
		@user = user
		@pass = pass
		@db =db
		begin 
			@connection=Mysql::new(@host, @user, @pass, @db)
			$stdout.flush

		rescue

			puts "Could not connect to the database with default settings!"

			user_in = user_input("Please enter host name [empty for default]: ").chomp
			@host = user_in if !user_in.nil? && !user_in.empty?
			user_in = user_input("Please enter database name [empty for default]: ").chomp
			@db = user_in if !user_in.nil? && !user_in.empty?
			user_in = user_input("Please enter username [empty for default]: ").chomp
			@user = user_in if !user_in.nil? && !user_in.empty?
			user_in = user_input("Please enter password [empty for default]: ").chomp
			@pass = user_in if !user_in.nil? && !user_in.empty?
			@connection=Mysql::new(@host, @user, @pass, @db)		
			
		end
		
	end

	# close the mysql connection
	def close
		@connection.close
	end

	def user_input(string)
        	print string
                answer = STDIN.gets
                answer
	end

	#get if its a valid ncbi_taxon_id
	def ncbi_taxon_id_exists(tax_id)
		res = @connection.query("select taxon_id from taxon where ncbi_taxon_id=#{tax_id}")
		return res.fetch_row[0] if res.num_rows==1
		return nil if res.num_rows==0
		raise "\nMultiple ncbi_taxon_id: #{tax_id}\n" if res.num_rows>1
	end

	def scientific_name(tax_id)
	
		raise "\nInvalid argument in scientific_name\n" if tax_id.nil? || tax_id.empty?

		#first get the taxon_id from given ncbi_taxon_id
		taxon_id = get_taxon_id(tax_id)

		return nil if taxon_id.nil?
		
		#raise "Couldnot find taxon_id for #{tax_id}" if taxon_id.nil?

		res = @connection.query("select name from taxon_name where taxon_id=#{taxon_id} and name_class='scientific name'")
		return nil if res.num_rows==0
		raise "More than one parent..." if res.num_rows>1
		return res.fetch_row[0]
	end

	def all_name(tax_id)
	
		raise "\nInvalid argument in all_name\n" if tax_id.nil? || tax_id.empty?

		#first get the taxon_id from given ncbi_taxon_id
		taxon_id = get_taxon_id(tax_id)

		raise "Couldnot find taxon_id for #{tax_id}" if taxon_id.nil?

		res = @connection.query("select name from taxon_name where taxon_id=#{taxon_id}")
		return nil if res.num_rows==0
		names = []
		res.each {|row| names << row[0] }
		return names
	end

	def rank(tax_id)
		# return "root" if tax_id is 1, its useful to
		# differentiate with other "no rank" nodes
		return nil if tax_id.nil?
		return "root" if tax_id.to_i==1
		res = @connection.query("select node_rank from taxon where ncbi_taxon_id='#{tax_id}'")
		raise "More than one ranks..." if res.num_rows > 1
		return nil if res.num_rows == 0
		return res.fetch_row[0]
	end

	# returns taxonomic id given the any name
	def tax_id(name)
	
		return nil if name.empty? || name.nil?
		#escape the name
		name = @connection.escape_string(name)
		res = @connection.query("select taxon_id from taxon_name where name='#{name}'")
		#raise "Could not find tax_id for: #{name}" if res.num_rows == 0
		return nil if res.num_rows == 0
                #raise "Found multiple tax_id for: #{name}" if res.num_rows > 1
		taxon_id = res.fetch_row[0]

		#get corresponding ncbi_taxon_id
		t = get_ncbi_taxon_id(taxon_id)
		return t
		
	end

	def tax_id_rlike(name)
		return nil if name.empty? || name.nil?
		#escape the name
		name = @connection.escape_string(name)
		res = @connection.query("select taxon_id from taxon_name where name r like '#{name}'")
		
		return nil if res.num_rows == 0
		return nil if res.num_rows > 1
		
		taxon_id = res.fetch_row[0]
		
		#get corresponding ncbi_taxon_id
		t = get_ncbi_taxon_id(taxon_id)
		return t
	end
	
	#get ncbi_taxon_id given taxon_id
	def get_taxon_id(ncbi_taxon_id)
		res = @connection.query("select taxon_id from taxon where ncbi_taxon_id=#{ncbi_taxon_id}")
                #raise "Could not find taxon_id for ncbi_taxon_id: #{ncbi_taxon_id}" if res.num_rows == 0
		return nil if res.num_rows == 0
                raise "Found multiple taxon_id for ncbi_taxon_id: #{ncbi_taxon_id}" if res.num_rows > 1
		return res.fetch_row[0]
	end

        #get ncbi_taxon_id given taxon_id
	def get_ncbi_taxon_id(taxon_id)
		res = @connection.query("select ncbi_taxon_id from taxon where taxon_id=#{taxon_id}")
		raise "Could not find ncbi_taxon_id for taxon_id: #{taxon_id}" if res.num_rows == 0
		raise "Found multiple ncbi_taxon_id for taxon_id: #{taxon_id}" if res.num_rows > 1
		return res.fetch_row[0]
	end


	# get immediate parent
	def parent(tax_id)
		#first get the taxon_id
		taxon_id = get_taxon_id(tax_id)	

		res = @connection.query("select parent_taxon_id from taxon where taxon_id=#{taxon_id}")
                raise "Could not find parent for: #{tax_id}" if res.num_rows == 0
		#return nil if res.num_rows == 0
		raise "Found multiple parents for: #{tax_id}" if res.num_rows > 1
		parent_taxon_id = res.fetch_row[0]

		#now get the ncbi_taxon_id
		ncbi_taxon_id = get_ncbi_taxon_id(parent_taxon_id)

		return ncbi_taxon_id
	end

	def children(tax_id)
		#get the taxon_id
		taxon_id = get_taxon_id(tax_id)

		res = @connection.query("select taxon_id from taxon where parent_taxon_id=#{taxon_id}")
		children = Array.new
		res.each {|child| children << get_ncbi_taxon_id(child)}
		children
	end

	# this method can be made more efficient and understandable
	def parent_atrank(tax_id,rank)

		return nil if tax_id.to_i < 0
		return "1" if rank == "root"

		#get taxon_id
		taxon_id = get_taxon_id(tax_id)

		return nil if taxon_id.nil?

		my_rank = ""
		parent_taxon_id = taxon_id # set self as parent
		while TRUE
			# rank is own rank
			res = @connection.query("select parent_taxon_id,node_rank from taxon where taxon_id=#{taxon_id}")
			return nil if res.num_rows == 0 # no parent at this level
			#raise "No parent found: #{tax_id}" if res.num_rows == 0
			raise "Multiple parent found: #{tax_id}" if res.num_rows > 1
			row = res.fetch_row
			parent_taxon_id = row[0]
			my_rank = row[1]
			if parent_taxon_id == "1"
				#reached highest level
				return nil
			elsif "#{my_rank}" == "#{rank}"
				return get_ncbi_taxon_id(taxon_id)
			else #continue to search upwards
				taxon_id = parent_taxon_id
			end
		end
		
		#get ncbi_taxon_id
		return get_ncbi_taxon_id(taxon_id)

	end

	def path_toroot(tax_id)
		path = Array.new
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
	# ranks are allowed ranks in the taxonomy
	# allranks=T gives equal branch lengths
	# add_nodes adds given nodes to tha path irrespective of their rank
	def tree_from_nodes(nodes, ranks, allranks, add_nodes)

		phylotree = Bio::Tree.new
		#get all ranks if not provided
		if ranks.empty?
			res = @connection.query("SELECT DISTINCT node_rank FROM taxon")
			ranks = []
			res.each do |row|
				ranks << row[0]
			end
		end

		# for each leaf add corresponding parents
		rank_cache = Hash.new(nil)
		path_to_root = []
		dummy_node = 0
		nodes.each do |leaf|
			#some of the nodes could have been already added through
			#path of some node
			#node=phylotree.get_node_by_name(leaf.to_s)
			#node = Bio::Tree::Node.new(leaf.to_s) if node.nil?
			#phylotree.add_node(node) #if the node exists then this does nothing

			path_to_root = []
			path_to_root << leaf if add_nodes
			ranks.each do |r|
				parent = parent_atrank(leaf, r)
				if parent.nil? && allranks
					dummy_node = dummy_node-1
					parent = dummy_node
				end
				path_to_root << parent if !parent.nil?
			end

			# take care of nil paths
			# path has length of 1 if the node has no parents or is not present
			next if path_to_root.length == 1 || path_to_root.length == 0
			child =  path_to_root[0]
			path_to_root.each do |p|

				#if !allranks
				#	rank=rank_cache[p]
				#		
				#	if rank.nil?
				#		rank = rank(p).to_s
				#		rank_cache[p] = rank
				#	end
				#
				#	#print  "\n"+p.to_s+" has rank "+rank
				#	raise "Nil rank for #{p}" if rank.nil?
				#	if !ranks.member?(rank)
				#		#puts "Skipping: #{p} (#{rank})"
				#		next
				#	end
				#end

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
	def remove_nonsense_nodes_but_not_root(tree,except)
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

end # end of BioSqlUtils module

