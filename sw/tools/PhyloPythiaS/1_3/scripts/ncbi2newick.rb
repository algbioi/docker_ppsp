#!/usr/bin/ruby1.8
#Author: Kaustubh Patil (patil@mpi-inf.mpg.de)
#Copyright (c) Kaustubh Patil 2010
#Disclaimer: This program is provided WITHOUT ANY WARRANTY.
#Author accepts no responsibility to anyone for the consequences of using it or for whether it serves any particular purpose or works at all. 
#Permission is hereby granted to copy and redistribute this software for non-commercial use under the restriction that the program is not 
#modified in any way and all disclaimers and acknowledgments remain intact.
#About:
#outputs a tree in newick format using ncbi taxonomy
#input is the node names either in a file with one node per line

#if you want to pass except_file with default taxonomic_ranks then pass "default" as taxonomic_ranks
#also make sure that taxonomic ranks include root as the last element
Usage = "\nUsage: ncbi2newick.rb node_file sqlite_db_file [taxonomic_ranks] [except_file]\nPass 'default' as taxonomic_ranksif you want to use default ranks with except_file\n"

require File.dirname(__FILE__)+"/lib/sqlite_utils.rb"

node_file = ARGV[0]
db_file = ARGV[1]

raise Usage if node_file.nil? || db_file.nil?

# create an object for the Sqlite Taxonomy class
t = SqliteUtils::Taxonomy.new(db_file)

#possibility to add nodes in the tree irrespective of their rank
addnodes = FALSE
allranks = FALSE

#taxonomic ranks
#make sure to have "root" as an allowed level if you want a rooted tree
levels = ["genus","family","order","class","phylum","superkingdom","root"]
levels = ARGV[2].split(",") if !ARGV[2].nil? || ARGV[2].downcase!="default"

raise "\nRoot must be the last element of taxonomic ranks\n" if levels.last!="root"

#allowed leaf nodes are read from the pipe
leaves = Array.new
fp = File.new(node_file, "r")
while(leaf = fp.gets)
	leaves << leaf.chomp if !leaf.nil?
end
fp.close()

#read in the nodes that must be there
file = ARGV[3]
except = []
if !file.nil?
	fp = File.new(file, "r")
	while(leaf = fp.gets)
		except << leaf.chomp if !leaf.nil?
	end
	fp.close()
end

raise "\nUsage: ncbi2newick.rb node_file sqlite_db_file [except_file]\n" if leaves.empty?

#only process unique nodes
leaves = leaves.uniq

# now construct the tree

phylotree = t.tree_from_nodes(leaves,levels, allranks, addnodes)

#raise "\nUnrooted tree" if phylotree.root.nil?
print phylotree.newick(options={:indent=>false})

if !except.empty?
	phylotree = t.remove_nonsense_nodes_but_not_root(phylotree,except)
	puts "Pruned tree: #{phylotree.number_of_nodes}"
	print phylotree.newick(options={:indent=>false})
end

t.close

