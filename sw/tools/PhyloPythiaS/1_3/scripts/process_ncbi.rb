#!/usr/bin/ruby -w
#process ncbi download by creating fasta files with appropriate names
#supports gbk and gbff files

mydir = File.dirname(__FILE__)

require mydir+"/lib/utils.rb"

ncbi = NCBI.new
extra = Extra.new

Usage = "USAGE: process_ncbi.rb ncbi_download_dir output_dir [exclude_plasmids/phage:t/f]"

ncbi_dir = ARGV[0]
output_dir = ARGV[1]
exclude_plasmids = FALSE
exclude_plasmids = TRUE if !ARGV[2].nil? && ARGV[2].downcase=="t"

raise Usage if ncbi_dir.nil? || output_dir.nil?

#process the following files
process_extensions = ["gbk","gbff"]

action_file_exist = "A" ##A=Ask, S=Skip, O=Overwrite, SA=Skip All, OA=Overwrite All
valid_answers_file_exist = ["S","O","SA","OA"]
question_file_exist = "The file exists what should we do [S:Skip(default),O:Overwrite,SA:Skip All,OA:Overwrite All]? "

print "\nPreparing ..."

ncbi_dir_contents = Dir.new(ncbi_dir).entries
#collect a list of files to be processed
files_to_process = Hash.new(nil)
ncbi_dir_contents.each do |ndc|
	next if ndc == "." || ndc == ".."
	if File::directory?(ncbi_dir+"/"+ndc) 
		Dir.new(ncbi_dir+"/"+ndc).entries.each do |c|
			next if c == "." || c == ".."
			files_to_process[ncbi_dir+"/"+ndc+"/"+c] = ndc if process_extensions.member? c.split(/\./).last
		end
		next
	end

	files_to_process[ncbi_dir+"/#{ndc}"] = ncbi_dir if process_extensions.member? ndc.split(/\./).last
end

raise "\nNo files to process.\n" if files_to_process.empty?

print "done\nProcessing (this will take a while)..."

processed_taxon = []
processed_dir = []
n_file=0
files_to_process.each do |f,d|
	
	n_file+=1
	#print something every 100 files so that user is informed about progress
	print "*" if n_file%100==0
	
	parsed = ncbi.parse_gbk_gbff(f)
	taxon =  parsed["taxon"][0]
	definition = parsed["definition"][0]
	next if exclude_plasmids && (!definition.scan(/plasmid/i).empty? || !definition.scan(/phage/i).empty?)
	
	file = output_dir+"/"+taxon+".1.fna"
	fo = nil
	if File.exist?(file) && !processed_taxon.member?(taxon)
		if action_file_exist=="SA"
			next
		elsif action_file_exist=="OA"
		#this is meant to take care of multiple gbk files in one ncbi dir
			if !processed_dir.member? d
				fo = File.open( file, "w" )
			else
				fo = File.open( file, "a" )
			end
		else
			action_file_exist = extra.get_user_input(question_file_exist+taxon+".1.fna").chomp
			action_file_exist = "S" if !valid_answers_file_exist.member?(action_file_exist)
			if action_file_exist=="S" || action_file_exist=="SA"
				next
			else
				fo = File.open( file, "w" )		
			end
			action_file_exist="A" if action_file_exist=="S" || action_file_exist=="O"
		end
	else
		fo = File.open( file, "a" )
	end

	raise "Something went wrong in file opening" if fo.nil?

	#write to file
	accession = parsed["accession"]
	gi = parsed["gi"]
	seq = parsed["seq"]

	gi.each_index do |i|
	        fo.puts ">gi|#{gi[i]}|gb|#{accession[i]}|#{definition}\n#{seq[i]}"
	end
	fo.close

	processed_taxon << taxon
	processed_dir << d
end

print "done\nProcessing of raw NCBI data finished."

