#!/usr/bin/ruby
#Author: Kaustubh Patil (patil@mpi-inf.mpg.de)
#Copyright (c) Kaustubh Patil 2010
#Disclaimer: This program is provided WITHOUT ANY WARRANTY.
#Author accepts no responsibility to anyone for the consequences of using it or for whether it serves any particular purpose or works at all. 
#Permission is hereby granted to copy and redistribute this software for non-commercial use under the restriction that the program is not 
#modified in any way and all disclaimers and acknowledgments remain intact.
#About:
#prediction with already constructed models
#arguments: fasta file to predict labels for and configuration file (same one used for training) [model directory]

mydir = File.dirname(__FILE__)

require "#{mydir}/lib/sqlite_utils.rb"
require "#{mydir}/lib/utils.rb"

$LOAD_PATH.unshift("#{mydir}/lib/rubygems/svg_graph_0.6.1")
require 'SVG/Graph/Pie'

fasta_file=ARGV[0]
config_file=ARGV[1]
model_dir_user = ARGV[2]

raise "\nUsage: predict.rb fasta_file config_file\n" if fasta_file.nil? || config_file.nil?

#read configuraion file
print "Predicting .../#{File.basename(fasta_file)} according to #{config_file} (#{Time.now})"
print "\nReading configuration file and checking if things are in place..."

#read configuration file and get all information
project_dir=""
ncbi_tax_dir=""
frag_len = []
kmer=[]
taxonomy_ranks = []
n_classifiers=0
pie_charts=FALSE

fp = File.open(config_file, "r")
while(line = fp.gets)
        next if line.nil? || line[0].chr=="#"
        line.chomp!
        line = line.split(":")
        if line[0]=="PROJECT_DIR"
                project_dir = line[1]
	elsif line[0]=="NCBI_TAX_DIR"
		ncbi_tax_dir=line[1]	
        elsif line[0]=="FRAGMENT_LEN"
                frag_len = line[1].split(",")
		frag_len.each_index {|i| frag_len[i]=frag_len[i].to_i}
        elsif line[0]=="KMER"
                kmer = line[1].split("-").flatten
        elsif line[0]=="N_CLASSIFIERS"
                n_classifiers = line[1].to_i
	elsif line[0]=="REV_COMPLEMENT"
		rev_complement = line[1]
	elsif line[0]=="RM_REV_COMPLEMENT"
		rm_rev_complement = line[1]
	elsif line[0]=="KMER_NORMALIZATION"
		normalization = line[1]
	elsif line[0]=="PIE_CHARTS"
		pie_charts=TRUE if line[1]=~/TRUE/
        elsif line[0]=="TAXONOMY_RANKS"
		taxonomy_ranks = line[1].split(",")
        end
end

model_dir = "#{project_dir}/models"
model_dir = model_dir_user if !model_dir_user.nil?

#check if configuration file was ok
raise "Invalid PROJECT_DIR entry in the configuration file" if project_dir.empty? || !File.directory?(project_dir)
raise "Invalid FRAGMENT_LEN entry in the configuration file" if frag_len.empty?
raise "Invalid KMER entry in the configuration file" if kmer.empty?
raise "Invalid TAXONOMY_RANKS entry in the configuration file" if taxonomy_ranks.empty?
raise "Invalid N_CLASSIFIERS entry in the configuraion file" if n_classifiers==0
raise "Invalid REV_COMPLEMENT entry in the configuration file" if rev_complement.nil? || rev_complement.empty?
raise "Invalid KMER_NORMALIZTAION entry in the configuration file" if normalization.nil? || normalization.empty?
raise "Invalid NCBI_TAX_DIR entry in the configuration file" if ncbi_tax_dir.empty? || !File.directory?(ncbi_tax_dir)

#check if the directories exist
raise "\n\nThe model directory either doesn't exist or it's not a directory.\n\n" if !File.directory?(model_dir)
raise "\n\nThe project directory either doesn't exist or it's not a directory.\n\n" if !File.directory?(project_dir) && model_dir_user.nil?

#check if all the model files exist in the project directory
model_files = Dir.new(model_dir).entries
models = Hash.new(nil)
frag_len.each do |fl|
	model_files.each do |mf|
		next if mf == ".." || mf == "."
		mf_fl = mf.split(/[^a-zA-Z0-9]/).flatten.first.to_i
		next if mf_fl.nil? || mf_fl==0
		if fl == mf_fl
			raise "Multiple models for fragment length #{fl}: #{mf} and #{models[fl]}" if !models[fl].nil?
			models[fl] = mf
		end
	end
end

#check if all fragment lengths have a model file
frag_len.each {|fl| raise "No model file for fragment length #{fl}" if models[fl].nil? }

#check if the sqlite database exists
raise "\nThe SQLite database doesn't exist.\n" if !File.exists?("#{ncbi_tax_dir}/ncbitax_sqlite.db")
sqlite_taxonomy = SqliteUtils::Taxonomy.new("#{ncbi_tax_dir}/ncbitax_sqlite.db")

print "done"

if n_classifiers>frag_len.length
	print "\nMore classifiers requested than available...changing to #{frag_len.length}"
	n_classifiers=frag_len.length
end

kmer = kmer.collect {|kk| kk.to_i}

#create a non-ACTG filtered fasta file
fasta_file_new = fasta_file + ".nox.fna"
raise "The filtered fasta file exists please remove or move it!\n#{fasta_file_new}" if File.exist?(fasta_file_new) 
fp = File.open(fasta_file_new,"w")
fasta_seq = Bio::FastaFormat.open(fasta_file)
fasta_seq.each do |seq|
	rep = "X" * kmer.max
	pat = "[^A,C,G,T]{#{kmer.max},}"
        pat = Regexp.new(pat, Regexp::IGNORECASE)
        seq.data = seq.data.gsub(pat,rep)
	dna = Bio::Sequence::NA.new(seq.data)
	fp.puts dna.to_fasta(seq.definition, 1000)
end
fp.close
fasta_file = fasta_file_new

print "\nGenerating k-mer features (#{kmer.join("-")})..."

#now generate kmers
#
fasta2kmers_command = "#{mydir}/../bin/fasta2kmers -s 1 -o 1 -h 1 -l 0 -C #{rm_rev_complement} -t #{normalization} -r #{rev_complement}"
fasta2kmers2_command = "#{mydir}/../bin/fasta2kmers2 -a w -s 1 -l 2 -o 1 -b 1 -R #{rm_rev_complement} -h 1 -n #{normalization} -r #{rev_complement}"
use_fasta2kmers2 = TRUE
use_fasta2kmers2 = FALSE if normalization.to_i>1
#identify what kmers
if kmer.length==1
	k = kmer[0]
	fasta2kmers_command = fasta2kmers_command + " -k #{k}"
	fasta2kmers2_command = fasta2kmers2_command + " -k #{k} -j #{k}"
else
	k = kmer.max
	j = kmer.min
	fasta2kmers_command = fasta2kmers_command + " -k #{k} -j #{j}"
	fasta2kmers2_command = fasta2kmers2_command + " -k #{k} -j #{j}"
end
basename=File.basename(fasta_file)
fasta2kmers_command_final = fasta2kmers_command + "  #{fasta_file} | sed 's/^/1 /' > #{fasta_file}.sl"
fasta2kmers2_command_final = fasta2kmers2_command + " -i #{fasta_file} -f #{fasta_file}.sl"
fasta2kmers_command_use = fasta2kmers_command_final if !use_fasta2kmers2
fasta2kmers_command_use = fasta2kmers2_command_final if use_fasta2kmers2
success=system(fasta2kmers_command_use)
raise "\nError in creating k-mer features.\n" if !success

print "done"
print "\nCreating files for different ensembles..."

#go through the fasta file and separate the entries to be classified with different classifiers
seq_classifiers = Hash.new{|hash, key| hash[key] = Array.new;}
fasta_seq = Bio::FastaFormat.open(fasta_file)
n_seq=0
fasta_seq.each do |seq|
	n_seq+=1
	len = seq.length
	diff_classifier=[]
	#get the closet classifier
	frag_len.each {|fl| diff_classifier << (fl-len).abs }
	#find the closest classifier
	minim = diff_classifier.min
	ind = diff_classifier.index(minim)
	seq_classifiers[ind] << n_seq
end

#create different files for different classifiers
#use the already created kmer file for this
#this is slightly inefficient as we have to go through all the file for every classifier
seq_classifiers.each do |ind,lines|
	total_lines = lines.length
	len=frag_len[ind]
	fp = File.open("#{fasta_file}.#{len}.sl","w")
	fp_kmer = File.open("#{fasta_file}.sl","r")
	n_line=0
	printed_lines = 0
	while(line=fp_kmer.gets)
		next if line.nil?
		n_line+=1
		next if !lines.member?(n_line)
		fp.print(line)
		printed_lines += 1
		break if printed_lines >= total_lines
	end
	fp.close
end

print "done"

print "\nPredicting..."
classifier_command = "#{mydir}/../bin/svm_phylo_classify -v 0 "
ensemble_command = "#{mydir}/../bin/svm_phylo_classify_ensemble -v 0 -e 1 "
#now predict with the models
out_files = []
seq_classifiers.each do |ind,lines|
	next if lines.length==0
	len=frag_len[ind]
	test_file = "#{fasta_file}.#{len}.sl"
	out_file = "#{fasta_file}.#{len}.out"
	out_files << out_file	
	#get what classifiers to use
	classifiers_to_use = []
	len_ind = frag_len.index(len)
	(0...n_classifiers).each {|i| classifiers_to_use << frag_len[i+len_ind] if !frag_len[i+len_ind].nil?}
	print "\n\tfragments close to length #{len} with classifiers #{classifiers_to_use.join(",")}"
	n_classifiers_to_use = classifiers_to_use.length
	if n_classifiers_to_use==1
		extra_command = "#{test_file} #{model_dir}/#{models[len]} #{out_file}"
		final_command=classifier_command+extra_command
	else
		extra_command = "-m #{n_classifiers_to_use} #{test_file} "
		classifiers_to_use.each do |i|
			extra_command += "#{model_dir}/#{models[i]} "
		end
		extra_command += "#{out_file}"
		final_command=ensemble_command+extra_command
	end
	success=system(final_command)
	raise "\nError in classification: couldn't run the system command.\n" if !success
end

#join all the outputs
out_file = "#{fasta_file}.out"
success=system("cat #{out_files.join(" ")} > #{out_file}")
raise "\nError in concatenating the outputs.\n" if !success

print "\ndone\n"

#do post processing to create PP like output and pie charts
print "Post-processing... "

fp = File.open(out_file,"r")
paths = Hash.new{|hash, key| hash[key] = Array.new;}
preds = []
heads = []
while(line=fp.gets)
	next if line.nil?
	line = line.chomp.split(/\t/)
	heads << line[0]
	pred = line[4]
	preds << pred
	if paths[pred].empty?
		taxonomy_ranks.each do |rank| 
			paths[pred] << sqlite_taxonomy.parent_atrank(pred,rank)
			puts "#{pred} parent at #{rank}" if sqlite_taxonomy.parent_atrank(pred,rank)=="cellular organisms"
		end
	end
end
fp.close

#create a PP style output file
out_file_pp =  "#{fasta_file}.PP.out"

#write the file
#after converting taxonomy ids to names
fp = File.open(out_file_pp,"w")
fp.puts("#Output in PhyloPythia format (#{Time.now})")
fp.puts("#Taxonomic Hierarchy: NCBI Taxonomy")
fp.puts("#Prediction method: PhyloPythiaS")
fp.puts("#Contact: PhyloPythiaS mailing list (phylopythias@uni-duesseldorf.de)")
fp.puts("\n")

fp.puts("#ID\t#{taxonomy_ranks.reverse.join("\t")}")
paths_name = Hash.new{|hash, key| hash[key] = Array.new;}
preds.each_index do |i|
	pred = preds[i]
	head = heads[i]
	path_name = paths_name[pred]
	if path_name.empty?
		path = paths[pred]
		path_name = []
		path.each do |p|
			name = sqlite_taxonomy.scientific_name(p)
			name = "\t" if name.nil?
			path_name << name
		end
		paths_name[pred] = path_name.clone
	end
	path_for_output = []
	path_for_output << path_name.reverse
	path_for_output << head
	path_for_output.reverse!
	fp.puts(path_for_output.join("\t"))
end
fp.close

#create the pie charts
if pie_charts
	taxonomy_ranks.each_index do |i|
		rank = taxonomy_ranks[i]
		next if rank=="root"
		names = []
		preds.each do |pred| 
			name = paths_name[pred][i]
			name = "Unassigned" if name.nil? || name.empty? || name=="\t"
			names << name
		end
		names=names.sort
		counts = names.count
		#create graph
		graph = SVG::Graph::Pie.new({
			:height => 500,
			:width  => 500,
			
			:show_shadow => true,
			:shadow_offset => 3,
			:show_key_percent => true,
			:expand_greatest => true,
			:expand_gap => 3,

			:fields => counts.keys
			})
		
		graph.add_data({
				:data => counts.values,
				:title => "Taxonomic rank: #{rank}"
				})

		pie_chart_file = "#{fasta_file}.pie_#{rank}.svg"
		fp = File.open(pie_chart_file,"w")
		fp.puts(graph.burn())
		fp.close

	end
end

sqlite_taxonomy.close

print "done\n"
