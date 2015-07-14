#!/usr/bin/ruby
# given a directory it reads each file (assuming fasta format) and outputs 
# fragments of given length
#Usage
# fasta2fragments.rb frag_length fasta_file [max_frag]

#use own bioruby
mydir=File.dirname(__FILE__)

$LOAD_PATH.unshift("#{mydir}/lib/rubygems/bioruby/lib")
require "bio"

textwidth = 100 # to write fasta file !!CURRENTLY NOT USED!!
step_size = 0

#read-input parameters
len = ARGV[0].to_i
ifile = ARGV[1]
max = ARGV[2].to_i

raise "\nFiven argument (#{ifile}) is not a file\n" if !FileTest.file?(ifile)

fasta = Bio::FastaFormat.open(ifile)
raise "Invalid fasta file: #{ff}" if fasta.nil?


#get number of sequences and max_frag for each sequence
nseq = 0
fasta.each {|seq| nseq+=1}

max = (max/nseq).ceil if !max.nil?

#rewind fasta
fasta.rewind

seq_concat = ""
head = ""
n_seq=0
n_frag=0
fasta.each do |seq|
	raise "Invalid sequence (nil): #{ifpath}" if seq.nil?
	next if (seq.data).nil?

	n_seq+=1
		
	seq_concat = seq.data
	head = seq.definition

        #remove any newlines
	#NOTE seq.data[0] is '\n' character
	seq_concat = seq_concat.gsub("\n","")

	if seq_concat.length < len
		STDERR.puts "#{head} is too short to fragement"
		next
	end

	no_frag = (seq_concat.length/len).floor.to_i
	max_frag = max
	max_frag = no_frag if max_frag==0 #take all if none asked


	ind = (1..no_frag).to_a
	if no_frag > max_frag
		ind = ind.sort_by{rand}
		ind = ind[0...max_frag]
	end

	#print "\nIndex: " + ind.join(" ") + "\n"
	dna = Bio::Sequence::NA.new(seq_concat)
	
	if step_size==0
        	step = len
	else
        	step = step_size
	end

	i = 0
	start = 0

	remainder = dna.window_search(len, step) do |s|
		start = len*i			
		i += 1
		next if max_frag<no_frag && !ind.member?(i)
		puts s.to_fasta("#{head}|seq#{n_seq}|frag#{i}_at#{start}")
	end

	n_frag += no_frag
	break if n_frag >= max_frag
end

