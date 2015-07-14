#!/usr/bin/ruby -w
#Author: Kaustubh Patil (patil@mpi-inf.mpg.de)
#Copyright (c) Kaustubh Patil 2010
#Disclaimer: This program is provided WITHOUT ANY WARRANTY.
#Author accepts no responsibility to anyone for the consequences of using it or for whether it serves any particular purpose or works at all. 
#Permission is hereby granted to copy and redistribute this software for non-commercial use under the restriction that the program is not 
#modified in any way and all disclaimers and acknowledgments remain intact.
#About:
#generate training data and train models on it
#input file is the configuration file

mydir = File.dirname(__FILE__)

require "fileutils" #for recursively removing directory

require mydir+"/lib/utils.rb"
require mydir+"/lib/sqlite_utils.rb"

usage = "USAGE: train.rb configuration_file"

config_file = ARGV[0]
raise usage if config_file.nil?

extra = Extra.new

#######Advanced training options##############

#learning
loss_function = 1 #0:0/1 loss, 1:path loss
z_standardization = 1 #0:no standardization, 1:z standardization
misc_nodes = 1 #0:no misc nodes, 1:add misc nodes

#############################################

puts "Processing according to #{config_file} (#{Time.now})"

#read configuration file and get all information
ncbi_dir  = ""
project_dir=""
tree_file = ""
sample_specific_dir=""
frag_len = []
taxonomy_ranks = []
n_examples = ""
sample_specific_step = []
kmer=[]
c_grid=[]
only_models=FALSE
clean_up_train=FALSE
ncbi_tax_dir=""
genomes_exclude_file = nil
extensions = []
parallel_models = FALSE
balance_classes = FALSE
n_min_genomes_generic=nil

fp = File.open(config_file, "r")

while(line = fp.gets)
  next if line.nil? || line[0].chr=="#"
  line = line.chomp
  line = line.split(":")
  next if line[0].nil? || line[1].nil?
  if line[0]=="NCBI_PROCESSED_DIR"
    ncbi_dir = line[1]
  elsif line[0]=="PROJECT_DIR"
    project_dir = line[1]
  elsif line[0]=="TREE_FILE"
    tree_file = line[1]
  elsif line[0]=="SAMPLE_SPECIFIC_DIR"
    sample_specific_dir = line[1]
  elsif line[0]=="FRAGMENT_LEN"
    frag_len = line[1].split(",")
  elsif line[0]=="TAXONOMY_RANKS"
    taxonomy_ranks = line[1].split(",")
  elsif line[0]=="NUMBER_EXAMPLES"
    n_examples = line[1].to_f
  elsif line[0]=="KMER"
    kmer = line[1].split("-").flatten
  elsif line[0]=="C_GRID"
    c_grid = line[1].split(",").flatten
  elsif line[0]=="ONLY_MODELS"
    only_models = TRUE if line[1]=~/TRUE/
  elsif line[0]=="CLEAN_UP_TRAIN"
    clean_up_train = TRUE if line[1]=~/TRUE/
  elsif line[0]=="REV_COMPLEMENT"
    rev_complement = line[1]
  elsif line[0]=="RM_REV_COMPLEMENT"
    rm_rev_complement = line[1]
  elsif line[0]=="KMER_NORMALIZATION"
    normalization = line[1]
  elsif line[0]=="LOSS_ACTION"
    loss_action = line[1]
  elsif line[0]=="NCBI_TAX_DIR"
    ncbi_tax_dir = line[1]
  elsif line[0]=="GENOMES_EXCLUDE"
    genomes_exclude_file = line[1]
  elsif line[0]=="EXTENSIONS"
    extensions = line[1].split(",")
  elsif line[0]=="SAMPLE_SPECIFIC_STEP"
    sample_specific_step = line[1].split(",")
  elsif line[0]=="PARALLEL_MODELS"
    parallel_models = TRUE if line[1].downcase=="true"
  elsif line[0]=="BALANCE_CLASSES"
    balance_classes = TRUE if line[1].downcase=="true"
  elsif line[0]=="N_MIN_GENOMES_GENERIC"
    n_min_genomes_generic = line[1].to_i
  elsif line[0]=="KERNEL"
    kernel = line[1].to_i
  elsif line[0]=="KERNEL_POLYNOMIAL_DEGREE"
    kernel_polynomial_degree = line[1].to_i
  elsif line[0]=="KERNEL_POLYNOMIAL_S"
    kernel_polynomial_s = line[1].to_f
  elsif line[0]=="KERNEL_RBF_GAMMA"
    kernel_rbf_gamma = line[1].to_f
  else
    puts "Unknown configuration file option: #{line[0]}...ignoring"
  end
end #while

#check if configuration file was ok
raise "Invalid NCBI_PROCESSED_DIR entry in the configuration file" if ncbi_dir.nil? || ncbi_dir.empty? || !File.directory?(ncbi_dir)
raise "Invalid PROJECT_DIR entry in the configuration file" if project_dir.nil? || project_dir.empty? || !File.directory?(project_dir)
raise "Invalid TREE_FILE entry in the configuration file" if !tree_file.nil? && !tree_file.empty? && !File.exists?(tree_file)
raise "Invalid FRAGMENT_LEN entry in the configuration file" if frag_len.nil? || frag_len.empty?
raise "Invalid TAXONOMY_RANKS entry in the configuration file" if taxonomy_ranks.empty?
raise "Invalid NUMBER_EXAMPLES entry in the configuration file" if  n_examples.nil? || n_examples==0.0
raise "Invalid KMER entry in the configuration file" if kmer.nil? || kmer.empty?
raise "Invalid C_GRID entry in the configuration file" if c_grid.nil? || c_grid.empty?
raise "Invalid REV_COMPLEMENT entry in the configuration file" if rev_complement.nil? || rev_complement.empty?
raise "Invalid KMER_NORMALIZTAION entry in the configuration file" if normalization.nil? || normalization.empty?
raise "Invalid LOSS_NORMALIZTAION entry in the configuration file" if loss_action.nil? || loss_action.empty?
raise "Invalid NCBI_TAX_DIR entry in the configuration file" if ncbi_tax_dir.nil? || ncbi_tax_dir.empty? || !File.directory?(ncbi_dir)
raise "Invalid SAMPLE_SPECIFIC_DIR entry in the configuration file" if !sample_specific_dir.nil? && !sample_specific_dir.empty? && !File.directory?(sample_specific_dir)
raise "Invalid GENOMES_EXCLUDE entry in the configuration file" if !genomes_exclude_file.nil? && !File.exists?(genomes_exclude_file)

raise "Invalid N_MIN_GENOMES_GENERIC entry, must be a positive integer" if (tree_file.nil? || tree_file.empty?) && (n_min_genomes_generic.nil? || n_min_genomes_generic<=0)

puts "Tree file not provided, will build a generic model with #{n_min_genomes_generic} min clade genomes" if tree_file.nil? || tree_file.empty?

if sample_specific_step.length==1
  sample_specific_step = sample_specific_step*frag_len.length
elsif sample_specific_step.length>1 && sample_specific_step.length!=frag_len.length
  raise "Not same number of SAMPLE_SPECIFIC_STEP for each FRAGMENT_LEN"
end

#get the genomes to exclude
genomes_exclude = []
if !genomes_exclude_file.nil?
  fp = File.open(genomes_exclude_file,"r")
  while(line=fp.gets)
    next if line.nil?
    genomes_exclude << line.chomp	
  end
end

#check if the project directory is empty
entries =  Dir.new(project_dir).entries
if entries.length>2 && !only_models
  print "\nThe project directory #{project_dir} is not empty, this can result in unpredictable behavior."
  print "\nRemove? [Y/N] (default=N)\n"
  answer=STDIN.getc.chr
  answer="n" if answer=="\n" || answer==10 #10 is ascii
  if answer.downcase=="y"
    FileUtils.rm_rf(project_dir)
    Dir.mkdir(project_dir)
  else
    print "\nPlease provide an empty project directory. Quiting...\n"
    exit 1
  end
end

#check if the sqlite database exists otherwise create one
#assuming a constant name for the sqlite database
ncbi_tax_db = ncbi_tax_dir+"/ncbitax_sqlite.db"
if !File.exists?(ncbi_tax_db)
  print "\nThe SQLite database doesn't exist we will create one ..."
  success=system("#{mydir}/ncbitax2sqlite.rb #{ncbi_tax_dir} #{ncbi_tax_db} #{mydir}/sqlite_tax_schema.txt")
  raise "\nError in creating the SQLite database.\n" if !success
end

#create an object for using the taxonomy
sqlite_taxonomy = SqliteUtils::Taxonomy.new(ncbi_tax_db)

if !only_models
  
  #put excluded genomes in the project directory, might be useful for later reference
  if !genomes_exclude.empty?
    fp = File.open(project_dir+"/excluded.txt","w")
    genomes_exclude.each {|ge| fp.puts(ge)}
    fp.close	
  end
  
  print "Processing NCBI data... "
  genomes_excluded = [] #this is excluded genomes and genomes_exclude is the genomes to exclude
  
  entries =  Dir.new(ncbi_dir).entries
  if entries.length==2
    #the ncbi directiry is empty and new data can be downloaded
    puts "The NCBI_PROCESSED_DIR is empty, it is possible to download data from NCBI."
    puts "I can automatically download Bacterial & Archael data (for more possibilities see INSTALL.txt)."
    puts "Download sequence data from NCBI? [Y/N] (default=Y, timeout 2 minutes)\n"
    answer = extra.get_user_input_timeout("",120,"y","...timeout...")[0].chr
    if answer.downcase=="y"
      print "\nDownload may take some time ..."
      #create temporary directory in the project directory and download data there
      Dir.mkdir(project_dir+"/tmp")
      success=system("wget -O #{project_dir}/tmp/all.gbk.tar.gz ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gbk.tar.gz")
      raise "\nError in downloading sequence data from NCBI.\n" if !success
      print "done\nProcessing ..."
      #unpack data
      success=system("tar xfz #{project_dir}/tmp/all.gbk.tar.gz -C #{project_dir}/tmp/")
      raise "\nError in unpacking the downloaded sequence data.\n" if !success
      #process the data and create the fasta files in ncbi_dir
      success=system("#{mydir}/process_ncbi.rb #{project_dir}/tmp/ #{ncbi_dir}")
      raise "\nError in processing the downloaded sequence data.\n" if !success
      #clean the dowloaded NCBI data
      #Dir.rmdir will not work here
      FileUtils.rm_rf("#{project_dir}/tmp")
    else
      puts "There is no training data available, provide some and run the program again."
      puts "Read INSTALL.txt for details on how this can be done."
      exit
    end
    
  end
  
  #get all the organism names from the files in the ncbi_dir
  #this can be used for genetaing generic clades
  organisms = []
  organism_file_map = Hash.new{|hash, key| hash[key] = Array.new;} #one organism can have multiple files
  n_sequences = 0
  Dir.new(ncbi_dir).entries.each do |e|
    next if e=="." || e==".."
    #check for valid extensions
    ext = e.split(".").last
    next if !extensions.nil? && !extensions.empty? && !extensions.member?(ext)

    organism = e.scan(/^(\d+).*/).flatten.first
    if organism.nil?
	    puts "WARN: invalid processed file: #{e}"
	    next
    end
    #exclude this genome if asked to
    if genomes_exclude.member? organism
      genomes_excluded << organism
      next
    end
    n_sequences += 1
    organisms << organism
    organism_file_map[organism] << e
  end
  print "(excluded #{genomes_excluded.length} sequences)" if !genomes_excluded.empty?
  puts "done"
  
  #tree processing
  nodes = []
  clades_file = ""
  #if tree is not provided then create clades list
  if tree_file.nil? || tree_file.empty?
    #create clades list
    descendants = Hash.new(0)
    organisms.each do |organism|
       taxonomy_ranks.each do |rank|
         parent = sqlite_taxonomy.parent_atrank(organism, rank).to_s  #changed I
	 descendants[parent]+=1
       end
    end
    descendants.each {|parent,dsc| nodes << parent if dsc>=n_min_genomes_generic}
    nodes=nodes.uniq
    nodes-=[nil]
    nodes-=[""]
  else
    #read the tree string
    fp = File.open(tree_file, "r")
      while(line=fp.gets)
        nodes << line.chomp if !line.nil?
      end
    fp.close
    tree_string = nodes[0]
    raise "\nFist line is nil in the tree file: #{tree_file}\n" if tree_string.nil?
  end

  #check if the tree is a newick string or a node list
  #assuming fixed tree_file name as "tree.newick" in the project directory
  if !(tree_string=~/;/) && nodes.length>=1
    print "Generating tree from the clades list (#{nodes.length} clades)... "
    clades_file = "#{project_dir}/clades.txt"
    fp = File.open(clades_file,"w")
    fp.print(nodes.join("\n"))
    fp.close
    #run script to create tree
    success=system("#{mydir}/ncbi2newick.rb #{clades_file} #{ncbi_tax_db} #{taxonomy_ranks.join(",")} > #{project_dir}/tree.newick")
    raise "\nTree creation from nodes list failed!\n" if !success
  else
    nodes = tree_string.split(/[^0-9]/).flatten
    print "Copying tree to the project directory... "
    fp = File.open(project_dir+"/tree.newick","w")
    fp.print(tree_string)
    fp.close
  end
  #change tree_file to the tree in the project directory
  #get back the tree_string  and nodes, this is necessary for further processing
  tree_file = "#{project_dir}/tree.newick"
  fp = File.open(tree_file, "r")
  tree_string = fp.gets
  fp.close
  raise "\nFist line is nil in the newick file: #{tree_file}\n" if tree_string.nil?
  nodes = tree_string.split(/[^0-9]/).flatten
  print "done\n"
  #processing for tree is done


  #map them on the tree
  print "Mapping genomes on the tree... "
  organism_tree_map = Hash.new(nil)
  tree_organism_map = Hash.new{|hash, key| hash[key] = Array.new;}
  mappings = [] #this makes some computations easier
  organisms_invalid = []
  organisms.each do |organism|
    #get a mapping
    next if !organism_tree_map[organism].nil? #the mapping already exists
    
    mapped = nil
    mapped = organism if nodes.member? organism
    #assuming that taxonomy_ranks starts with more specific ranks
    if mapped.nil?
      taxonomy_ranks.each do |rank|
        parent = sqlite_taxonomy.parent_atrank(organism, rank).to_s #changed I
        if nodes.member? parent and parent != ""
          mapped = parent
          break
        end
      end
    end
    raise "\nCould not map #{organism} on the tree\n" if mapped.nil?
    organisms_invalid << organism if mapped=="1"
    organism_tree_map[organism] = mapped
    tree_organism_map[mapped] << organism
    mappings << mapped
    if mapped == ""
      puts "Empty mapping #{organism} #{mapped}"
      exit
    end

  end #entries.each
  
  mappings = mappings.uniq
  organisms = organisms.uniq
  n_frag_per_node = (n_examples/mappings.length.to_f).ceil
 
  if !organisms_invalid.empty? 
    print "\n\tthese #{organisms_invalid.length} organisms will not be processed due to lack of mapping: \n\t#{organisms_invalid.join(",")}\n"
    fp=File.open(project_dir+"/organisms_not_used.txt","w")
    fp.puts organisms_invalid.join("\n")
    fp.close
  end
  
  #create directory structure in the project_dir
  directories = ["labels","sampled_fasta","train_data","models"]
  directories.each { |d| Dir.mkdir(project_dir+"/"+d) if !File.exists?(project_dir+"/"+d) }
  
  frag_len.each { |fl| Dir.mkdir(project_dir+"/sampled_fasta/"+fl)  }
  
  puts " (#{n_sequences} sequences) done"
  
  print "Generating sequence fragments... "
  $stdout.flush
  #chop and sample fasta
  #tree_organism_map contains node as keys and organism array as values
  tree_organism_map.each do |node,os|
    
    #get list of files for this node
    files = []
    os.each do |o| 
      next if organisms_invalid.member? o
      files << organism_file_map[o]
    end
    files = files.flatten
    next if files.empty?
    
    n_frag_per_file = (n_frag_per_node.to_f/files.length.to_f).ceil
    files.each do |file|
      organism = file.scan(/^(\d+).*/).flatten.first
      mapped = organism_tree_map[organism]
      fasta = Bio::FastaFormat.open(ncbi_dir+"/"+file)
      definition = ""
      seq_concat = ""
      fasta.each do |seq|
        definition = seq.definition if definition.empty?
        seq_concat << seq.data
      end
      seq_concat = seq_concat.gsub(/[^a-zA-Z]/,"")
      #remove non-ACGT to get a better normalized k-mer vector
      kmer_i = kmer.collect {|km| km.to_i}
      rep = "X" * kmer_i.max
      pat = "[^A,C,G,T]{#{kmer_i.max},}"
      pat = Regexp.new(pat, Regexp::IGNORECASE)
      seq_concat = seq_concat.gsub(pat,rep)
      
      frag_len.each do |fl|
        fl = fl.to_i
        next if seq_concat.length < fl
        #determine how many fragments to take
        n_frag = (seq_concat.length.to_f/fl.to_f).floor
        sampled_frag = (0...n_frag).to_a
        if n_frag > n_frag_per_file
          #choose some random fragments
          sampled_frag = sampled_frag.sort_by{rand}
          sampled_frag = sampled_frag[0...n_frag_per_file]
        end
        
        fp = File.open(project_dir+"/sampled_fasta/#{fl}/"+file,"a")
        
        n_frag_processed = 0
        start = 0
        step = fl #non-verlapping fragments
        dna = Bio::Sequence::NA.new(seq_concat)
        remainder = dna.window_search(fl, step) do |s|
          start = fl*n_frag_processed
          fp.puts s.to_fasta("#{definition}|at#{start}|label:#{mapped}", 1000) if sampled_frag.member?(n_frag_processed) #1000 characters per line
          n_frag_processed += 1
        end
        
        fp.close
        
      end #frag_len.each
      
    end #files.each
    
  end #tree_organism_map.each 
  
  puts " done"
  
  #deal with sample specific data if provided
  #this is kinda redundant and code can be reduced by possibly combining this with ncbi data processing
  genomes_excluded_ss = []
  n_sequences_ss = 0
  if !sample_specific_dir.nil? && !sample_specific_dir.empty?
    print "Processing sample specific data (all data will be used)... "
    entries = Dir.new(sample_specific_dir).entries
    entries.each do |e|
      next if e=="." || e==".."
      #check for valid extensions
      ext = e.split(/\./).last
      next if !extensions.nil? && !extensions.empty? && !extensions.member?(ext)
      
      organism = e.scan(/^(\d+).*/).flatten.first
      raise "\nInvalid sample specific file: #{e}\n\n" if organism.nil?
      
      #exclude this genome if asked to
      #not from sample specific directory
      #if genomes_exclude.member? organism
      #	genomes_excluded_ss << organism
      #       next
      #end
      
      mapped = nil
      #first try to map it directly to one of the nodes
      mapped = organism if nodes.member? organism
      
      #get mapping on the tree
      #assuming that taxonomy_ranks starts with more specific ranks
      if mapped.nil?
        taxonomy_ranks.each do |rank|
          parent = sqlite_taxonomy.parent_atrank(organism, rank).to_s
          if nodes.member? parent
            mapped = parent
            break
          end
        end
      end
      raise "\nCould not map #{organism} on the tree\n" if mapped.nil?
      if mapped=="1"
        print "\n\tSkipping #{organims} due to lack of mapping"
        next
      end
      
      #generate fragments
      fasta = Bio::FastaFormat.open(sample_specific_dir+"/"+e)
      definition = ""
      seq_concat = ""
      fasta.each do |seq|
        definition = seq.definition if definition.empty?
        seq_concat << seq.data
      end
      seq_concat = seq_concat.gsub(/[^a-zA-Z]/,"")
      
      frag_len.each_index do |ii|
        fl = frag_len[ii].to_i
        step = sample_specific_step[ii]
        step = fl if step.nil? || step==0
        if seq_concat.length < fl
          puts "No sample specific data for organism #{organism} at frag_len #{fl}"
          next
        end
        #take all fragments from sample specific data
        fp = File.open(project_dir+"/sampled_fasta/#{fl}/"+e,"a")
        n_frag_processed = 0
        start = 0
        step = step.to_i
        dna = Bio::Sequence::NA.new(seq_concat)
        remainder = dna.window_search(fl, step) do |s|
          start = fl*n_frag_processed
          fp.puts s.to_fasta("#{definition}|at#{start}|label:#{mapped}",1000)
          n_frag_processed += 1
        end
        
        fp.close
        
      end #frag_len.each
      
      n_sequences_ss+=1
      
    end
    puts "WARNING!!! no data processed in SAMPLE_SPECIFIC_DIR" if n_sequences_ss==0
  
  print "(excluded from SS #{genomes_excluded_ss.length})" if !genomes_excluded_ss.empty?
  print " (#{n_sequences_ss} SS sequences) done"

  end #sample specific
  
  print "\nGenerating k-mer features... "
  
  #now generate kmers
  #
  fasta2kmers_command = "#{mydir}/../bin/fasta2kmers -s 1 -o 1 -h 1 -l 1 -C #{rm_rev_complement} -t #{normalization} -r #{rev_complement}"
  fasta2kmers2_command = "#{mydir}/../bin/fasta2kmers2 -a a -s 1 -l 1 -o 1 -b 1 -R #{rm_rev_complement} -n #{normalization} -r #{rev_complement}"
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
  
  frag_len.each do |fl|
    print "#{fl} "
    entries = Dir.new(project_dir+"/sampled_fasta/#{fl}").entries
    pids = []
    entries.each do |e|
      next if e=="." || e==".."
      fasta2kmers_command_final = fasta2kmers_command + "  #{project_dir}/sampled_fasta/#{fl}/#{e} >> #{project_dir}/train_data/#{fl}.sl"
      fasta2kmers2_command_final = fasta2kmers2_command + " -i #{project_dir}/sampled_fasta/#{fl}/#{e} -f #{project_dir}/train_data/#{fl}.sl"
      fasta2kmers_command_use = fasta2kmers_command_final if !use_fasta2kmers2
      fasta2kmers_command_use = fasta2kmers2_command_final if use_fasta2kmers2 
      success=system(fasta2kmers_command_use)
      raise "\nError in creating k-mer features: #{fasta2kmers_command_use}\n" if !success
    end #entries
    
    #here the runtime is bounded by the max amount of data for one of the fragment lengths
  end #frag_len
  
  print "done"
else
   tree_file = "#{project_dir}/tree.newick" 
end #only_models

print "\nBuilding models (#{Time.now})..." 

#now as the training data is ready get the models
#if no grid was given then just build models
#kernel options
kernel_opt = "-t #{kernel} -g #{kernel_rbf_gamma} -d #{kernel_polynomial_degree} -s #{kernel_polynomial_s}"
loss_opt = "-l #{loss_function} --L #{loss_action}"
other_opt = "--z #{z_standardization} --v #{misc_nodes} --t #{tree_file}"
learn_command = "#{mydir}/../bin/svm_phylo_learn #{kernel_opt} #{loss_opt} #{other_opt} -v 1 -o 2"
cv_command    = "#{mydir}/../bin/svm_phylo_cv #{kernel_opt} #{loss_opt} #{other_opt} -x 3 -v 1 -o 2 --r 1 --S 1"
cv_command    = "svm_phylo_cv #{kernel_opt} #{loss_opt} #{other_opt} -x 3 -v 1 -o 2 --r 1" ###DELETE ME!!!!

learn_command << " --c 1" if balance_classes
cv_command << " --c 1" if balance_classes

if c_grid.length==1
  c_val = c_grid[0]
  frag_len.each do |fl|
    puts "\tBuilding #{fl} length model #{Time.now}"
    learn_command_final = learn_command + " -c #{c_val} #{project_dir}/train_data/#{fl}.sl #{project_dir}/models/#{fl}_c#{c_val}.svm"	
    if parallel_models
      pid = fork do
        exec(learn_command_final) if !fork.nil?
        exit 99 #this is just in case the forked process doesnt exit
      end
      Process.detach(pid)
    else
      success=system(learn_command_final)
      raise "\nError in learning the #{fl} fragment length model: #{learn_command_final}\n" if !success
    end
    puts "\t#{fl} length model ready #{Time.now}"
  end
else #do cross validation
   
  frag_len.each do |fl|
    puts "\n\tCross-validating #{fl} length model #{Time.now}"
    cv_loss = []
    cv_zeroone = []
    c_grid.each do |c_val|
      print "\t\tC=#{c_val} #{Time.now}..."
      cv_command_final = cv_command + " -c #{c_val} #{project_dir}/train_data/#{fl}.sl"
      cv_out = nil
      IO.popen cv_command_final do |cv_out|
        while(line = cv_out.gets)
          loss = line.scan(/Average loss in cross-validation: (.+)$/)
          cv_loss << loss.flatten.first.to_f if !loss.nil? && !loss.empty?
          zeroone = line.scan(/one-error in cross-validation: (.+)$/)
          cv_zeroone << zeroone.flatten.first.to_f if !zeroone.nil? && !zeroone.empty?
        end #while
      end #do
      puts " done #{Time.now}\n"
    end #c_grid

    if cv_loss.length != c_grid.length
      puts cv_loss.join("\t")
      raise "Error, something went wrong with cross-validation of #{fl} length fragment model. Quitting"
    end

    puts "C grid: " + c_grid.join("\t")
    puts "CV loss: " + cv_loss.join("\t")
    puts "CV  0-1: " + cv_zeroone.join("\t")
 
    #get the c_val with minimum loss
    loss_min = cv_loss.min
    cv_loss.each_index do |i|
      if cv_loss[i]==loss_min
        #build model
        c_val = c_grid[i]
        puts "\n\t#{fl} length model with C=#{c_val} and CV-loss=#{loss_min}"
        learn_command_final = learn_command + " -c #{c_val} #{project_dir}/train_data/#{fl}.sl #{project_dir}/models/#{fl}_c#{c_val}.svm"
        success=system(learn_command_final)
        raise "\nError in learning the model: #{learn_command_final}\n" if !success
        puts " done #{Time.now}\n"
	break
      end
    end
    
  end
end

puts "done #{Time.now}"

if clean_up_train
  print "\nCleaning up..."
  FileUtils.rm_rf("#{project_dir}/train_data")
  FileUtils.rm_rf("#{project_dir}/sampled_fasta")
  print "done"
end

sqlite_taxonomy.close

if parallel_models
  puts "!!!YOU NEED TO CHECK YOURSELF IF THE MODELS ARE READY!!!"
  puts "Model directory: #{project_dir}/models."
else
  puts "Processing finished (#{Time.now})...models are ready in #{project_dir}/models."
end
