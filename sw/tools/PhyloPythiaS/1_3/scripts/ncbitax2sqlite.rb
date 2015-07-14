#!/usr/bin/ruby
#create SQLite database from NCBI taxonomy dumps

mydir = File.dirname(__FILE__)

$LOAD_PATH.unshift(mydir+"/lib/rubygems/sqlite3-ruby-1.2.5/lib")
require "sqlite3"

require "#{mydir}/lib/utils.rb"
extra = Extra.new

Usage = "\nncbitax2sqlite.rb ncbi_dmp_dir sqlite_db [schema_file]\n"

ncbi_dmp_dir=ARGV[0]
sqlite_db=ARGV[1]
schema_file = mydir+"/sqlite_tax_schema.txt"
schema_file = ARGV[2] if !ARGV[2].nil?

raise Usage if ncbi_dmp_dir.nil? || sqlite_db.nil?
raise "\nThe schema file does not exist\n" if !File.exists?(schema_file)

if !File.exists?(ncbi_dmp_dir+"/nodes.dmp") || !File.exists?(ncbi_dmp_dir+"/names.dmp")
	print "\nNCBI taxonomy dump files are not present."
	print "\nDownload taxonomy data from NCBI? [Y/N] (default=Y, timeout 2 minutes)\n"
	answer = extra.get_user_input_timeout("",120,"y","...timeout...")[0].chr
	if answer.downcase=="y"
		success=system("wget -O #{ncbi_dmp_dir}/taxdump.tar.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
		raise "\nDownloading failed\n" if !success
		success=system("tar xfz #{ncbi_dmp_dir}/taxdump.tar.gz -C #{ncbi_dmp_dir}")
		raise "\nUnpacking failed\n" if !success
		#delete the tar
		File.delete("#{ncbi_dmp_dir}/taxdump.tar.gz")
	else
		print "\nPlease provide correct directory name with NCBI taxonomy dump files.\n"
		exit
	end
end

if File.exists?(sqlite_db)
	#check if database is corrupted	
	db = SQLite3::Database.new(sqlite_db)
	db_check_passed=FALSE
	begin
		res1=db.execute("SELECT name FROM taxon_name WHERE taxon_id=?",1)
		res2=db.execute("SELECT ncbi_taxon_id FROM taxon WHERE taxon_id=?",1)
		db_check_passed=TRUE if res1.length!=0 && res2.length!=0
	rescue SQLite3::NotADatabaseException
		#this means the database is not a valid database
		db_check_passed=FALSE
	end
	db.close

	print "\nThe SQLite database file '#{sqlite_db}' exists"
	print "\nthis file passed the check and is a valid database (possibly outdated)" if db_check_passed
	print "\nthis file did not pass the check and is corrupted" if !db_check_passed
	print "\nRemove it? [Y/N] (default=Y, timeout 2 minutes)\n"
	answer = extra.get_user_input_timeout("",120,"y","...timeout...")[0].chr
	if answer.downcase=="y"
		File.delete(sqlite_db)
	else
		if db_check_passed
			puts "The database will not be touched"
			exit
		else
			raise "Please run again with valid options, quiting"
		end
	end

end

#create a new database
db = SQLite3::Database.new(sqlite_db)

#create the tables from the schema
schema = File.readlines(schema_file).join("")
db.execute_batch(schema)

print "Uploading nodes.dmp ..."

#read the ncbi dumps and populate the database
#first load the nodes.dmp as we will need it for names.dmp
#using parameterized queries for performance
#examples at http://rosettacode.org/wiki/Parametrized_SQL_statement
fp = File.open("#{ncbi_dmp_dir}/nodes.dmp")
n_line=0
insert_taxon = "INSERT INTO taxon (ncbi_taxon_id,parent_taxon_id,node_rank) values (?,?,?)"
db.transaction
while(line=fp.gets)
	
	next if line.nil?
	n_line+=1
	line = line.chomp.split(/\s*\|\s*/)
	
	#we are using only first three columns here
	db.execute(insert_taxon, line[0...3])
end
fp.close
db.commit

raise "\nThe previous transaction is active\n" if db.transaction_active?
print "done\nUploading names.dmp ..."

#load the names.dmp file
fp = File.open("#{ncbi_dmp_dir}/names.dmp")
n_line=0
insert_taxon_name = "INSERT INTO taxon_name (taxon_id,name,name_class) values (?,?,?)"
db.transaction
while(line=fp.gets)
	next if line.nil?
	line = line.chomp.split(/\s*\|\s*/)
	ncbi_taxon_id = line[0]
	taxon_id = db.execute("SELECT taxon_id from taxon WHERE ncbi_taxon_id=?",ncbi_taxon_id)
	line[0] = taxon_id[0]
	#only the following entries are used (2nd entry is "unique name" and not necessary)
	line=[line[0],line[1],line[3]]
	db.execute(insert_taxon_name, line)
end
fp.close
db.commit

print "done\nThe database is ready to use."

