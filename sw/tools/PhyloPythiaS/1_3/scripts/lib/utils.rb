# utilities
# some requires

mydir = File.dirname(__FILE__)

$LOAD_PATH.unshift("#{mydir}/rubygems/bioruby/lib")
require "bio"

require "timeout"

class File_util

	def sample_lines(file_name, n_lines)
		fp = File.open(file_name, "r")
		lines = []
		while(line=fp.gets)
			next if line.nil?
			lines << line
		end
		fp.close

		lines_to_sammple = (0...lines.length).to_a
		lines_to_sammple = (0...lines.length).to_a.sort_by{rand}[0...n_lines] if n_lines < lines.length 

		extracted = []
		lines_to_sammple.each { |i| extracted << lines[i] }

		extracted
	
	end

end


# extend ruby's array class to calculate mean etc
class Array
	def sum
		inject( 0 ) { |sum,x| sum+x }
	end

	def sumSq
		inject(0) {|r,i| r + i**2}
	end
	
	def mean
		(size > 0) ? sum.to_f / size : 0
	end


	#ivariance, standatd deviation
	def sd
		n=self.length.to_f
		mean=self.mean
		sum=0
		self.each { |d| sum += (d-mean)**2 }
		var = sum / (n-1).to_f
		sd = Math.sqrt(var)
		sd
	end

	def z_score
		mean = self.mean
		sd = self.sd
		sd = 1 if sd==0
		z = []
		self.each {|v| z << (v-mean)/sd }
		z
	end

	def pearson(y)
   		n=self.length 

   		sumx=self.sum
   		sumy=y.sum

   		sumxSq=self.sumSq
   		sumySq=y.sumSq

   		prods=[]; self.each_with_index{|this_x,i| prods << this_x*y[i]}
   		pSum=prods.sum

   		# Calculate Pearson score 
   		num=pSum-(sumx*sumy/n) 
   		den=((sumxSq-(sumx**2)/n)*(sumySq-(sumy**2)/n))**0.5 
   		if den==0
     			return 0 
   		end
   		r=num/den 
   		return r 
 	end

	def count
		k=Hash.new(0)
		self.each{|x| k[x]+=1 }
		k
	end

	# Divides an array into chunks of the same size.
        def chunk(pieces)
                len = self.length;
                mid = (len/pieces)
                chunks = []
                start = 0
                1.upto(pieces) do |i|
                        last = start+mid
                        last = last-1 unless len%pieces >= i
                        chunks << self[start..last] || []
                        start = last+1
                end
                chunks
        end


end

#compare strings/arrays: returns number of similar characters
module Cmpr
	def cmpr(str,start = 0)
		#start tells where to start comparison
		cmp = 0
		(start...self.length).each do |i|
			cmp += 1 if self[i]==str[i]			
		end
		return cmp
	end
end

#extend string and array with rotate
module Rotate
	def rotate(steps = 1)
		steps = steps % self.length
		self[-steps..-1].concat self[0...-steps]
	end
end

class String
	include Rotate
	include Cmpr
end

class Array
	include Rotate
	include Cmpr
end

class NCBI
	#processings on ncbi data
	#add functionality to identify phage source etc.
	def parse_gbk_gbff(file)
		#given a gbk/gbff file returns SOURCE and all LOCI concatenated
		fp = File.open(file,"r")
		raise "ERROR: Could not open file to read #{file}" if fp.nil?
		taxon = []
		definition = nil
		accession = []
		gi = []
		seq = []
		source = nil #there must be a single source
		#read entire file at once
		n_seq=0
		while(line=fp.gets)
			next if line.nil?
			line.chomp!

			if !(match = line.scan(/^DEFINITION\s+(.+)/) ).empty?
				definition = match[0]
			elsif line =~ /^VERSION/
				accession << line.scan(/^VERSION\s+([A-Za-z0-9\.]+)/)[0]
				gi << line.scan(/.+GI:(.+)/)[0]
			elsif !(match = line.scan(/^SOURCE\s+(.+)/) ).empty?
				raise "Different sources" if !source.nil? && source!=match[0]
				source = match[0]
			elsif !(match = line.scan(/.+db_xref.+taxon:(\d+).+/) ).empty?
				taxon << match[0]
			elsif line =~ /^ORIGIN/
				#this is the sequence
				seq[n_seq] = ""
				while(line2=fp.gets)
					break if line2 =~ /^\/\//o
					#dont use += its hopelessly slow
					#http://sob.apotheon.org/?p=309
					seq[n_seq] << line2
				end #while
				seq[n_seq].gsub!(/[^a-zA-Z]/,"")
				n_seq+=1
			end

			
		end	#while

		raise "\nSomething went wrong in parsing" if gi.length!=seq.length

		#create hash and return
		foo = Hash.new(nil)
		foo = {
			"definition" => definition,
			"source" => source,
			"taxon" => taxon.flatten,
			"accession" => accession.flatten,
			"gi" => gi.flatten,
			"seq" => seq.flatten
		}
		return foo

	end
end #class NCBI


class Extra

	def read_stdin
		lines = []
		unless STDIN.tty?   #we are in a pipeline
		        while((line = STDIN.gets))
				next if line.nil?
				line = line.chomp
				lines << line
			end
		end
		return lines
	end

	def get_user_input(string)
		print string
		answer = STDIN.gets
		STDOUT.flush
		answer
	end

	def get_user_input_timeout(string,time,default,timeout_text="")
		#string is the question you want to ask
		#time is timeout in seconds
		#default is the default return value
		#timeout_text is something to print after timeout
		print string
		answer=default
		begin
                	answer = Timeout::timeout(time) do
				STDIN.gets
			end
		rescue Timeout::Error
			print timeout_text
			answer = default
		end
		answer=default if answer==10 || answer=="\n" #10 is ascii for \n
		STDOUT.flush
		return answer
	end


end

