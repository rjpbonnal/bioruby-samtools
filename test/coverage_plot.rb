$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')

require "bio/db/sam"
require "bio/db/sam/sam"

bam_file = ARGV[0]
fasta_file = ARGV[1]
chromosmes = ARGV[2]


sam       = Bio::DB::Sam.new({:bam=>bam_file, :fasta=>fasta_file})
sam.open
File.open(chromosmes, "r") do |file|
    file.each_line{|line|
      
      
      
     
      fetching = line.split()[0]
	    last = line.split()[1].to_i
	    covs = sam.chromosome_coverage(fetching, 0, last)
	    puts "POS\tCOV"
	   # puts "#{fetching}\t#{last}"
	    covs.each_with_index{ |cov, i| puts "#{i}\t#{cov}" }
    }
end
sam.close