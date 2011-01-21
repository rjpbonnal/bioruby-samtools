$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')

require "bio/db/sam"
require "bio/db/sam/sam"

bam_file = ARGV[0]
fasta_file = ARGV[1]
chromosmes = ARGV[2]
i = 0
sam       = Bio::DB::Sam.new({:bam=>bam_file, :fasta=>fasta_file})
sam.open
File.open(chromosmes, "r") do |file|
  file.each_line{|line|
    $stderr.puts i.to_s + " done" if i%10 == 0
    ObjectSpace.each_object(Bio::DB::Sam) {|x|  $stderr.puts x } if i%10 == 0
    count =  ObjectSpace.each_object() {|x|  } if i%10 == 0
    $stderr.puts count if i%10 == 0
    i = i + 1
    fetching = line.split(' ')[0]
	  
	  cov = sam.average_coverage(fetching, 0, 16000)
	  puts fetching + " " + cov.to_s
  }
end
sam.close