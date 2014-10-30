#!/usr/bin/env ruby

require 'bio'

require 'optparse'
require 'set'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bio-samtools.rb')
require path


def log(msg)
  time=Time.now.strftime("%Y-%m-%d %H:%M:%S.%L")
  puts "#{time}: #{msg}"
end


options = {}
options[:min_cov] = 5
options[:min_percentage] = 0.5
options[:output_file] = "-"
OptionParser.new do |opts|
  
  opts.banner = "Usage: bam_consensus.rb [options]"

  opts.on("-b", "--bam_file FILE", "BAM File to call for the consensus") do |o|
    options[:bam] = o
  end
  opts.on("-r", "--reference FASTA", "reference with the contigs") do |o|
    options[:reference] = o
  end
  opts.on("-p", "--min_percentage FLOAT", "Minimum percentage to call for a base. When more than one base gets the percentage, an ambiguty code is produced") do |o|
    options[:min_percentage] = o.to_f / 100
  end

   opts.on("-m", "--min_cov INT", "Minimum percentage to call for a base. When more than one base gets the percentage, an ambiguty code is produced") do |o|
    options[:min_cov] = o.to_i 
  end

  opts.on("-f", "--filter_entries FILE", "File with a list of entries to process") do |o|
    options[:filter_entries] = o
  end

  opts.on("-o" "--output_file FILE", "Output of the program, in fasta format") do |o|
    options[:output_file] = o
  end
  
end.parse!

bam =  @parental_1_sam =  Bio::DB::Sam.new({:fasta=>options[:reference], :bam=>options[:bam]})
region_set = nil
  if options[:filter_entries]
  region_set = Set.new
  File.foreach(options[:filter_entries]) do |line| 
    region_set << line.chomp
  end
end

fasta_db = Bio::DB::Fasta::FastaFile.new(:fasta=> options[:reference])
fasta_db.load_fai_entries

output = $stdout


output = File.open(options[:output_file], "w") if options[:output_file] != "-"

fasta_db.index.entries.each do | r |
  process = true
	region=r.get_full_region

  process = region_set.include? region.entry if region_set
  if process
    reg = bam.fetch_region({:region=>region,  :min_cov=>options[:min_cov],:min_per=>options[:min_percentage],  :A=>1})
    cons = reg.consensus
    org = fasta_db.fetch_sequence(region)
    if cons.upcase != org.upcase 
      output.puts ">#{region.entry}"
      tmp = cons.scan /.{1,80}/
      output.puts  tmp.join("\n")
    end
  end
end

output.close if options[:output_file] != "-"