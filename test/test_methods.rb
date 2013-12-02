$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require 'rubygems'
require 'bio/db/sam'

testReference = "/Users/ethering/temp/samtools_data/contigs.fa"
testBAMFile = "/Users/ethering/temp/samtools_data/Map.bam"
testSamFile = "/Users/ethering/temp/samtools_data/Map2.sam"

@sam = Bio::DB::Sam.new(
        :fasta => testReference, 
        :bam => testBAMFile
    )
#opts={:chr=>'scaffold_1', :start=>1, :stop=>5928454}
opts={}
@sam.faidx(opts)
#@sam.calmd(:C => '30')
