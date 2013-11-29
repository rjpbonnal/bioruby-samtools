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

#@sam.faidx("scaffold_1", 1, 5928454)
@sam.reheader(testSamFile)