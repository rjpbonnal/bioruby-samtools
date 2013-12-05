$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require 'rubygems'
require 'bio/db/sam'
require 'bio/db/alignment'
require 'pp'

@test_folder = "samples/small/"
@testTAMFile = @test_folder + "test.tam"
@testBAMFile = @test_folder + "testu.bam"
@testReference = @test_folder + "test_chr.fasta"
@outfile = @test_folder + "test.out"
    @sam = Bio::DB::Sam.new(
        :fasta => @testReference, 
        :bam => @testBAMFile
    )

#@sam.faidx()
#@sam.faidx(:chr=>'scaffold_1', :start=>1, :stop=>5928454)
#@sam.calmd(:u=>'true')
#@sam.calmd()
#@sam.phase(:b => "phased", :q=> 39)

@sam.fetch("chr_1", 10,1000) do |sam|
  pp sam.class
end