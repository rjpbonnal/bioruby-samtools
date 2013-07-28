$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require 'rubygems'
require 'bio/db/sam'
require "test/unit"


class TestBioDbSam < Test::Unit::TestCase
  
  def setup
    @test_folder                = "test/samples/small"
    @testTAMFile                = @test_folder + "/test.tam"
    @testBAMFile                = @test_folder + "/testu.bam"
    @testReference              = @test_folder + "/test_chr.fasta"
  end
  
  def test_new
    sam = Bio::DB::Sam.new(
        :fasta => @testReference, 
        :bam => @testBAMFile
    )
    assert_kind_of Bio::DB::Sam, sam
  end
  
end