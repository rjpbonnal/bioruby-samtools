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
    @sam = Bio::DB::Sam.new(
        :fasta => @testReference, 
        :bam => @testBAMFile
    )
  end
  
  def test_new
    assert_kind_of Bio::DB::Sam, @sam
  end
  
  def test_view
    #how to get Bio::DB::Alignment objects ..
    @sam.view() do |sam|
      pp sam
    end
    
    #how to get binary 
    #f = File.open("view.bam", "w")
    #@sam.view( :b => true ) do |line|
    #  print line
    #end
    #f.close
    
  end
  
  def test_fetch
    
  end
  
end