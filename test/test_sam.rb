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
      #pp sam
    end
    
    #how to get binary 
    #f = File.open("view.bam", "w")
    #@sam.view( :b => true ) do |line|
    #  print line
    #end
    #f.close
    
  end
  
  def test_fetch
    @sam.fetch("chr_1", 10,1000) do |sam|
      #pp sam
    end
    #puts @sam.last_command
  end
  
  def test_fetch_with_function
    #block = Proc.new {|a| pp a}
    #@sam.fetch_with_function("chr_1", 10,1000, &block)
  end
  
  def test_chromosome_coverage
    
  end
  
  def test_pileup
    @sam.mpileup(:u => true) do |pileup|
      pp pileup
    end
    puts @sam.last_command
  end
  
end