$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require 'rubygems'
require 'bio/db/sam'
require "test/unit"


class TestBioDbSam < Test::Unit::TestCase
  
  def setup
    @test_folder                = "test/samples/small"
    @testTAMFile                = @test_folder + "/test.tam"
    @testBAMFile                = @test_folder + "/testu_fixedmates.bam"
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
      assert_equal(sam.class, Bio::DB::Alignment)
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
      assert_equal(sam.class, Bio::DB::Alignment)
    end
    #puts @sam.last_command
  end
  
  def test_fetch_with_function
    block = Proc.new {|a| assert_equal(a.class, Bio::DB::Alignment)}
    @sam.fetch_with_function("chr_1", 10,1000, &block)
  end
  
  def test_chromosome_coverage
    pp @sam.chromosome_coverage("chr_1", 322, 5)
  end
  
  def test_average_coverage
    #pp @sam.average_coverage("chr_1", 322, 5)
  end
   
  def test_index_stats
    @sam.index_stats.each_pair do |seq, stat|
      assert_send([['chr_1' , '*'], :member?, seq])
    end
  end
  
  def test_fetch_reference
    seq_expected = "CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
    seq_fetched = @sam.fetch_reference("chr_1", 1, 70, :as_bio => false)
    assert_equal(seq_fetched, seq_expected)
  end
  
  def test_sort
    @sam.sort(:prefix=>"test_sorted")
    sortedsam = "test_sorted.bam"
    @sortsam = Bio::DB::Sam.new(
        :fasta => @testReference, 
        :bam => sortedsam
    )
    pos = 0
    @sortsam.view()do |sam|
      assert(sam.pos > pos, "Not sorted by position")
      pos = sam.pos
    end
  end

  def test_mpileup
    @sam.mpileup(:g => false) do |pileup|
      assert_equal(pileup.class, Bio::DB::Pileup)
      assert_kind_of(Bio::DB::Pileup, pileup)
      assert_equal(pileup.ref_name, 'chr_1')
    end

    @sam.mpileup(:u => true) do |pileup|
      assert_kind_of(Bio::DB::Vcf, pileup)
      assert_equal(pileup.chrom, 'chr_1')
    end
  end

  def test_depth
    #the depth of coverage should be '1' at all positions
    @sam.depth(:r=>"chr_1:25-42") do |al|
    assert_equal(al[2].to_i, 1)
    end
  end

  def test_index
    @sam.index()
    test_bai_file = @testBAMFile+".bai"
    assert_nothing_thrown do
      File.open(test_bai_file, "r")
    end
    test_bai_file = @test_folder+"/different_index.bam.bai"
    @sam.index(:out_index=> test_bai_file)
    assert_nothing_thrown do
      File.open(test_bai_file, "r")
    end
  end

  def test_fixmate

  end

  def test_flagstat
     @sam.flag_stats
  end

  def test_merge

  end

  def test_rmdup

  end
end