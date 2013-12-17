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
      #test that all the objects are Bio::DB::Alignment objects
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
      #test that all the objects are Bio::DB::Alignment objects
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
    #this is the first 70 nucleotides of the test seqeunce
    seq_expected = "CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
    #fetch the first 70 nuclotides
    seq_fetched = @sam.fetch_reference("chr_1", 1, 70, :as_bio => false)
    #test they're the same
    assert_equal(seq_fetched, seq_expected)
  end
  
  def test_sort
    #sort the bam file
    sortedsam = @test_folder + "/test_sorted.bam"
    @sam.sort(:prefix=>@test_folder + "/test_sorted")
    #create a new Bio::DB::Sam from the sorted bam
    @sortsam = Bio::DB::Sam.new(
        :fasta => @testReference, 
        :bam => sortedsam
    )
    pos = 0
    #iterate over the sorted sam file and make sure that the it's sorted by checking the order of the start positions for each read.
    @sortsam.view()do |sam|
      assert(sam.pos > pos, "Not sorted by position")
      pos = sam.pos
    end
  end

  def test_mpileup
    #create an mpileup
    @sam.mpileup(:g => false) do |pileup|
      #test that all the objects are Bio::DB::Pileup objects
      assert_kind_of(Bio::DB::Pileup, pileup)
      #test that the reference name is 'chr_1' for all objects
      assert_equal(pileup.ref_name, 'chr_1')
    end
    #do the same for Vcf output
    @sam.mpileup(:u => true) do |pileup|
      assert_kind_of(Bio::DB::Vcf, pileup)
      assert_equal(pileup.chrom, 'chr_1')
    end
  end

  def test_depth
    #the depth of coverage should be '1' at all given positions
    @sam.depth(:r=>"chr_1:25-42") do |al|
      assert_equal(al[2].to_i, 1)
    end
  end

  def test_index
    #index the bam file
    @sam.index()
    test_bai_file = @testBAMFile+".bai"
    #make sure the .bai file exists
    assert_nothing_thrown do
      File.open(test_bai_file, "r")
    end
    #as above, but give the output a different name
    test_bai_file = @test_folder+"/different_index.bam.bai"
    @sam.index(:out_index=> test_bai_file)
    assert_nothing_thrown do
      File.open(test_bai_file, "r")
    end
  end

  def test_fixmate
    mates_fixed_bam = @test_folder + "/mates_fixed.bam"
    @sam.fix_mates(:out_bam=>mates_fixed_bam)
    assert_nothing_thrown do
      File.open(mates_fixed_bam, "r")
    end
    
  end

  def test_flagstat
    #get the stats
    stats = @sam.flag_stats()
    #the number of reads mapped will be the first character on the first line.
    no_reads_mapped = stats[0][0].to_i
    #check that it's '9'
    assert_equal(no_reads_mapped, 9)
  end

  def test_merge
    bam1 = @test_folder + "/map_to_merge1.bam"
    bam2 = @test_folder + "/map_to_merge2.bam"
    bam_to_merge1 = Bio::DB::Sam.new(:fasta => @testReference, :bam => bam1)
    bam_to_merge2 = Bio::DB::Sam.new(:fasta => @testReference, :bam => bam2)
    bam_files = [bam_to_merge1, bam_to_merge2]
    
    merged_bam_file = @test_folder + "/maps_merged.bam"
    
    @sam.merge(:out=>merged_bam_file, :bams=>bam_files, :n=>true)
    merged_bam = Bio::DB::Sam.new(:fasta => @testReference, :bam => merged_bam_file)
    no_reads_mapped = 0;
    
    merged_bam.view() do |al|
      assert_kind_of(Bio::DB::Alignment, al)
      no_reads_mapped+=1
    end
    assert_equal(no_reads_mapped, 10)
  end

  def test_rmdup

  end
end