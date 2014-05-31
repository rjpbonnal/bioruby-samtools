$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require 'rubygems'
require 'bio/db/sam'
require "test/unit"
#gem 'ruby-prof'
gem 'test-unit'
#require "ruby-prof"


class TestBioDbSam < Test::Unit::TestCase
#  include RubyProf::Test
  
  class << self

    def shutdown
      File.delete("test/samples/small/different_index.bam.bai")
      File.delete("test/samples/small/dupes_rmdup.bam")
      File.delete("test/samples/small/mates_fixed.bam")
      File.delete("test/samples/small/reheader.bam")
      File.delete("test/samples/small/test_chr.fasta.fai")
      File.delete("test/samples/small/test_sorted.bam")
      File.delete("test/samples/small/maps_merged.bam")
      File.delete("test/samples/small/maps_cated.bam")
    end
  end
  
  
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
    assert_kind_of(Bio::DB::Sam, @sam)
  end
  
  def test_index
    test_bai_file = @testBAMFile+".bai"
    #test to see if the index file exists. If so, delete it
    if File.exist?(test_bai_file) == true
      puts "bam index exists....deleting..."
      File.delete(test_bai_file)
    end
 
    #No bam file 
    assert_equal(@sam.indexed?, false)
    #index the bam file
    @sam.index()
    assert_equal(@sam.indexed?, true)
    #make sure the .bai file exists
    assert_nothing_thrown do
      File.open(test_bai_file, "r")
    end
    assert(File.size(test_bai_file) > 0, "From test_index: .bai file is empty")
    #as above, but give the output a different name
    test_bai_file = @test_folder+"/different_index.bam.bai"
    @sam.index(:out_index=> test_bai_file)
    assert_nothing_thrown do
      File.open(test_bai_file, "r")
    end
    assert(File.size(test_bai_file) > 0, "From test_index: .bai file is empty")
  end

  def test_view
    #how to get Bio::DB::Alignment objects ..
    @sam.view() do |sam|
      #test that all the objects are Bio::DB::Alignment objects and their reference is 'chr_1'
      assert_equal(sam.class, Bio::DB::Alignment)
      assert_equal(sam.rname, "chr_1")
    end
  end
  
  def test_fetch
#puts    @sam.inspect
    i = 0
    @sam.index
    @sam.fetch("chr_1", 10,1000) do |sam|
      #test that all the objects are Bio::DB::Alignment objects
      assert_equal(sam.class, Bio::DB::Alignment)
      assert_equal(sam.rname, "chr_1")
      i += 1
    end
    assert(i>0)
    assert_equal(i,9)
    
  end
  
  def test_fetch_with_function
    #pass the assert to method
    block = Proc.new {|a| assert_equal(a.class, Bio::DB::Alignment)}
    @sam.fetch_with_function("chr_1", 10,1000, &block)
  end
  
  def test_chromosome_coverage
    #the coverage should only be 1.0 or 2.0
    cov = @sam.chromosome_coverage("chr_1", 33, 19)
    cov.each do |pu|
      assert_send([[1.0 , 2.0], :member?, pu])
    end
  end
  
  def test_average_coverage
    #there should be 10 positions with cov of 1.0 and 10 with cov of 2.0, so average of 1.5
    test_bai_file = @testBAMFile+".bai"
    if File.exist?(test_bai_file) == false
      @sam.index()
    end
    avcov = @sam.average_coverage("chr_1", 33, 19)
    assert_equal(avcov, 1.5)
    File.delete(test_bai_file)
  end
  
  def test_faidx
    @sam.faidx()
    test_fai_file = @testReference+".fai"
    #test that the .fai file exists
    assert_nothing_thrown do
      File.open(test_fai_file, "r")
    end
    #test that the file is not empty
    assert(File.size(test_fai_file) > 0, "From test_faidx: .fai file is empty")
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
  
  def test_reheader
    sam_header = @test_folder + "/map_for_reheader.sam"
    outfile = @test_folder + "/reheader.bam"
    
    @sam.reheader(sam_header, :out=>outfile)
    reheader_bam = Bio::DB::Sam.new(:fasta => @testReference, :bam => outfile)
    #check that the reference is 'chr_2'
    reheader_bam.view()do |sam|
      assert_equal(sam.rname, "chr_2")
    end
  end

  def test_calmd
    no_md_sam = @test_folder + "/no_md.sam"
    md = Bio::DB::Sam.new(:fasta => @testReference, :bam => no_md_sam)
    block = Proc.new {|a| assert(a.tags.has_key?('MD'), "From test_calmd: couldn't find the MD tag")}
    md.calmd(:S=>true, &block)
    
  end
  
  def test_mpileup
    #create an mpileup
  #  @sam.index
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
  
  def test_mpileup_reg
    #create an mpileup
    reg = Bio::DB::Fasta::Region.new
    reg.entry = "chr_1"
    reg.start = 1
    reg.end = 334
    
    @sam.mpileup_cached(:r=>reg,:g => false, :min_cov => 1, :min_per =>0.2) do |pileup|
      #test that all the objects are Bio::DB::Pileup objects
      assert_kind_of(Bio::DB::Pileup, pileup)
      #test that the reference name is 'chr_1' for all objects
      #puts pileup
      assert_equal(pileup.ref_name, 'chr_1')
    end
  
    region = @sam.cached_regions[reg.to_s]
    #puts "cahced_region: #{region.inspect}"
    puts "AVG COV: #{region.average_coverage}"
    puts "Reference: #{region.reference}"
    puts "Consensus: #{region.consensus}"
    puts "called: #{region.called}"
    #, :snps, :reference, :base_ratios, :consensus, :coverages
    snps_tot = Bio::Sequence.snps_between(region.reference, region.consensus)
    assert_equal(snps_tot, 5)
    assert_equal(region.called, 213)
  end
  
  def test_mpileup_reg_05
    #create an mpileup
    reg = Bio::DB::Fasta::Region.new
    reg.entry = "chr_1"
    reg.start = 1
    reg.end = 334
    @sam.mpileup_cached(:r=>reg,:g => false, :min_cov => 1, :min_per =>0.4) do |pileup|
      #test that all the objects are Bio::DB::Pileup objects
      assert_kind_of(Bio::DB::Pileup, pileup)
      #test that the reference name is 'chr_1' for all objects
      #puts pileup
      assert_equal(pileup.ref_name, 'chr_1')
  
    end
  
    region = @sam.cached_regions[reg.to_s]
    
    #, :snps, :reference, :base_ratios, :consensus, :coverages
    snps_tot = Bio::Sequence.snps_between(region.reference, region.consensus)
    assert_equal(snps_tot, 1)
    assert_equal(region.called, 213)
  end

  def test_depth
    #the depth of coverage should be '1' at all given positions
    @sam.depth(:r=>"chr_1:25-42") do |al|
      assert_equal(al[2].to_i, 1)
    end
  end

  def test_fixmate
    mates_fixed_bam = @test_folder + "/mates_fixed.bam"
    @sam.fix_mates(:out_bam=>mates_fixed_bam)
    assert_nothing_thrown do
      File.open(mates_fixed_bam, "r")
    end
    assert(File.size(mates_fixed_bam) > 0, "From test_fixmate: .bam file is empty")
  end

  def test_flagstats
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
    File.delete merged_bam_file if File.exists?(merged_bam_file)
#    File.delete("test/samples/small/maps_merged.bam")
    @sam.merge(:out=>merged_bam_file, :bams=>bam_files, :n=>true)
    merged_bam = Bio::DB::Sam.new(:fasta => @testReference, :bam => merged_bam_file)
    no_reads_mapped = 0;
    
    merged_bam.view() do |al|
      assert_kind_of(Bio::DB::Alignment, al)
      no_reads_mapped+=1
    end
    assert_equal(no_reads_mapped, 10)
  end
  
  def test_cat
    #same files used for merge, but we'll cat them instead
    bam1 = @test_folder + "/map_to_merge1.bam"
    bam2 = @test_folder + "/map_to_merge2.bam"
    
    bam_files = [bam1, bam2]
    
    cat_bam_file = @test_folder + "/maps_cated.bam"
    File.delete cat_bam_file if File.exists?(cat_bam_file)
    @sam.merge(:out=>cat_bam_file, :bams=>bam_files)
    cated_bam = Bio::DB::Sam.new(:fasta => @testReference, :bam => cat_bam_file)
    
    no_reads_mapped = 0;
    cated_bam.view() do |al|
      assert_kind_of(Bio::DB::Alignment, al)
      no_reads_mapped+=1
    end
    #there should be 10 reads in the cat'd maps
    assert_equal(no_reads_mapped, 10)
  end

  def test_rmdup
    #dupes contains 4 reads mapped once and one read mapped to the same place 268 times.
    dupes = @test_folder + "/dupes.bam"
    unduped = @test_folder + "/dupes_rmdup.bam"
    bam_with_dupes = Bio::DB::Sam.new(:fasta => @testReference, :bam => dupes)
    bam_with_dupes.remove_duplicates(:s=>true, :out=>unduped)

    unduped_bam = Bio::DB::Sam.new(:fasta => @testReference, :bam => unduped)
    #rmdup should remove 267 of the 268 reads mapping to the same place, so producing a bam file with 5 reads
    readcount = 0
    unduped_bam.view()do |sam|
        readcount +=1  
    end
    assert_equal(readcount, 5)
  end
  
  def test_targetcut
    sorted_bam = @test_folder + "/sorted.bam"
    cut = Bio::DB::Sam.new(:fasta => @testReference, :bam => sorted_bam)
    assert_nothing_thrown do
      cut.targetcut
    end
  end
  
  def test_docs
    #force an error (use 'samtool' instead of 'samtools')
    output = Bio::DB::Sam.docs('samtool', 'tview')
    assert_equal(output, "program must be 'samtools' or 'bcftools'")
  end   
end