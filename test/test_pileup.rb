$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')

require "bio/db/sam/pileup"
require "test/unit"

class TestPileup < Test::Unit::TestCase
  
  def setup
    @six_col = Pileup.new("seq1\t272\tT\t24\t,.$.....,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&")
    @ten_col = Pileup.new("seq2\t151\tG\tG\t36\t0\t99\t12\t...........A\t:9<;;7=<<<<<")
    @snp = Pileup.new("seq1\t272\tT\t24\t,.$.....,,.gGgGgGgGgGgGg^+.\t<<<+;<<<<<<<<<<<=<;<;7<&") 
    @snp_2 = Pileup.new("seq1\t272\tT\t24\t......aaaaaaggggggcccccc$^+\t<<<+;<<<<<<<<<<<=<;<;7<&")
    @snp_3 = Pileup.new("seq1\t272\tT\t24\t......aaaaaaaagggggccccc$^+\t<<<+;<<<<<<<<<<<=<;<;7<&")
  end
  
  def test_new_from_6_column
    assert_equal("seq1", @six_col.ref_name)
    assert_equal(272, @six_col.pos)
    assert_equal("T", @six_col.ref_base)
    assert_equal(24, @six_col.coverage)
    assert_equal(",.$.....,,.,.,...,,,.,..^+.", @six_col.read_bases)
    assert_equal("<<<+;<<<<<<<<<<<=<;<;7<&", @six_col.read_quals)
  end
  
  def test_new_from_10_column
    assert_equal("seq2", @ten_col.ref_name)
    assert_equal(151, @ten_col.pos)
    assert_equal("G", @ten_col.ref_base)
    assert_equal("G", @ten_col.consensus)
    assert_equal(36, @ten_col.consensus_quality)
    assert_equal(0, @ten_col.snp_quality)
    assert_equal(99, @ten_col.rms_mapq)
    assert_equal(12, @ten_col.coverage)
    assert_equal("...........A", @ten_col.read_bases)
    assert_equal(":9<;;7=<<<<<", @ten_col.read_quals)
  end
  
  def test_non_refs
    assert_equal({:A => 1, :C => 0, :T => 0, :G => 0}, @ten_col.non_refs)
    assert_equal({:A => 0, :C => 0, :T => 0, :G => 0}, @six_col.non_refs)
  end
  
  def test_consensus
    assert_equal("G", @snp.consensus)
    assert_equal("ACGT", @snp_2.consensus)
    assert_equal("A", @snp_3.consensus)
  end
  
  def test_non_ref_count
    assert_equal(13,@snp.non_ref_count)
    assert_equal(18,@snp_2.non_ref_count)
    assert_equal(18,@snp_3.non_ref_count)
  end
  
  def test_ref_count
    assert_equal(11,@snp.ref_count)
    assert_equal(6,@snp_2.ref_count)
    assert_equal(6,@snp_3.ref_count)
  end
  
  def test_ref_plus_non_ref_equal_to_coverage
    assert_equal(@snp.coverage,@snp.ref_count + @snp.non_ref_count)
   assert_equal(@snp_2.coverage,@snp_2.ref_count + @snp_2.non_ref_count)
   assert_equal(@snp_3.coverage,@snp_3.ref_count + @snp_3.non_ref_count)
  end
 
end
