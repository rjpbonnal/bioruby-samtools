$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require 'rubygems'
require 'bio/db/vcf'
require "test/unit"


class TestVcf < Test::Unit::TestCase
  
  def setup
    @vcf1 = Bio::DB::Vcf.new("19	111	.	A	C	9.6	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3",["a","b","c"]) #from a 3.3 vcf file
    @vcf2 = Bio::DB::Vcf.new("20	14370	rs6054257	G	A	29	0	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:-1,-1") #from a 3.3 vcf file
    @vcf3 = Bio::DB::Vcf.new("19	111	.	A	C	9.6	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3") #from a 4.0 vcf file
    @vcf4 = Bio::DB::Vcf.new("20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,") #from a 4.0 vcf file
  end
  
  
  def test_parse
    assert_equal("19", @vcf1.chrom)
    assert_equal(111, @vcf1.pos)
    assert_equal(nil, @vcf1.id)
    assert_equal("A", @vcf1.ref)
    assert_equal("C",@vcf1.alt)
    assert_equal(9.6,@vcf1.qual)
    assert_equal(nil, @vcf1.filter)
    assert_equal(nil, @vcf1.info)
    assert_equal({"a"=>{"GT"=>"0|0", "HQ"=>"10,10"},
                 "b"=>{"GT"=>"0|0", "HQ"=>"10,10"},
                 "c"=>{"GT"=>"0/1", "HQ"=>"3,3"}}, @vcf1.samples)

    assert_equal("20", @vcf2.chrom)
    assert_equal(14370, @vcf2.pos)
    assert_equal('rs6054257', @vcf2.id)
    assert_equal("G", @vcf2.ref)
    assert_equal("A",@vcf2.alt)
    assert_equal(29,@vcf2.qual)
    assert_equal("0", @vcf2.filter)
    assert_equal({"DP"=>"14", "AF"=>"0.5", "NS"=>"3", "DB"=>nil, "H2"=>nil}, @vcf2.info)
    assert_equal({"1"=>{"DP"=>"1", "GT"=>"0|0", "HQ"=>"51,51", "GQ"=>"48"},
     "2"=>{"DP"=>"8", "GT"=>"1|0", "HQ"=>"51,51", "GQ"=>"48"},
     "3"=>{"DP"=>"5", "GT"=>"1/1", "HQ"=>"-1,-1", "GQ"=>"43"}}, @vcf2.samples)

     assert_equal("19", @vcf3.chrom)
     assert_equal(111, @vcf3.pos)
     assert_equal(nil, @vcf3.id)
     assert_equal("A", @vcf3.ref)
     assert_equal("C",@vcf3.alt)
     assert_equal(9.6,@vcf3.qual)
     assert_equal(nil, @vcf3.filter)
     assert_equal(nil, @vcf3.info)
     assert_equal({"1"=>{"GT"=>"0|0", "HQ"=>"10,10"},
                  "2"=>{"GT"=>"0|0", "HQ"=>"10,10"},
                  "3"=>{"GT"=>"0/1", "HQ"=>"3,3"}}, @vcf3.samples)

     assert_equal("20", @vcf4.chrom)
     assert_equal(14370, @vcf4.pos)
     assert_equal('rs6054257', @vcf4.id)
     assert_equal("G", @vcf4.ref)
     assert_equal("A",@vcf4.alt)
     assert_equal(29,@vcf4.qual)
     assert_equal("PASS", @vcf4.filter)
     assert_equal({"DP"=>"14", "AF"=>"0.5", "NS"=>"3", "DB"=>nil, "H2"=>nil}, @vcf4.info)
     assert_equal({"1"=>{"DP"=>"1", "GT"=>"0|0", "HQ"=>"51,51", "GQ"=>"48"},
                   "2"=>{"DP"=>"8", "GT"=>"1|0", "HQ"=>"51,51", "GQ"=>"48"},
                   "3"=>{"DP"=>"5", "GT"=>"1/1", "HQ"=>".,", "GQ"=>"43"}}, @vcf4.samples)
  end
  
  def test_int_or_raw
    @vcf1.int_or_raw(1)
  end
  
end