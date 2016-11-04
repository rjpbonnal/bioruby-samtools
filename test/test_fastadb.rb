$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')

require 'rubygems'
require 'bio/db/sam'
require "test/unit"
#gem 'ruby-prof'
gem 'test-unit'
#require "ruby-prof"


class TestBioDbfastaDB < Test::Unit::TestCase

	def setup
		@test_folder                = "test/samples/small"
		@testReference  = @test_folder + "/test_fastadb.fasta"
		@fasta_local = Bio::DB::Fasta::FastaFile.new(
        	fasta: @testReference,
        	samtools: false
    	)

    	@fasta_samtools = Bio::DB::Fasta::FastaFile.new(
        	:fasta => @testReference,
        	samtools: true
    	)
	end
     class << self    
        def shutdown
            File.delete("test/samples/small/test_fastadb.fasta.fai")
        end
    end

	def test_faidx
		@fasta_samtools.index()
    	test_fai_file = @testReference+".fai"
    	#test that the .fai file exists
    	assert_nothing_thrown do
      		File.open(test_fai_file, "r")
    	end
    	#test that the file is not empty
    	assert(File.size(test_fai_file) > 0, "From test_faidx: .fai file is empty")
	end



	def test_fetch_reference
    	#this is the first 70 nucleotides of the test seqeunce
    	@fasta_samtools.faidx()

    	test_regions = { "chr_1:1-70" => "CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA", 
    		"chr_1:68-74" => "CTAACCC",
    		"chr.2:300-450" => "GGTGCAGTCACGGCTCGCAGTTATACTCCGGGAGAATGGAAAAGATTGTCCAAGGACCAA\
CAGGAAAAAGTGCGATCGCTGCGTAATAAAAAGAAGCAAGGAG\
GGAAACCCGAGGAATCAGAGAGGAGTGTTGACAGTGTAGCACGGGATG",
			"chr_3:536-542" => "GTATACG", 
			"chr.2:765-770" => "TCACAT"
    	}

    	test_regions.each_pair do |region, seq_expected|  
    		seq_fetched = @fasta_samtools.fetch_sequence(region)
    		assert_equal(seq_expected, seq_fetched.upcase)

    		seq_fetched_local = @fasta_local.fetch_sequence(region)
    		assert_equal(seq_expected, seq_fetched_local.upcase)
    	end

    	test_regions = {
    		"chr_3:772-780" => Bio::DB::Fasta::FastaDBException,
            "chr_3:-1-10" => Bio::DB::Fasta::FastaDBException,
            "chr_3:0-10" => Bio::DB::Fasta::FastaDBException,
            "bfafdaads" => Bio::DB::Fasta::FastaDBException,
            "chr.2:765-771" => Bio::DB::Fasta::FastaDBException,
            "chr_2:765-770" => Bio::DB::Fasta::FastaDBException,

    	}
    	test_regions.each_pair do |region, expected_eception|  
            assert_raise expected_eception do
                seq_fetched = @fasta_samtools.fetch_sequence(region)
            end

            assert_raise expected_eception do
    	       seq_fetched_local = @fasta_local.fetch_sequence(region)
    	   end
        end
    end

end