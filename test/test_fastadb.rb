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

	def test_faidx
		@fasta_samtools.faidx()
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

    	# { :dog => "charlie", :cat => "kwiki", :mouse => "squeaky" }
    	test_regions = { "chr_1:1-70" => "CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA", 
    		"chr_1:68-74" => "CTAACCC",
    		"chr.2:300-450" => "GGTGCAGTCACGGCTCGCAGTTATACTCCGGGAGAATGGAAAAGATTGTCCAAGGACCAA\
CAGGAAAAAGTGCGATCGCTGCGTAATAAAAAGAAGCAAGGAGGGAAACCCGAGGAATCAGAGAGGAGTGTTGACAGTGTAGCACGGGATG"
    	}

    	test_regions.each_pair do |region, seq_expected|  
    		seq_fetched = @fasta_samtools.fetch_sequence(region)
    		assert_equal(seq_expected, seq_fetched.upcase)

    		seq_fetched_local = @fasta_local.fetch_sequence(region)
    		assert_equal(seq_expected, seq_fetched_local.upcase)
    	end
    end

end