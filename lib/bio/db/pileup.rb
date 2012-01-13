# :title:Pileup
# = Bio::DB::Pileup 
# A class representing information in SAMTools pileup format
# Author:: Dan MacLean (dan.maclean@tsl.ac.uk)
# Pileup is described at http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ#I_do_not_understand_the_columns_in_the_pileup_output.
# Briefly (when you invoke pileup with the -c option):
# * 1 reference sequence name
# * 2 reference coordinate
# * (3) reference base, or `*' for an indel line
# * (4) genotype where heterozygotes are encoded in the IUB code: M=A/C, R=A/G, W=A/T, S=C/G, Y=C/T and K=G/T; indels are indicated by, for example, */+A, -A/* or +CC/-C. There is no difference between */+A or +A/*.
# * (5) Phred-scaled likelihood that the genotype is wrong, which is also called `consensus quality'.
# * (6) Phred-scaled likelihood that the genotype is identical to the reference, which is also called `SNP quality'. Suppose the reference base is A and in alignment we see 17 G and 3 A. We will get a low consensus quality because it is difficult to distinguish an A/G heterozygote from a G/G homozygote. We will get a high SNP quality, though, because the evidence of a SNP is very strong.
# * (7) root mean square (RMS) mapping quality
# * 8 # reads covering the position
# * 9 read bases at a SNP line (check the manual page for more information); the 1st indel allele otherwise
# * 10 base quality at a SNP line; the 2nd indel allele otherwise
# * (11) indel line only: # reads directly supporting the 1st indel allele
# * (12) indel line only: # reads directly supporting the 2nd indel allele
# * (13) indel line only: # reads supporting a third indel allele
# If pileup is invoked without `-c', indel lines and columns between 3 and 7 inclusive will not be outputted.
# 
# NB mpileup uses the 6 column output format eg
# "seq2\t151\tG\tG\t36\t0\t99\t12\t...........A\t:9<;;7=<<<<<"
# Pileup provides accessors for all columns (6 or 10 column format) and a few other useful methods
# 
# 
module Bio
  class DB
class Pileup
  attr_accessor :ref_name, :pos, :ref_base, :coverage, :read_bases, :read_quals, :consensus, :consensus_quality, :snp_quality, :rms_mapq, :ar1, :ar2, :ar3
  
  #creates the Pileup object
  #    pile_up_line = "seq2\t151\tG\tG\t36\t0\t99\t12\t...........A\t:9<;;7=<<<<<"
  #    pile = Bio::DB::Pileup.new(pile_up_line)
  def initialize(pile_up_line)
    cols = pile_up_line.split(/\t/)
    if cols.length == 6 ##should only be able to get 6 lines from mpileup
      @ref_name, @pos, @ref_base, @coverage, @read_bases, @read_quals = cols
    elsif (10..13).include?(cols.length) ##incase anyone tries to use deprecated pileup with -c flag we get upto 13 cols...
      @ref_name, @pos, @ref_base, @consensus, @consensus_quality, @snp_quality, @rms_mapq, @coverage, @read_bases, @read_quals, @ar1, @ar2, @ar3 = cols
      @consensus_quality = @consensus_quality.to_f
      @snp_quality = @snp_quality.to_f
      @rms_mapq = @rms_mapq.to_f
    else
      #raise RuntimeError, "parsing line '#{pile_up_line.chomp}' failed"
    end
      
    @pos = @pos.to_i
    @coverage = @coverage.to_f
    @ref_count = nil
    @non_ref_count_hash = nil 
    @non_ref_count = nil
  end
  
  # Calculate the total count of each non-reference nucleotide and return a hash of all 4 nt counts, returns a hash
  #    pile.non_refs #{:A => 1, :C => 0, :T => 0, :G => 0}
  def non_refs
    if @non_ref_count_hash.nil?
       @non_ref_count_hash = {:A => self.read_bases.count("Aa"), :C => self.read_bases.count("Cc"), :G => self.read_bases.count("Gg"), :T => self.read_bases.count("Tt")}
    end
      @non_ref_count_hash
  end
  
  # returns the total non-reference bases in the reads at this position
  def non_ref_count
    if @non_ref_count.nil?
      @non_ref_count = @read_bases.count("ATGCatgc").to_f
    end
    @non_ref_count
  end
  
  # returns the count of reference-bases in the reads at this position
  def ref_count
    if @ref_count.nil?
      @ref_count = self.read_bases.count(".,")
    end
    @ref_count
  end
  
  # returns the consensus (most frequent) base from the pileup, if there are equally represented bases returns a string containing all equally represented bases in alphabetical order   
  def consensus
      if @consensus.nil?
        max = self.non_refs.values.max
        if (self.ref_count / self.coverage) > 0.5
          @consensus = self.ref_base 
        elsif self.ref_count > max
          @consensus = self.ref_base
        else
          arr = self.non_refs.select {|k,v| v == max }
          bases = arr.collect {|b| b[0].to_s }
          bases << self.ref_base if self.ref_count == max
          @consensus = bases.sort.join
        end
      end
      @consensus
  end
  
end
end
end
