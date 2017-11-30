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
      attr_accessor :ref_name, :pos, :ref_base, :coverage, :read_bases, :read_quals, :consensus_quality, :snp_quality, :rms_mapq, :ar1, :ar2, :ar3, :indel_1, :indel_2
      
      #creates the Pileup object
      #    pile_up_line = "seq2\t151\tG\tG\t36\t0\t99\t12\t...........A\t:9<;;7=<<<<<"
      #    pile = Bio::DB::Pileup.new(pile_up_line)
      def initialize(pile_up_line)
        cols = pile_up_line.split(/\t/)
        @consensus = nil
        @consensus_quality = nil
        @read_quals = nil
        @bases = nil
        @allele_frequency = nil
        @consensus_iuap = nil
        if cols.length == 6 ##should only be able to get 6 lines from mpileup
          @ref_name, @pos, @ref_base, @coverage, @read_bases, @read_quals = cols
        elsif (10..13).include?(cols.length) ##incase anyone tries to use deprecated pileup with -c flag we get upto 13 cols...
          if cols[2] == '*' #indel
            @ref_name, @pos, @ref_base, @consensus, @consensus_quality, @snp_quality, @rms_mapq, @coverage, @indel_1, @indel_2, @ar1, @ar2, @ar3 = cols
          else #snp / identity
            @ref_name, @pos, @ref_base, @consensus, @consensus_quality, @snp_quality, @rms_mapq, @coverage, @read_bases, @read_quals = cols
          end
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
      
      #Calculate the total count of each non-reference nucleotide and return a hash of all 4 nt counts
      #returns a hash pile.non_refs #{:A => 1, :C => 0, :T => 0, :G => 0}
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
      
      # returns the consensus (most frequent) base from the pileup, if there are equally represented bases returns a string of all equally represented bases in alphabetical order   
      def consensus
          if @consensus.nil?
            max = self.non_refs.values.max
            #if the ref base is in more than half the coverage..
            if (self.ref_count / self.coverage) > 0.5
              #..then the ref base is the concensus
              @consensus = self.ref_base
            ##not sure if the following will ever apply as the non_refs method also returns the ref base count, hence can never be over the max count  
            #elsif self.ref_count > max
            #  @consensus = self.ref_base
            else
              #get the base(s) and count(s) that has the max count
              arr = self.non_refs.select {|k,v| v == max }
              #just get the bases (remove the counts)
              bases = arr.collect {|b| b[0].to_s }
              #add the ref base if the ref base has a max count (commenting this out as it should already be in)
              #bases << self.ref_base if self.ref_count == max
              @consensus = bases.sort.join
            end
          end
          @consensus
      end
      
      #returns basic VCF string as per samtools/misc sam2vcf.pl except that it scrimps on the ref for indels, returning a '*' instead of the reference allele
      def to_vcf
        alt,g = self.genotype_list 
        alt = self.consensus.split(//).join(',') unless self.ref_base == '*'
        alt = '.' if alt == self.ref_base
        alt = alt.split(',')
        #if the reference base is in alt, remove it
        alt.delete(self.ref_base.to_s)
        alt = alt.join(',')
        [self.ref_name, self.pos, '.', self.ref_base, alt, self.snp_quality.to_i, "0", "DP=#{self.coverage.to_i}", "GT:GQ:DP", "#{g}:#{self.consensus_quality.to_i}:#{self.coverage.to_i}" ].join("\t")
      end
      
      private
      def Pileup.vcf_header
        %{##fileformat=VCFv3.3\n##INFO=DP,1,Integer,"Total Depth"\n##FORMAT=GT,1,String,"Genotype"\n##FORMAT=GQ,1,Integer,"Genotype Quality"\n##FORMAT=DP,1,Integer,"Read Depth"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDATA\n}
      end
      
      
 
      #returns the genotype of the indel
      def indel_gt
        return "undef" if self.consensus.instance_of?(Array)
        al1, al2 = self.consensus.split(/\//)
        if al1 == al2 && al1 == '*'   
          al1=self.indel_1 
          al2=self.indel_2
        end
        alt1 = parse_indel(al1)
        alt2 = parse_indel(al2)
        alt,gt = nil,nil
        
        return nil if !alt1 and !alt2
        if !alt1 
          alt = alt2
          gt = '0/1'
        elsif !alt2
          alt = alt1
          gt - '0/1'
        elsif alt1 == alt2
          alt = alt1
          gt = '1/1'
        else
          alt="#{alt1},#{alt2}"
          gt= '1/2'
        end
        return [alt, gt]
        
      end
      #returns the genotype of the snp
      def snp_gt
        return ['.','0/0'] if self.ref_base == self.consensus
        bases = Pileup.iupac_to_base(self.consensus)
        if bases[0] == self.ref_base
           return [bases[1],'0/1']
        elsif bases[1] == self.ref_base
           return [bases[0],'0/1']
        else
          return ["#{bases[0]},#{bases[1]}",'1/1']
        end 
      end
      
      #identifies the reference base and returns the indel or snp genotype as applicable
      public
      def genotype_list
        if self.ref_base == '*'
          return indel_gt
        else
         return snp_gt
        end
      end
      
      #returns the two bases for the corresponding iupac code
      public
      def Pileup.iupac_to_base(alt_base)
          case alt_base
                when 'K' then ['G','T']
                when 'M' then ['A','C']
                when 'S' then ['C','G']
                when 'R' then ['A','G']
                when 'W' then ['A','T']
                when 'Y' then ['C','T']
                else alt_base.split(//)
          end
      end
      
      #identifies if the indel is an insertion or a deletion  
      def parse_indel(alt)
        return "D#{$'.length}" if alt =~/^-/ 
        if alt=~/^\+/
          return "I#{$'}"
        elsif alt == '*' 
          return nil
        end
      end
      
      
      #returns pileup format line
      def to_s
        if @read_quals and !@consensus_quality #6col
          [@ref_name, @pos, @ref_base, @coverage.to_i, @read_bases, @read_quals].join("\t")    
        elsif @indel_1 #13 cols
          [@ref_name, @pos, @ref_base, @consensus, @consensus_quality.to_i, @snp_quality.to_i, @rms_mapq.to_i, @coverage.to_i, @indel_1, @indel_2, @ar1, @ar2, @ar3].join("\t")
        else #10 cols
          [@ref_name, @pos, @ref_base, @consensus, @consensus_quality.to_i, @snp_quality.to_i, @rms_mapq.to_i, @coverage.to_i, @read_bases, @read_quals].join("\t")
        end
          
      end
      
      
      def bases
         return @bases if @bases
         @bases = self.non_refs
         #puts self.ref_count
         @bases[self.ref_base.upcase.to_sym] = self.ref_count 
         @bases
      end

       def base_coverage
         total = 0
         @bases.each do |k,v|
           total += v  
         end
         total
       end

      #returns the frequency of all bases in pileup position
       def allele_freq
         return @allele_frequency if @allele_frequency
         bases = self.bases
         @allele_frequency = Hash.new
         bases.each do |k,v| 
           @allele_frequency[k] = v.to_f/self.base_coverage.to_f 
         end
         @allele_frequency
       end

       # returns the consensus (most frequent) base from the pileup, if there are equally represented bases returns a string of all equally represented bases in alphabetical order   
       def consensus_iuap(minumum_ratio_for_iup_consensus)
         
         tmp = []
         if @consensus_iuap.nil?
           @consensus_iuap = self.ref_base.downcase
           bases = self.bases
           #tmp = String.new
           bases.each do |k,v|
             tmp << k[0].to_s if v/self.coverage.to_f > minumum_ratio_for_iup_consensus
           end
           if tmp.length > 0
             tmp = tmp.collect{ |x| Bio::Sequence::NA.new(x) }
              # creates alignment object
              a = Bio::Alignment.new(tmp)
              # shows IUPAC consensus
             @consensus_iuap = a.consensus_iupac
           end
         end 
         @consensus_iuap
       end
    end
  end
end
