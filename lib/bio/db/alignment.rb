module Bio
  class DB
    
    #a class to represent the SAM OPT values, presented in SAM as TAG:VTYPE:VALUE 
    class Tag
      attr_accessor :tag, :type, :value
      def set(str)
        @tag   = str[0..1]
        @type  = str[3]
        @value = str[5..-1]
      end
    end
    
    #Attrobites frp, the flag field (see chapter 2.2.2 of the sam file documentation)
    #query_strand and mate_strand are true if they are forward. It is the opposite to
    #the definition in the BAM format for clarity.
    #primary is the negation of is_negative from the BAM format
    class Alignment
      attr_accessor :qname, :flag, :rname,:pos,:mapq,:cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags, :al, :samstr, :calend, :qlen

      attr_accessor :is_paired, :is_mapped, :query_unmapped, :mate_unmapped, :query_strand, :mate_strand, :first_in_pair,:second_in_pair, :primary, :failed_quality, :is_duplicate
      
      #parses the SAM string into its constituents and set its attributes
      def initialize(sam_string)
        s = sam_string.chomp.split("\t")
        @qname = s[0]
        @flag  = s[1].to_i
        @rname = s[2]
        @pos   = s[3].to_i
        @mapq  = s[4].to_i
        @cigar = s[5]
        @mrnm  = s[6]
        @mpos  = s[7].to_i
        @isize = s[8].to_i
        @seq   = s[9]
        @qual =  s[10]
        @tags = {} 
        11.upto(s.size-1) {|n| 
          t = Bio::DB::Tag.new 
          t.set(s[n])
          tags[t.tag] = t
        }    
      
        @is_paired  = (@flag & 0x0001) > 0
        @is_mapped             = @flag & 0x0002 > 0
        @query_unmapped        = @flag & 0x0004 > 0
        @mate_unmapped         = @flag & 0x0008 > 0
        @query_strand          = !(@flag & 0x0010 > 0)
        @mate_strand           = !(@flag & 0x0020 > 0)
        @first_in_pair         = @flag & 0x0040 > 0
        @second_in_pair        = @flag & 0x0080 > 0
        @primary               = !(@flag & 0x0100 > 0)
        @failed_quality        = @flag & 0x0200 > 0
        @is_duplicate          = @flag & 0x0400 > 0
      
      end
    end
  end
end
