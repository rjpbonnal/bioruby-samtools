module Bio
  class DB
    class Sam
      attr_accessor :sams, :bams, :fasta
      
      # Creates a new Bio::DB::Sam object
      # @param fasta [String] the path to the Fasta reference sequence
      # @param bam [String,Array] path or list of paths to bam files
      # @param sam [String,Array] path or list of paths to bam files
      # @param samtools [String] path to alternative installation of samtools
      # @return [Bio::DB::Sam] a new `Bio::DB::Sam` object
      def initialize(args)
        @fasta = args[:fasta]
        @sams = args[:sam]
        @bams = args[:bam]
        @samtools = args[:samtools]
        
        @sams = [@sams] if @sams.instance_of?(String)
        @bams = [@bam] if @bams.instance_of?(String)
        
        raise ArgumentError, "Need Fasta and at least one BAM or SAM" if not @fasta and (@bam or @sam)
        raise IOError, "File not found" if not files_ok?
        
      end
      
      #backward compatibility method, returns true if file exists otherwise, complains and quits.
      def open
        files_ok?
      end
      
      private
      # checks existence of file
      def files_ok?
        [@fasta, @sam, @bam].flatten.compact.each {|f| return false unless File.exists? f }
        true
      end
      
      
    end
  end
end
