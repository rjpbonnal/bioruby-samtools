require 'bio/db/sam/library'
require 'bio/db/sam/bam'
require 'bio/db/sam/faidx'
require 'bio/db/sam/sam'

module LibC
  extend FFI::Library
  ffi_lib FFI::Library::LIBC
  attach_function :free, [ :pointer ], :void
  # call #attach_function to attach to malloc, free, memcpy, bcopy, etc.
end

module Bio
  class DB
    class Sam
      attr_reader :sam_file

      def initialize(optsa={})
        opts =  { :fasta => nil,  :bam => nil,:tam => nil, :compressed => true, :write => false }.merge!(optsa)



        @fasta_path = opts[:fasta]
        @compressed = opts[:compressed]
        @write = opts[:write]
        bam = opts[:bam]
        tam = opts[:tam]

        if bam == nil && tam == nil && @fasta_path == nil then
          raise SAMException.new(), "No alignment or reference file"
        elsif bam != nil && tam != nil then
          raise SAMException.new(), "Alignment has to be in either text or binary format, not both"
        elsif bam != nil then
          @binary = true
          @sam = bam     
        elsif tam != nil then
          @sam = tam     
          @binary = false

        end
        @fasta_file = nil
        @sam_file   = nil

        ObjectSpace.define_finalizer(self,  self.class.method(:finalize).to_proc)
      end

      def open()

        raise SAMException.new(), "Writing not supported yet" if @write
        raise SAMException.new(), "No SAM file specified" unless @sam 

        opts = @write ? "w" : "r"
        if @binary then  
          opts += "b" 
          if @write then
            unless @compressed then 
              opts += "u"
            end
          end
        end
        valid = ["r", "w", "wh", "rb", "wb" , "wbu"]
        unless valid.include?(opts) then
          raise SAMException.new(), "Invalid options for samopen: " + opts 
        end

        samFile = Bio::DB::SAM::Tools.samopen(@sam, opts, nil)
        if samFile.null? then
          @sam_file = nil
          raise SAMException.new(), "File not opened:  " + @sam
        end
        @sam_file = Bio::DB::SAM::Tools::SamfileT.new(samFile)

      end

      def to_s()
        (@binary ? "Binary" : "Text") + " file: " + @sam + " with fasta: " + @fasta_path
      end

      def close()
        Bio::DB::SAM::Tools.fai_destroy(@fasta_index) unless @fasta_index.nil? || @fasta_index.null?
        Bio::DB::SAM::Tools.bam_index_destroy(@sam_index) unless @sam_index.nil? || @sam_index.null?
        Bio::DB::SAM::Tools.samclose(@sam_file) unless @sam_file.nil? 
        @sam_file = nil
        @fasta_index = nil
      end

      def Sam.finalize(id)
        id.close()
        puts "Finalizing #{id}  at #{Time.new}"       
      end

      def load_index()
        raise SAMException.new(), "Indexes are only supported by BAM files, please use samtools to convert your SAM file" unless @binary
        @sam_index = Bio::DB::SAM::Tools.bam_index_load(@sam)
        if @sam_index.null? then
          p "Generating index for: " + @sam
          Bio::DB::SAM::Tools.bam_index_build(@sam)
          @sam_index = Bio::DB::SAM::Tools.bam_index_load(@sam)
          raise SAMException.new(), "Unable to generate bam index for: " + @sam if @sam_index.nil? || @sam_index.null?
        end
      end

      def load_reference()
        raise SAMException.new(), "No path for the refernce fasta file. " if @fasta_path.nil?

        @fasta_index = Bio::DB::SAM::Tools.fai_load(@fasta_path)

        if @fasta_index.null? then
          p "Generating index for: " + @fasta_path
          Bio::DB::SAM::Tools.fai_build(@fasta_path)
          @fasta_index =  Bio::DB::SAM::Tools.fai_load(@fasta_path)
          raise SAMException.new(), "Unable to generate fasta index for: " + @fasta_path if @fasta_index.nil? ||  @fasta_index.null?
        end

      end

      def average_coverage(chromosome, qstart, len)

        #reference = fetch_reference(chromosome, qstart,len)
        # len = reference.length if len > reference.length


        coverages = chromosome_coverage(chromosome, qstart, len)
        total = 0
        len.times{ |i| total= total + coverages[i] }
        avg_cov = total.to_f / len
        #LibC.free reference
        avg_cov
      end

      def chromosome_coverage(chromosome, qstart, len)
        #  reference = fetch_reference(chromosome, qstart,len)
        #  len = reference.length if len > reference.length
        #p qend.to_s + "-" + qstart.to_s + "framesize " + (qend - qstart).to_s
        coverages = Array.new(len, 0)

        chr_cov_proc = Proc.new do |alignment|
          #last = qstart + len
          #first = qstart
          #last = alignment.calend if last > alignment.calend
          #first = alignment.pos if first < alignment.pos
          # p first
          last = alignment.calend - qstart
          first = alignment.pos - qstart
          if last < first
            tmp = last
            last = first
            first = last
          end

          #  STDERR.puts "#{first} #{last}\n"
          first.upto(last-1) { |i|

            coverages[i-1] = 1 + coverages[i-1]  if i-1 < len && i > 0
          }
        end

        fetch_with_function(chromosome, qstart, qstart+len,  chr_cov_proc)
        #p coverages
        coverages
      end

      def fetch_reference(chromosome, qstart,qend)
        load_reference if @fasta_index.nil? || @fasta_index.null? 
        query = query_string(chromosome, qstart,qend)
        len = FFI::MemoryPointer.new :int
        reference = Bio::DB::SAM::Tools.fai_fetch(@fasta_index, query, len)
        raise SAMException.new(), "Unable to get sequence for reference: "+query if reference.nil?

        reference
      end

      def query_string(chromosome, qstart,qend)
        query = chromosome + ":" + qstart.to_s + "-" + qend.to_s 
        query
      end

      def fetch(chromosome, qstart, qend)
        als = Array.new
        fetchAlignment = Proc.new do |alignment|
          als.push(alignment.clone)   
          0  
        end
        fetch_with_function(chromosome, qstart, qend, fetchAlignment)
        als
      end  

      def fetch_with_function(chromosome, qstart, qend, function)
        load_index if @sam_index.nil? || @sam_index.null?
        chr = FFI::MemoryPointer.new :int
        beg = FFI::MemoryPointer.new :int
        last = FFI::MemoryPointer.new :int
        query = query_string(chromosome, qstart,qend)
        qpointer = FFI::MemoryPointer.from_string(query)
        header = @sam_file[:header]
        Bio::DB::SAM::Tools.bam_parse_region(header,qpointer, chr, beg, last) 
        #raise SAMException.new(), "invalid query: " + query  if(chr.read_int < 0)
        count = 0;

        fetchAlignment = Proc.new do |bam_alignment, data|
          alignment =  Alignment.new
          alignment.set(bam_alignment, header)
          function.call(alignment)
          count = count + 1
          0  
        end
        Bio::DB::SAM::Tools.bam_fetch(@sam_file[:x][:bam], @sam_index,chr.read_int,beg.read_int, last.read_int, nil, fetchAlignment)
        #LibC.free chr
        #LibC.free beg
        #LibC.free last
        #LibC.free qpointer
        count
      end

    end

    class Tag
      attr_accessor :tag, :type, :value
      def set(str)
        v = str.split(":")
        @tag   = v[0]
        @type  = v[1]
        @value = v[2]
      end
    end

    class Alignment

      def initialize
        ObjectSpace.define_finalizer(self,
        self.class.method(:finalize).to_proc)
      end
      def Alignment.finalize(object_id)

        #           puts "Object #{object_id} dying at #{Time.new}"
        #             p  "?" . object_id.al
        #            p object_id.al
        LibC.free object_id.al
        LibC.free object_id.sam
        LibC.free object_id.calend
        LibC.free object_id.qlen

        LibC.free object_id.samstr
      end

      #Attributes from the format
      attr_accessor :qname, :flag, :rname,:pos,:mapq,:cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags, :al, :samstr
      #Attributes pulled with the C library
      attr_accessor  :calend, :qlen
      #Attrobites frp, the flag field (see chapter 2.2.2 of the sam file documentation)
      #query_strand and mate_strand are true if they are forward. It is the opposite to the definition in the BAM format for clarity.
      #primary is the negation of is_negative from the BAM format
      attr_accessor :is_paired, :is_mapped, :query_unmapped, :mate_unmapped, :query_strand, :mate_strand, :first_in_pair,:second_in_pair, :primary, :failed_quality, :is_duplicate

      def set(bam_alignment, header)
        #Create the FFI object
        @al = Bio::DB::SAM::Tools::Bam1T.new(bam_alignment) 

        #set the raw data
        tmp_str =  Bio::DB::SAM::Tools.bam_format1(header,al)
        #self.sam =  tmp_str
        #ObjectSpace.define_finalizer(self, proc {|id| puts "Finalizer one on #{id}" })
        self.sam = String.new(tmp_str)
        #LibC.free tmp_str
        #Set values calculated by libbam
        core = al[:core]
        cigar = al[:data][core[:l_qname]]#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname)) 
        @calend = Bio::DB::SAM::Tools.bam_calend(core,cigar)
        @qlen = Bio::DB::SAM::Tools.bam_cigar2qlen(core,cigar)

        #process the flags
        @is_paired             = @flag & 0x0001 > 0
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


      def sam=(sam)
        #p sam
        s = sam.split("\t")
        self.qname = s[0]
        self.flag  = s[1].to_i
        self.rname = s[2]
        self.pos   = s[3].to_i
        self.mapq  = s[4].to_i
        self.cigar = s[5]
        self.mrnm  = s[6]
        self.mpos  = s[7].to_i
        self.isize = s[8].to_i
        self.seq   = s[9]
        self.qual =  s[10]
        self.tags = {} 
        11.upto(s.size-1) {|n| 
          t = Tag.new 
          t.set(s[n])
          tags[t.tag] = t
        }


        #<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> \
        #[<TAG>:<VTYPE>:<VALUE> [...]] 

      end

    end

    class SAMException < RuntimeError
      #we can add further variables to give information of the excpetion
      def initialize()

      end
    end
  end
end


