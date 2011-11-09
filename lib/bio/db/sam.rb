require 'bio/db/sam/library'
require 'bio/db/sam/bam'
require 'bio/db/sam/faidx'
require 'bio/db/sam/sam'
require 'bio/db/sam/pileup'
require 'bio/db/sam/vcf'

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

      # To make a new sam object. Initialize expects a hash optsa with the following elemets:
      # fasta:: The fasta file with the reference. (nil)
      # bam:: path to a binary SAM file (nil)
      # tam:: path to a text SAM file (nil) 
      # compressed:: If the binary file is compressed (true)
      # write:: If the file is to be writen (false). Not supported yet. 
      # *NOTE:* you can't use binary and text formats simultaneusly. To make queries, the file has to be a sorted binary. 
      # This function doesn't actually open the file, it just prepares the object to be opened in a later stage. 
      #
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
      
      #Function that actually opens the sam file
      #Throws a SAMException if the file can't be open.
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

      #Prints a description of the sam file in a text format containg if it is binary or text, the path
      #and the fasta file of the reference
      def to_s()
        (@binary ? "Binary" : "Text") + " file: " + @sam + " with fasta: " + @fasta_path
      end

      #Closes the sam file and destroys the C pointers using the functions provided by libbam
      def close()
        Bio::DB::SAM::Tools.fai_destroy(@fasta_index) unless @fasta_index.nil? || @fasta_index.null?
        Bio::DB::SAM::Tools.bam_index_destroy(@sam_index) unless @sam_index.nil? || @sam_index.null?
        Bio::DB::SAM::Tools.samclose(@sam_file) unless @sam_file.nil? 
        @sam_file = nil
        @fasta_index = nil
      end

      # Destructor method that closes the file before letting the object be garbage collected.
      def Sam.finalize(id)
        id.close()
        puts "Finalizing #{id}  at #{Time.new}"       
      end

      #Loads the bam index to be used for fetching. If the index doesn't exists the index is built provided that
      #the user has writing access to the folder where the BAM file is located. If the creation of the file fails
      #a SAMException is thrown. 
      #If the index doesn't exist, loading it will take more time. It is suggested to generate the index separatedly
      #if the bam file sits on a server where the executing user may not have writing permissions in the server.
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

      #Loads the reference file to be able to query regions of it. This requires the fai index to exist in the same
      #folder than the reference. If it doesn't exisits, this functions attempts to generate it. If user doesn't
      #have writing permissions on the folder, or the creation of the fai fails for any reason, a SAMException is thrown.
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

      #Returns the average coverage of a region in a bam file.
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

      #Returns an array with the coverage at each possition in the queried region
      #This is a simple average coverage just calculated with the first and last
      #possition of the alignment, ignoring the gaps.
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

      #Returns the sequence for a given region.
      def fetch_reference(chromosome, qstart,qend)
        load_reference if @fasta_index.nil? || @fasta_index.null? 
        query = query_string(chromosome, qstart,qend)
        len = FFI::MemoryPointer.new :int
        reference = Bio::DB::SAM::Tools.fai_fetch(@fasta_index, query, len)
        raise SAMException.new(), "Unable to get sequence for reference: "+query if reference.nil?

        reference
      end

      #Generates a query sting to be used by the region parser in samtools. 
      #In principle, you shouldn't need to use this function.
      def query_string(chromosome, qstart,qend)
        query = chromosome + ":" + qstart.to_s + "-" + qend.to_s 
        query
      end

      #Returns an array of Alignments on a given region.
      def fetch(chromosome, qstart, qend)
        als = Array.new
        fetchAlignment = Proc.new do |alignment|
          als.push(alignment.clone)   
          0  
        end
        fetch_with_function(chromosome, qstart, qend, fetchAlignment)
        als
      end  

      #Executes a function on each Alignment inside the queried region of the chromosome. The chromosome
      #can be either the textual name or a FixNum with the internal index. However, you need to get the
      #internal index with the provided API, otherwise the pointer is outside the scope of the C library. 
      #Returns the count of alignments in the region. 
      #WARNING: Accepts an index already parsed by the library. It fails when you use your own FixNum (FFI-bug?)
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
      
      #Merges n BAM files. This doesn't require to create a SAM object
      #files:: An array with the paths to the files.
      #merged_file:: The path to the merged file
      #headers:: The BAM file containing the header
      #add_RG:: If true, the RG tag is added (infered from the filenames)
      #by_qname:: If true, the bamfiles should by ordered by query name, if false, by coordinates. 
      def self.merge(files, merged_file, headers, add_RG, by_qname)
        strptrs = []
        strptrs << FFI::MemoryPointer.from_string("merge")
        files.each do |file|
          strptrs << FFI::MemoryPointer.from_string(file)
        end
        strptrs << nil

        # Now load all the pointers into a native memory block
        argv = FFI::MemoryPointer.new(:pointer, strptrs.length)
        strptrs.each_with_index do |p, i|
           argv[i].put_pointer(0,  p)
        end
        #void bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int add_RG)
        Bio::DB::SAM::Tools.bam_merge_core(by_qname, merged_file, headers, strptrs.length, argv, add_RG)
      end

      #calls the mpileup function, opts is a hash of options identical to the command line options for mpileup.
      #is an iterator that yields a Pileup object for each postion
      #the command line options that generate/affect BCF/VCF are ignored ie (g,u,e,h,I,L,o,p)
      #call the option as a symbol of the flag, eg -r for region is called :r => "some SAM compatible region"
      #eg bam.mpileup(:r => "chr1:1000-2000", :q => 50) gets the bases with quality > 50 on chr1 between 1000-5000 
      def mpileup( opts )

              raise SAMException.new(), "No BAMFile provided" unless @sam and @binary
              raise SAMException.new(), "No FastA provided" unless @fasta_path
              #long option form to short samtools form..
              long_opts = {
              :region => :r,
              :illumina_quals => :six,
              :count_anomalous => :A,
              :no_baq => :B,
              :adjust_mapq => :C,
              :max_per_bam_depth => :d,
              :extended_baq => :E,
              :exclude_reads_file => :G,
              :list_of_positions => :l,
              :mapping_quality_cap => :M,
              :ignore_rg => :R,
              :min_mapping_quality => :q,
              :min_base_quality => :Q
              }

              ##convert any long_opts to short opts 
              opts.each_pair do |k,v|
                if long_opts[k]
                  opts[long_opts[k]] = v 
                  opts.delete(k)
                end
              end

              ##remove any calls to -g or -u for mpileup, bcf output is not yet supported
              ##and also associated output options
              [:g, :u, :e, :h, :I, :L, :o, :p].each {|x| opts.delete(x) }

              strptrs = []
              strptrs << FFI::MemoryPointer.from_string("mpileup")
              opts.each do |k,v|
                next unless opts[k] ##dont bother unless the values provided are true.. 
                k = '6' if k == :six
                k = '-' + k.to_s
                strptrs << FFI::MemoryPointer.from_string(k)
                strptrs << FFI::MemoryPointer.from_string(v.to_s) unless ["-R", "-B", "-E", "-6", "-A"].include?(k) #these are just flags so don't pass a value...
              end
              strptrs << FFI::MemoryPointer.from_string('-f')
              strptrs << FFI::MemoryPointer.from_string(@fasta_path)
              strptrs << FFI::MemoryPointer.from_string(@sam)
              strptrs << nil

              # Now load all the pointers into a native memory block
              argv = FFI::MemoryPointer.new(:pointer, strptrs.length)
              strptrs.each_with_index do |p, i|
                 argv[i].put_pointer(0,  p)
              end

              old_stdout = STDOUT.clone
              read_pipe, write_pipe = IO.pipe()
              STDOUT.reopen(write_pipe)
                #int bam_mpileup(int argc, char *argv[])
                Bio::DB::SAM::Tools.bam_mpileup(strptrs.length - 1,argv)
                if fork
                  write_pipe.close
                  STDOUT.reopen(old_stdout) #beware .. stdout from other processes eg tests calling this method can get mixed in...
                  begin
                    while line = read_pipe.readline
                      yield Pileup.new(line)
                    end
                  rescue EOFError
                    read_pipe.close
                    Process.wait
                  end
                end
            end
      
      # utility method that does not use the samtools API, it calls samtools directly as if on the command line and catches the output,
      # to use this method you must have a version of samtools that supports the pileup command (< 0.1.17)
      # otherwise the command will fail.
      # mpileup is the preferred method for getting pileups. 
      # With this method the sam object should be created as usual, but you need to pass this method a string of options for samtools
      # you don't need to provide the call to samtools pileup itself or -f <fasta file> or the bam file itself, these are taken from the sam object
      def deprecated_pileup( cmd )
        
        system('samtools pileup > /dev/null 2>&1')
        ##assumes samtools is in the path...
        if $?.exitstatus > 1
          raise RuntimeError, "samtools is required on the path. A version of samtools with the pileup function is required"
        end
        
        raise SAMException.new(), "No BAMFile provided" unless @sam and @binary
        raise SAMException.new(), "No FastA provided" unless @fasta_path
        
        command = 'samtools pileup ' + cmd + " -f #{@fasta_path}" + " #{@sam}" 
        
        pipe = IO.popen(command)
        while line = pipe.gets
          yield Pileup.new(line)
        end
        pipe.close
      end
      
      
      def index_stats 
        raise SAMException.new(), "No BAMFile provided" unless @sam and @binary
        raise SAMException.new(), "No FastA provided" unless @fasta_path
        strptrs = []
        strptrs << FFI::MemoryPointer.from_string("idxstats")
        strptrs << FFI::MemoryPointer.from_string(@sam)
        strptrs << nil

        # Now load all the pointers into a native memory block
        argv = FFI::MemoryPointer.new(:pointer, strptrs.length)
        strptrs.each_with_index do |p, i|
           argv[i].put_pointer(0,  p)
        end
        
        index_stats = {}
        
        old_stdout = STDOUT.clone
        read_pipe, write_pipe = IO.pipe()
        STDOUT.reopen(write_pipe)
        
        #int bam_idxstats(int argc, char *argv[])
        Bio::DB::SAM::Tools.bam_idxstats(strptrs.length - 1,argv)
        if fork
          write_pipe.close
          STDOUT.reopen(old_stdout) #beware .. stdout from other processes eg tests calling this method can get mixed in...
          begin
            
            while line = read_pipe.readline #TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
                info = line.split(/\t/)
                next unless info.length == 4
                index_stats[ info[0] ] = {:length => info[1].to_i, :mapped_reads => info[2].to_i, :unmapped_reads => info[3].to_i } 
            end
            rescue EOFError
              read_pipe.close
              Process.wait
            end
        end #fork
        index_stats
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


