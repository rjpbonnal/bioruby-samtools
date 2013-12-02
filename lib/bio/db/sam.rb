module Bio
  class DB
    class Sam
      attr_accessor :bam, :fasta, :samtools, :bcftools, :last_command
      
      
      # Creates a new Bio::DB::Sam object
      # @param fasta [String] the path to the Fasta reference sequence
      # @param bam [String] path to bam files
      # @param samtools [String] path to alternative installation of samtools
      # @param bcftools [String] path to alternative installation of bcftools
      # @return [Bio::DB::Sam] a new `Bio::DB::Sam` object
      def initialize(args)
        @fasta = args[:fasta]
        @bam = args[:bam]
        @samtools = args[:samtools] || File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','samtools')
        @bcftools = args[:bcftools] || File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','bcftools')
        
        @files = [@files] if @files.instance_of?(String)
        @last_command = nil
        raise ArgumentError, "Need Fasta and at least one BAM or SAM" if not @fasta or not @bam
        raise IOError, "File not found" if not files_ok?
        @bams = [@bams] if @bams.instance_of? String
        
      end

      #backward compatibility method, returns true if file exists otherwise, complains and quits.
      def open
        files_ok?
      end
      
      # runs the samtools view command
      # @param opts [Hash] options for view as follows
      # :b => nil, # -b output BAM
      # :h => nil, # -h print header for the SAM output
      # :H => nil, # -H print header only (no alignments)
      # :S => nil, # -S input is SAM
      # :u => nil, # -u uncompressed BAM output (force -b)
      # :one => nil, # -1 fast compression (force -b)
      # :x => nil, # -x output FLAG in HEX (samtools-C specific)
      # :X => nil, # -X output FLAG in string (samtools-C specific)
      # :c => nil, # -c print only the count of matching records
      # :B => nil, # -B collapse the backward CIGAR operation
      # :at => nil, # -@ INT number of BAM compression threads [0]
      # :L => nil, # -L FILE output alignments overlapping the input BED FILE [null]
      # :t => nil, # -t FILE list of reference names and lengths (force -S) [null]
      # :T => nil, # -T FILE reference sequence file (force -S) [null]
      # :o => nil, # -o FILE output file name [stdout]
      # :R => nil, # -R FILE list of read groups to be outputted [null]
      # :f => nil, # -f INT required flag  0 for unset [0]
      # :F => nil, # -F INT filtering flag  0 for unset [0]
      # :q => nil, # -q INT minimum mapping quality [0]
      # :l => nil, # -l STR only output reads in library STR [null]
      # :r => nil, # -r STR only output reads in read group STR [null]
      # :s => nil # -s FLOAT fraction of templates to subsample; integer part as seed [-1]
      # :chr => nil # name of reference sequence to get alignments from
      # :start => nil # start position on reference sequence
      # :stop => nil # end postion on reference sequence
      def view(opts={},&block)
        region = String.new
        if opts[:chr] and opts[:start] and opts[:stop]
          region = "#{opts[:chr]}:#{opts[:start]}-#{opts[:stop]}"
          [:chr, :start, :stop].each {|o| opts.delete(o)}
        end
        if opts[:at]
          opts['@'] = opts[:at]
          opts.delete(:at)
        end
        
        if opts[:one]
          opts['1'] = opts[:one]
          opts.delete(:one)
        end
                
        command = form_opt_string(@samtools, 'view', opts, [:b, :h, :H, :S, :u, '1', :x, :X, :c, :B]) + " " + region
        @last_command = command
        type = (opts[:u] or opts[:b]) ? :binary : :text
        klass = (type == :binary) ? String : Bio::DB::Alignment
        yield_from_pipe(command, klass, type, &block)
      end
      
      def fetch(chr, start,stop, &block)
        view(
        :chr => chr,
        :start => start,
        :stop => stop, 
        &block  
        )
      end
      
      alias_method :fetch_with_function, :fetch
      
      def chromosome_coverage(chr,start,length)
        result = []
        region = "#{chr}:#{start}-#{start + length}"
        self.mpileup(:r => region) do |p|
          result << p.coverage
        end
        result
      end
      
      def average_coverage(chr,start,length)
        arr = self.chromosome_coverage(chr,start,length)
        arr.inject{ |sum, el| sum + el }.to_f / arr.size
      end
      
      def mpileup(opts={}, &block)
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
        :min_base_quality => :Q,
        ###following options are for the -g -u option
        :genotype_calling => :g,
        :uncompressed_bcf => :u,
        :extension_sequencing_probability => :e,
        :homopolymer_error_coefficient => :h,
        :no_indels => :I,
        :skip_indel_over_average_depth => :L,
        :gap_open_sequencing_error_probability => :o,
        :platforms => :P 
        }

        ##convert any long_opts to short opts 
        temp_opts = opts.dup
        opts.each_pair do |k,v|
          if long_opts[k]
            temp_opts[long_opts[k]] = v 
            temp_opts.delete(k)
          end
        end
        opts = temp_opts
        opts[:u] = true if opts[:g] #so that we always get uncompressed output
        opts.delete(:g)
        
        opts[:f] = @fasta
        
        if opts[:six]
          opts["6"] = nil
          opts.delete(:six)
        end
        
        command = form_opt_string(@samtools, "mpileup", opts, [:R, :B, :E, "6", :A, :g, :u, :I] )

        if opts[:u]
          command = command + " | #{@bcftools} view -cg -"
        end
        
        klass = opts[:u] ? Bio::DB::Vcf : Bio::DB::Pileup
        @last_command = command
        yield_from_pipe(command, klass, :text, &block)

      end
      
      def fetch_reference(chr,start,stop, opts={:as_bio => false})
        command = "#{@samtools} faidx #{@fasta} #{chr}:#{start}-#{stop}"
        @last_command = command
        seq = ""
        yield_from_pipe(command, String, :text ) {|line| seq = seq + line unless line =~ /^>/}
        if opts[:as_bio]
          seq = Bio::Sequence::NA.new(seq).to_fasta("#{chr}:#{start}-#{stop}")
        end
        seq
      end
      
  
    def faidx(opts={})
      if opts.has_key?(:chr) and opts.has_key?(:start) and opts.has_key?(:stop)
      opts={:as_bio => false}
      self.fetch_reference(:chr,:start,:stop,opts)
      else
        command = "#{@samtools} faidx #{@fasta}"
        @last_command = command
        system(command)
      end
    end

      #:out_index name of index
      def index(opts={})
        opts.merge!({:out_index=>nil})
        command = form_opt_string(@samtools, "index",opts) + " #{opts[:out_index]}"
        @last_command = command
        system(command)
      end

      #:out_bam name of outfile
      #:r  remove unmapped reads and secondary alignments
      def fix_mates(opts={})
        opts.merge!({:out_index=>nil})
        command = "#{form_opt_string(@samtools, "fixmate", opts, [:r])} #{opts[:out_bam]}"
        @last_command = command
        system(command)
      end

      alias_method :fixmate, :fix_mates

      def flag_stats(opts={})
        command = form_opt_string(@samtools, "flagstat", opts, [])
        @last_command = command
        strings = []
        yield_from_pipe(command,String) {|line| strings << line.chomp}
        strings
      end

      alias_method :flagstat, :flag_stats

      def index_stats
        stats = {}   
        command = form_opt_string(@samtools, "idxstats #{@fasta}", {}, [])
        @last_command = command
        yield_from_pipe(command, String, :text,skip_comments=true, comment_char="#") do |line|
          info = line.chomp.split(/\t/)
          stats[ info[0] ] = {:length => info[1].to_i, :mapped_reads => info[2].to_i, :unmapped_reads => info[3].to_i }
        end
        stats
      end

      alias_method :idxstats, :index_stats


      #:n       sort by read names
      #:r       attach RG tag (inferred from file names)
      #:u       uncompressed BAM output
      #:f       overwrite the output BAM if exist
      #:one       compress level 1
      #:l INT   compression level, from 0 to 9 [-1]
      #:at INT   number of BAM compression threads [0]
      #:R STR   merge file in the specified region STR [all]
      #:h FILE  copy the header in FILE to <out.bam> [in1.bam]
      #:out FILE     out file name
      #:bams FILES or Bio::DB::Sam list of input bams, or Bio::DB::Sam objects
      def self.merge(opts={})
        if opts[:one]
          opts['1'] = nil
          opts.delete(:one)
        end

        if opts[:at]
          opts['@'] = opts[:at]
          opts.delete(:at)
        end

        out = opts[:out]
        opts.delete(:out)

        bam_list = opts[:bams].collect do |b|
          b.bam rescue b
        end.join(' ')

        opts.delete(:bams)

        command = "#{form_opt_string(@samtools, "merge", opts, [:n, :r, :u, :f, '1'] )} #{out} #{bam_list}"
        @last_command = command
        system(command)

      end

      #:h  header.sam
      #:out FILE     out file name
      #:bams FILES or Bio::DB::Sam list of input bams, or Bio::DB::Sam objects
      def self.cat(opts)

        out = opts[:out]
        opts.delete(:out)

        bam_list = opts[:bams].collect do |b|
          b.bam rescue b
        end.join(' ')

        command = "#{form_opt_string(@samtools, "cat", opts, [:h])} -o #{out} #{bam_list}"
        @last_command = command
        system(command)

      end

      #program  one of 'samtools' 'bcftools'
      #command one of the commands relevant to the program
      def self.docs(program, command)
        return "program must be 'samtools' or 'bcftools' " if not ['samtools', 'bcftools'].include? program
        command = "#{program} #{command}"
        `#{command}`
      end

      #:s   rmdup for SE reads
      #:S   treat PE reads as SE in rmdup (force -s)
      #:out FILE output bam
      def remove_duplicates(opts={})
        out = opts[:out]
        opts.delete(:out)
        command = "#{form_opt_string(@samtools, "rmdup", opts, [:s, :S])} #{out} #{@bam}"
        @last_command = command
        system(command)
      end

      alias_method :rmdup, :remove_duplicates
      
      #        :n  sort by read name
      #        :f        use <out.prefix> as full file name instead of prefix
      #        :o        final output to stdout returns bio::db::alignment
      #        :l INT    compression level, from 0 to 9 [-1]
      #        :at INT    number of sorting and compression threads [1]
      #        :m INT    max memory per thread; suffix K/M/G recognized [768M]
      #        :prefix prefix for output bamfile
      def sort(opts={})
        opts.merge!({:prefix => "sorted"})
        prefix = opts[:prefix]
        opts.delete(:prefix)
        command = form_opt_string(@samtools, "sort", opts, [:n, :f, :o])
        command = command + " " + prefix
        @last_command = command
        if opts[:o]
          yield_from_pipe(command, Bio::DB::Alignment)
        else
          system(command)
        end
      end
      
      
      def tview(opts={})
        #to do
      end
       
      def reheader(header_sam)
        command = "#{@samtools} reheader #{header_sam} #{@bam}"
        @last_command = command
        system(command)
      end
      
      #	-A 	When used jointly with -r this option overwrites the original base quality.
      #	:e 	Convert a the read base to = if it is identical to the aligned reference base. Indel caller does not support the = bases at the moment.
      #	:u 	Output uncompressed BAM
      #	:b 	Output compressed BAM
      #	:S 	The input is SAM with header lines
      #	:C INT 	Coefficient to cap mapping quality of poorly mapped reads. See the pileup command for details. [0]
      #	:r 	Compute the BQ tag (without -A) or cap base quality by BAQ (with -A).
      #	:E 	Extended BAQ calculation. This option trades specificity for sensitivity, though the effect is minor. 
      def calmd(opts={})
        command = "#{form_opt_string(@samtools, "calmd", opts, [:E, :e, :u, :b, :S, :r] )} #{bam}"
        @last_command = command
        system(command)
      end
      
      def targetcut(opts={})
        #to do
      end
      
      def phase(opts={})
        #to do
      end


      # :b <bed>            list of positions or regions
      # :l <int>            minQLen
      # :q <int>            base quality threshold
      # :Q <int>            mapping quality threshold
      # :r <chr:from-to>    region
      #returns an array for each position with [sequence_name, position, depth]
      def depth(opts={})
        command = form_opt_string(@samtools, "depth", opts)
        @last_command = command
        puts command
        yield_from_pipe(command, String) do |line|
          yield line.split(/\t/)
        end

      end

      private
      
      # returns a command string from a program
      # @param program [Symbol] either `:samtools` or `:bcftools`
      # @param opts [Hash] the options hash
      # @param singles `flag` options [Array] the options in `opts` that are single options 
      def form_opt_string(prog, command, opts, singles=[])
        opts_string = commandify(opts, singles)
        "#{prog} #{command} #{opts_string} #{@bam}"
      end
      
      # turns an opts hash into a s
      def commandify(opts, singles)
        list = []
        opts.each_pair do |tag,value|
          value = "" if singles.include?(tag)
          list << "-#{tag.to_s} #{value}" 
        end
        list.join(" ")
      end
      
      # checks existence of files in instance
      def files_ok?
        [@fasta, @sam, @bam].flatten.compact.each {|f| return false unless File.exists? f }
        true
      end
      
      def yield_from_pipe(command, klass, type=:text, skip_comments=true, comment_char="#", &block)
        pipe = IO.popen(command)
        if type == :text
          while (line = pipe.gets)
            next if skip_comments and line[0] == comment_char
            yield klass.new(line.chomp)
          end
        elsif type == :binary
          while (c = pipe.gets(nil))
            yield c
          end
        end
        pipe.close
      end
      
    end
  end
end
