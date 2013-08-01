module Bio
  class DB
    class Sam
      attr_accessor :bams, :fasta, :samtools, :bcftools, :last_command
      
      # Creates a new Bio::DB::Sam object
      # @param fasta [String] the path to the Fasta reference sequence
      # @param bam [String,Array] path or list of paths to bam files
      # @param sam [String,Array] path or list of paths to bam files
      # @param samtools [String] path to alternative installation of samtools
      # @param samtools [String] path to alternative installation of bcftools
      # @return [Bio::DB::Sam] a new `Bio::DB::Sam` object
      def initialize(args)
        @fasta = args[:fasta]
        @bams = args[:bam]
        @samtools = args[:samtools] || File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','samtools')
        @bcftools = args[:bcftools] || File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','bcftools')
        
        @files = [@files] if @files.instance_of?(String)
        @last_command = nil
        raise ArgumentError, "Need Fasta and at least one BAM or SAM" if not @fasta or not @bams
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
        region = ""
        if opts[:chr] and opts[:start] and opts[:stop]
          region = "#{opts[:chr]}:#{opts[:start]}-#{opts[:stop]}"
          [:chr, :start, :stop].each {|o| opts.delete(o)}
        end
        if opts[:at]
          opts["@"] = opts[:at]
          opts.delete(:at)
        end
        
        if opts[:one]
          opts["1"] = opts[:one]
          opts.delete(:one)
        end
                
        command = form_opt_string(@samtools, "view", opts, [:b, :h, :H, :S, :u, "1", :x, :X, :c, :B]) + " " + region
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
        return result
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
      
      private
      
      # returns a command string from a program
      # @param program [Symbol] either `:samtools` or `:bcftools`
      # @param opts [Hash] the options hash
      # @param singles `flag` options [Array] the options in `opts` that are single options 
      def form_opt_string(prog, command, opts, singles)
        opts_string = commandify(opts, singles)
        "#{prog} #{command} #{opts_string} #{@bams.join(' ')}"
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
          while line = pipe.gets 
            next if skip_comments and line[0] == comment_char
            yield klass.new(line.chomp)
          end
        elsif type == :binary
          while c = pipe.gets(nil)
            yield c
          end
        end
        pipe.close
      end
      
    end
  end
end
