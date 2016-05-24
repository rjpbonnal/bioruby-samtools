module Bio
  class DB
    class Sam
      attr_accessor :bam, :fasta, :samtools, :bcftools, :last_command
      attr_accessor :minumum_ratio_for_iup_consensus
      attr_reader :cached_regions
      #attr_accessor :pileup_cache
      @minumum_ratio_for_iup_consensus = 0.20
      BASE_COUNT_ZERO =  {:A => 0, :C => 0, :G => 0,  :T => 0}

      #Creates a new Bio::DB::Sam object
      #* fasta [String] - the path to the Fasta reference sequence
      #* bam [String] - path to bam files
      #* samtools [String] - path to alternative installation of samtools
      #* bcftools [String] - path to alternative installation of bcftools
      #* returns [Bio::DB::Sam] a new `Bio::DB::Sam` object
      def initialize(args)
        @fasta = args[:fasta]
        @bam = args[:bam]
        @samtools = args[:samtools] || File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','samtools')
        @bcftools = args[:bcftools] || File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','bcftools')

        @files = [@files] if @files.instance_of?(String)

        @last_command = nil
        raise ArgumentError, "Need Fasta and at least one BAM or SAM" if not @fasta or not @bam
        raise IOError, "File not found #{@files}" if not files_ok?
        @bams = [@bams] if @bams.instance_of? String

      end

      #backward compatibility method, returns true if file exists otherwise, complains and quits.
      def open
        files_ok?
      end

      #runs the samtools view command   
      #* b - output BAM
      #* h - print header for the SAM output
      #* H - print header only (no alignments)
      #* S - input is SAM
      #* u - uncompressed BAM output (force -b)
      #* one - fast compression (force -b)
      #* x - output FLAG in HEX (samtools-C specific)
      #* X - output FLAG in string (samtools-C specific)
      #* c - print only the count of matching records
      #* B - collapse the backward CIGAR operation
      #* at - INT number of BAM compression threads [0]
      #* L - FILE output alignments overlapping the input BED FILE [null]
      #* t - FILE list of reference names and lengths (force -S) [null]
      #* T - FILE reference sequence file (force -S) [null]
      #* o - FILE output file name [stdout]
      #* R - FILE list of read groups to be outputted [null]
      #* f - INT required flag  0 for unset [0]
      #* F - INT filtering flag  0 for unset [0]
      #* q - INT minimum mapping quality [0]
      #* l - STR only output reads in library STR [null]
      #* r - STR only output reads in read group STR [null]
      #* s - FLOAT fraction of templates to subsample; integer part as seed [-1]
      #* chr - name of reference sequence to get alignments from
      #* start - start position on reference sequence
      #* stop - end postion on reference sequence
      def view(opts={},&block)
        region = String.new
        if opts[:chr] and opts[:start] and opts[:stop]
          has_e = self.has_entry? opts[:chr]
          raise Exception.new(), "[view] The sequence #{opts[:chr]} is not in the bam file" unless self.has_entry? opts[:chr] 
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
        command = String.new
        command = form_opt_string(@samtools, 'view', opts, [:b, :h, :H, :S, :u, '1', :x, :X, :c, :B]) 
        command = command + " '#{region}'" if region.size > 0
        @last_command = command
        type = (opts[:u] or opts[:b]) ? :binary : :text
        klass = (type == :binary) ? String : Bio::DB::Alignment
        yield_from_pipe(command, klass, type, &block)
      end

      #fetches a subsequence and calls code block
      #* chr - the reference name for the subsequence
      #* start - the start position for the subsequence
      #* stop - the stop position for the subsequence
      #* &block - the the block of code to execute
      def fetch(chr, start,stop, &block)
       
        view(
        :chr => chr,
        :start => start,
        :stop => stop, 
        &block  
        )
      end

      alias_method :fetch_with_function, :fetch

      #returns an array of coverage for each location for which there are mapped reads
      #* chr - the reference name
      #* start - the start position 
      #* length - the length of the region queried
      def chromosome_coverage(chr,start,length)
        result = []
        region = "#{chr}:#{start}-#{start + length}"
        self.mpileup(:r => region) do |p|
          result << p.coverage
        end
        result
      end


      #returns an svg file or object, plotting coverage for each location for which there are mapped reads
      #* chr - the reference name
      #* start - the start position 
      #* length - the length of the region queried
      #OPTIONS
      #* bin - the amount of bins to split the histogram into. The arithmetic mean score for each bin will be plotted. [default 30 bins]
      #* svg - a file to write the svg image to [default a String object containing the SVG]
      def plot_coverage(chr,start,length, opts={})
        chr = opts[:chr] if chr.nil?
        start = opts[:start] if start.nil?
        length = opts[:length] if length.nil?
        if opts[:bin]
          bin = length/opts[:bin]
        else
          bin = length/30
        end
        result = []
        region = "#{chr}:#{start}-#{start + length}"
        self.mpileup(:r => region) do |p|
          result << p.coverage
        end
        p = Bio::Graphics::Page.new(:width => 1000, 
        :height => 200, 
        :number_of_intervals => 10,
        :font_size => 14
        )
        default_options = {:glyph => :histogram, 
        :stroke => 'black',
        :fill_color => 'gold',
        :track_height => 150,
        :name => 'read coverage', 
        :label => true, 
        :stroke_width => '1', 
        :x_round => 1,
        :y_round => 1 }
        opts = default_options.merge(opts)
        
        data_track = p.add_track(opts)
        index = 0;        
        result.each_slice(bin) {|slice| 
          #result.each_with_index {|val, index|
          data_feature = Bio::Graphics::MiniFeature.new(:start => start + index,
          :end => (start + index + bin),
          :segment_height => slice.inject{|sum,x| sum + x }.to_f / slice.size)
          data_track.add(data_feature)
          index+=bin
        }
        if opts[:svg]
          svg = opts[:svg].to_s
          p.write(svg)
        else
          return p.get_markup
        end


      end

      #returns the average coverage over the region queried
      #* chr - the reference name
      #* start - the start position 
      #* length - the length of the region queried
      def average_coverage(chr,start,length)
        arr = self.chromosome_coverage(chr,start,length)
        arr.inject{ |sum, el| sum + el }.to_f / arr.size
      end

      #returns a Bio::DB::Pileup or Bio::DB::VCF object
      #* region - Only generate pileup in region [chrom:start-stop] 
      #* illumina_quals - Assume the quality is in the Illumina 1.3+ encoding
      #* count_anomalous - Do not skip anomalous read pairs in variant calling
      #* no_baq - Disable probabilistic realignment for the computation of base alignment quality (BAQ). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.
      #* adjust_mapq - [INT] Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. [0] 
      #* max_per_bam_depth - [INT] At a position, read maximally INT reads per input BAM. [250] 
      #* extended_baq - Extended BAQ computation. This option helps sensitivity especially for MNPs, but may hurt specificity a little bit.
      #* exclude_reads_file - [FILE] exclude read groups listed in FILE [null]
      #* list_of_positions - [FILE] BED or position list file containing a list of regions or sites where pileup or BCF should be generated [null]
      #* mapping_quality_cap - [INT] cap mapping quality at INT [60]
      #* ignore_rg - ignore read group tags
      #* min_mapping_quality - [INT] skip alignments with mapQ smaller than INT [0]
      #* min_base_quality - [INT] skip bases with baseQ/BAQ smaller than INT [13]   
      #* ##following options are for the -g -u option
      #* genotype_calling - generate BCF output (genotype likelihoods)
      #* uncompressed_bcf - generate uncompress BCF output
      #* extension_sequencing_probability - [INT] Phred-scaled gap extension seq error probability [20]
      #* homopolymer_error_coefficient - [INT] coefficient for homopolymer errors [100]
      #* no_indels - do not perform indel calling
      #* skip_indel_over_average_depth - [INT] max per-sample depth for INDEL calling [250]
      #* gap_open_sequencing_error_probability - [INT] Phred-scaled gap open sequencing error probability [40]
      #* platforms - [STRING] comma separated list of platforms for indels [all]
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
        opts = Hash.new
        #To remove any unwanted options. 
        long_opts.each_pair do |k,v|
          opts[v] = temp_opts[v] if temp_opts.has_key?(v)
        end

        #        opts = temp_opts
        opts[:u] = true if opts[:g] #so that we always get uncompressed output
        opts.delete(:g)

        opts[:f] = @fasta

        #TOODO: reduce the string handling
        query = opts[:r].to_s
        query = opts[:r].to_region.to_s if opts[:r].respond_to?(:to_region)
        if not query.nil? and query.size > 0
          raise Exception.new(), "The sequence #{query} is not in the bam file"  unless has_region? query 
        end
        opts[:r] = query
        
        if opts[:six]
          opts["6"] = nil
          opts.delete(:six)
        end

        command = form_opt_string(@samtools, "mpileup", opts, [:R, :B, :E, "6", :A, :g, :u, :I] )
        puts "Running: #{command}" if $VERBOSE
        if opts[:u]
          command = command + " | #{@bcftools} view -cg -"
        end
        
        klass = opts[:u] ? Bio::DB::Vcf : Bio::DB::Pileup
        @last_command = command
        yield_from_pipe(command, klass, :text, &block)

      end

      #fetches a subsequence from a reference genome and option returns it as a Bio::Sequence::NA object
      #* chr -  [STRING] the reference name for the subsequence
      #* start - [INT] the start position for the subsequence
      #* stop - [INT] the stop position for the subsequence
      #* as_bio - boolean stating if the returned object should be a Bio::Sequence::NA object
      def fetch_reference(chr,start,stop, opts={:as_bio => false})
        raise Exception.new(), "The sequence #{chr} is not in the bam file" unless has_entry? chr
        seq = ""
        unless @fasta #We return a string of Ns if we don't know the reference. 
          seq = "n" * (stop-start) 
        else
          command = "#{@samtools} faidx \"#{@fasta}\" '#{chr}:#{start}-#{stop}'"
          puts "Running: #{command}" if $VERBOSE
          @last_command = command
          seq = ""
          yield_from_pipe(command, String, :text ) {|line| seq = seq + line unless line =~ /^>/}
        end

        if opts[:as_bio]
          seq = Bio::Sequence::NA.new(seq).to_fasta("#{chr}:#{start}-#{stop}")
        end
        seq
      end

      #Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence. If no region is specified, faidx will index the file and create <ref.fasta>.fai on the disk. If regions are speficified, the subsequences will be retrieved and printed to stdout in the FASTA format.
      #Options - if a subsequence is required
      #* chr - [STRING] the reference name of the subsequence
      #* start - [INT] the start position for the subsequence
      #* stop - [INT] the stop position for the subsequence
      def faidx(opts={})
        if opts.has_key?(:chr) and opts.has_key?(:start) and opts.has_key?(:stop)
          opts={:as_bio => false}
          self.fetch_reference(:chr,:start,:stop,opts)
        else
          command = "#{@samtools} faidx \"#{@fasta}\""
          @last_command = command
          system(command)
        end
      end

      #Index sorted alignment for fast random access. Index file <aln.bam>.bai will be created of no out_index is provided.
      #* out_index - [STRING] name of index
      def index(opts={})
        command = "#{@samtools} index \"#{@bam}\" #{opts[:out_index]}"
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        system(command)
      end

      #Fill in mate coordinates, ISIZE and mate related flags from a name-sorted alignment
      #* out_bam name of outfile
      #* r - remove unmapped reads and secondary alignments
      def fix_mates(opts={})
        #opts.merge!({:out_index=>nil})
        remove_reads = ""
        if opts[:r]
          remove_reads = "-r"
        end
        command = "#{@samtools} fixmate #{remove_reads} \"#{@bam}\" #{opts[:out_bam]}"
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        system(command)
      end

      alias_method :fixmate, :fix_mates

      #generate simple stats with regard to the number and pairing of reads mapped to a reference
      def flag_stats(opts={})
        command = form_opt_string(@samtools, "flagstat", opts, [])
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        strings = []
        yield_from_pipe(command,String) {|line| strings << line.chomp}
        strings
      end

      alias_method :flagstat, :flag_stats

      #Retrieve and print stats in the index file. The output is TAB delimited with each line consisting of reference sequence name, sequence length, number of mapped reads and number unmapped reads.
      def index_stats
       return @stats if @stats
        stats = {}   
        command = form_opt_string(@samtools, "idxstats", {}, [])
        @last_command = command
        puts "Running: #{command}" if $VERBOSE
        yield_from_pipe(command, String, :text, true, "#") do |line|
          info = line.chomp.split(/\t/)
          stats[ info[0] ] = {:length => info[1].to_i, :mapped_reads => info[2].to_i, :unmapped_reads => info[3].to_i }
        end
        @stats = stats
        return @stats
      end

      alias_method :idxstats, :index_stats
      
      #Retrive a hash with all the regions, with the region id as index or runs the function on each region
      def each_region
        index_stats 
        if @regions 
          return @regions unless block_given? 
        else
          @regions = Hash.new
        end
        index_stats.each do |k,v|
          reg = Bio::DB::Fasta::Region.new
          reg.entry = k
          reg.start = 1
          reg.end = v[:length]
          reg.orientation = :forward
          @regions[k] = reg unless @regions[k]
          yield reg if block_given?
        end
        @regions
      end
      
      #Tells if the bam file contains the entry. It has to be indexed.
      def has_entry?(entry)
         index_stats.has_key?(entry)
    #    puts "#{entry} #{@stats.inspect}"
      #  index_stats
      end
      
      def has_region?(region)
        index_stats
        reg=Bio::DB::Fasta::Region::parse_region(region)
        return 0 unless has_entry? (reg.entry) 
         len = @stats[reg.entry][:length]
         reg.start > 0 and reg.end <= len
      end

      #Merge multiple sorted alignments
      #* n - sort by read names
      #* r - attach RG tag (inferred from file names)
      #* u - uncompressed BAM output
      #* f - overwrite the output BAM if exist
      #* one - compress level 1
      #* l  - [INT] compression level, from 0 to 9 [-1]
      #* at - [INT] number of BAM compression threads [0]
      #* R - [STRING] merge file in the specified region STR [all]
      #* h - [FILE] copy the header in FILE to <out.bam> [in1.bam]
      #* out - [FILE] out file name
      #* bams - [FILES] or Bio::DB::Sam list of input bams, or Bio::DB::Sam objects
      def merge(opts={})
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
        options = commandify(opts, [:n, :r, :u, :f, '1'] )
        command = "#{@samtools} merge #{options} #{out} #{bam_list}"

        @last_command = command
        puts "Running: #{command}" if $VERBOSE
        system(command)

      end

      #Concatenate BAMs. The sequence dictionary of each input BAM must be identical.
      #* h - header.sam
      #* out -[FILE] out file name
      #* bams -[FILES] or Bio::DB::Sam list of input bams, or Bio::DB::Sam objects
      def cat(opts={})
        bam_list = opts[:bams].collect do |b|
          b.bam rescue b
        end.join(' ')
        opts.delete(:bams)
        options = commandify(opts, [:h] )
        command = "#{@samtools} cat #{options} -o #{out} #{bam_list}"
        puts command if $VERBOSE
        @last_command = command
        system(command)

      end

      #* program  - one of 'samtools' 'bcftools'
      #* command - one of the commands relevant to the program
      def self.docs(program, command)
        return "program must be 'samtools' or 'bcftools'" if not ['samtools', 'bcftools'].include? program
        command = "#{program} #{command}"
        `#{command}`
      end

      #Remove potential PCR duplicates: if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality.
      #* s - rmdup for SE reads
      #* S - treat PE reads as SE in rmdup (force -s)
      #* out - [FILE] output bam
      def remove_duplicates(opts={})
        out = opts[:out]
        opts.delete(:out)
        command = "#{form_opt_string(@samtools, "rmdup", opts, [:s, :S])} #{out} \"#{@bam}\""
        @last_command = command
        system(command)
      end

      alias_method :rmdup, :remove_duplicates

      #Sort alignments by leftmost coordinates
      #* n - sort by read name
      #* f - use <out.prefix> as full file name instead of prefix
      #* o - final output to stdout returns bio::db::alignment
      #* l - [INT]  compression level, from 0 to 9 [-1]
      #* at - [INT] number of sorting and compression threads [1]
      #* m - [INT] max memory per thread; suffix K/M/G recognized [768M]
      #* prefix - [STRING] prefix for output bamfile
      def sort(opts={})
        if !opts.has_key?(:prefix)
          opts.merge!({:prefix => "sorted"})
        end
        prefix = opts[:prefix]
        opts.delete(:prefix)
        command = form_opt_string(@samtools, "sort", opts, [:n, :f, :o])
        command = command + " " + prefix
        @last_command = command
        puts "Running: #{command}" if $VERBOSE
        if opts[:o]
          yield_from_pipe(command, Bio::DB::Alignment)
        else
          system(command)
        end
      end

      #used to generate a text alignment viewer
      #* d - display, output as (H)tml or (C)urses or (T)ext 
      #* p - [chr:pos] go directly to this position
      #* s - [STR] display only reads from this sample or group
      def tview(opts={})
        if opts[:d]
          opts['d'] = opts[:d]
          opts.delete(:d)
        end
        if opts[:p]
          opts['p'] = opts[:p]
          opts.delete(:p)
        end
        if opts[:s]
          opts['s'] = opts[:s]
          opts.delete(:s)
        end
        command = "#{form_opt_string(@samtools, "tview", opts)}"
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        system(command)
      end

      #Replace the header of the current bam file with the header in header_sam
      #* header_sam - the sam file from which the new header will be taken
      #* out - [FILE] output bam file
      def reheader(header_sam, opts={})
        if opts.has_key?(:out)
          out=opts[:out]
          command = "#{@samtools} reheader #{header_sam} \"#{@bam}\" > #{out}"
        else
          command = "#{@samtools} reheader #{header_sam} \"#{@bam}\""
        end
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        system(command)
      end

      #Generate the MD tag. If the MD tag is already present, this command will give a warning if the MD tag generated is different from the existing tag. Output SAM by default.
      #* A - When used jointly with -r this option overwrites the original base quality.
      #* e - Convert a the read base to = if it is identical to the aligned reference base. Indel caller does not support the = bases at the moment.
      #* u - Output uncompressed BAM
      #* b - Output compressed BAM
      #* S - The input is SAM with header lines
      #* C - [INT] Coefficient to cap mapping quality of poorly mapped reads. See the pileup command for details. [0]
      #* r - Compute the BQ tag (without -A) or cap base quality by BAQ (with -A).
      #* E - Extended BAQ calculation. This option trades specificity for sensitivity, though the effect is minor. 
      def calmd(opts={}, &block)
        command = form_opt_string(@samtools, "calmd",  opts, [:E, :e, :u, :b, :S, :r] )+ " " + @fasta
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        type = :text 
        klass = Bio::DB::Alignment
        yield_from_pipe(command, klass, type, true, "@",&block)
      end

      #Identifies target regions by examining the continuity of read depth, computes haploid consensus sequences of targets and outputs a SAM with each sequence corresponding to a target. When option -f is in use, BAQ will be applied. 
      #* Q - [INT] Minimum base quality for a base to be considered [13] 
      #* i - in penalty
      #* 0 - em0
      #* 1 - em1
      #* 2 - em2
      #* f - reference
      def targetcut(opts={})
        if opts[:f]
          opts['f'] = @fasta
          opts.delete(:s)
        end

        command = "#{form_opt_string(@samtools, "targetcut", opts, [] )}"
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        system(command)
      end

      #Call and phase heterozygous SNPs
      #* A - Drop reads with ambiguous phase.
      #* b - [STR] Prefix of BAM output. When this option is in use, phase-0 reads will be saved in file STR.0.bam and phase-1 reads in STR.1.bam. Phase unknown reads will be randomly allocated to one of the two files. Chimeric reads with switch errors will be saved in STR.chimeric.bam. [null]
      #* F - Do not attempt to fix chimeric reads.
      #* k - [INT] Maximum length for local phasing. [13]
      #* q - [INT] Minimum Phred-scaled LOD to call a heterozygote. [40]
      #* Q - [INT] Minimum base quality to be used in het calling. [13] 
      def phase(opts={})
        command = "#{form_opt_string(@samtools, "phase", opts, [:A, :F] )}"
        puts "Running: #{command}" if $VERBOSE
        @last_command = command
        system(command)
      end


      #returns an array for each position with [sequence_name, position, depth]
      #* b - list of positions or regions in BED format
      #* l - [INT] minQLen
      #* q - [INT] base quality threshold
      #* Q - [INT] mapping quality threshold
      #* r - [chr:from-to] region
      def depth(opts={})
        command = form_opt_string(@samtools, "depth", opts)
        @last_command = command
        system(command)
      end

      #Returns the pipelup of a region, encapsulated as a Bio::DB::Fasta::Region object.  
      #The opts are the same as for mpileup
      def fetch_region(opts={})   
        region = opts[:r] ? opts[:r] : opts[:region]
        opts[:r] = region
        opts[:region] = region
        reg =  Bio::DB::Fasta::Region.parse_region(region.to_s)
        reg.reference = self.fetch_reference(region.entry, region.start, region.end).downcase
        tmp = Array.new
        mpileup(opts) do | pile | 
          #  puts pile
          tmp << pile 
          yield pile if block_given?
        end
        reg.pileup =  tmp
        reg.calculate_stats_from_pile(opts)
        reg
      end

      #Same as mpilup, but it caches the pileup, so if you want several operations on the same set of regions
      #the pile for different operations, it won't execute the mpilup command several times
      #Whenever you finish using a region, call mpileup_clear_cache to free the cache
      #The argument Region is required, as it will be the key for the underlying hash. 
      #We asume that the options (other than the region) are constant. If they are not, the cache mechanism may not be consistent. 
      #
      #TODO: It may be good to load partially the pileup
      def mpileup_cached (opts={})      
        raise Exception.new(), "A region must be provided" unless opts[:r] or opts[:region]
        @cached_regions = Hash.new unless @cached_regions
        region = opts[:r] ? opts[:r] : opts[:region]
        @cached_regions[region.to_s] = fetch_region(opts) unless @cached_regions[region.to_s]
        if block_given?
          @cached_regions[region.to_s].pileup.each do | pile |
            yield pile 
          end  
        end
        region.pileup
      end


      #Clears the pileup cache. If a region is passed as argument, just the specified region is removed
      #If no region is passed, the hash is emptied
      def mpileup_clear_cache (region)
        return unless @cached_regions
        if region
          @cached_regions[region.to_s] = nil
        else
          @cached_regions.clear
        end
      end

      def bedcov(opts={})
        bed = opts[:bed]
        #bam = opts[:bam]
        if opts.has_key?(:out)
          out=opts[:out]
          command = "#{@samtools} bedcov \"#{bed}\" \"#{@bam}\" > \"#{out}\""
        else
          command = "#{@samtools} bedcov \"#{bed}\" \"#{@bam}\""
        end
         puts "Running: #{command}" if $VERBOSE
        #puts command
        @last_command = command
        system(command)
      end


      #Extract the reads that align to a region
      #* region [String] - Region to extract (chromosome:start-end)
      #* fastq - [INT] fastq file where to print. If empty, prints to stdout
      #* q - [INT] base quality threshold
      # Not tested yet
      def extract_reads(opts={})
        opts[:region] = Bio::DB::Fasta::Region.parse_region( opts[:region] .to_s)  unless opts[:region].class == Bio::DB::Fasta::Region
        fastq_filename = opts[:fastq]

        out = $stdout 
        print_fastq = Proc.new do |alignment|
          out.puts "@#{alignment.qname}"
          out.puts "#{alignment.seq}"
          out.puts "+#{alignment.qname}"
          out.puts "#{alignment.qual}"
        end

        if fastq_filename
          out = File.open(fastq_filename, "w")
        end
        fetch_with_function(chromosome, qstart, qstart+len,  print_fastq)
        out.close if fastq_filename
      end
      
       # checks existence of files in instance
      def files_ok?
        [@fasta, @sam, @bam].flatten.compact.each {|f| return false unless File.exists? f }
        true
      end
      
      #Returns true if the .bai exists. It doesn't validate if it is valid. 
      def indexed?
        File.exists? @bam and File.exists? "#{@bam}.bai"
      end
        
      private
      #Returns Process::Status with the execution status. If run in a $VERBOSE environment, stderr of the process
      #is forwarded to the default stdout
      def yield_from_pipe(command, klass, type=:text, skip_comments=true, comment_char="#", &block)
        puts "[yield_from_pipe] #{command}" if $VERBOSE
        stdin, pipe, stderr, wait_thr = Open3.popen3(command)
        pid = wait_thr[:pid]  # pid of the started process.       
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
        exit_status = wait_thr.value  # Process::Status object returned.
        puts "Running: #{command}" if $VERBOSE 
        stdin.close
        pipe.close
        stderr.close
        return exit_status
      end


      # returns a command string from a program
      # @param program [Symbol] either `:samtools` or `:bcftools`
      # @param opts [Hash] the options hash
      # @param singles `flag` options [Array] the options in `opts` that are single options 
      def form_opt_string(prog, command, opts, singles=[])
        opts_string = commandify(opts, singles)
        "#{prog} #{command} #{opts_string} \"#{@bam}\""
      end

      # turns an opts hash into a string
      def commandify(opts, singles)
        list = []
        opts.each_pair do |tag,value|
          value = "\"#{value}\""
          value = "" if singles.include?(tag)

          list << "-#{tag.to_s} #{value}" 
        end
        list.join(" ")
      end
    end
  end
end
