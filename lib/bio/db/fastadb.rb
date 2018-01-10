#Module to hold the information about the fasta file

module Bio::DB::Fasta
  #This class contains the entries in a fasta, as generated by samtools faidx
  class Index
    include Enumerable
    attr_reader :entries

    def initialize
      @entries=[]
      @entries_map = Hash.new
    end

    #This doesnt validate if you are adding the same entry twice. I may add
    #a validation for that. 
    def <<(entry)
      @entries << entry
      @entries_map[entry.id] = entry
    end

    def each(&block)
      @entries.entries(&block)
    end  
    #Total number of entries
    def length
      @entries.length
    end
    alias_method :size, :length

    #Returns a new Index just with the specified range, as if it was an Array. 
    #The return object is of type Index. 
    def [](args)
      tmp = @entries[args]
      @new_index = Index.new
      tmp.each do | entry |
        @new_index << entry
      end       
    end

    #Gets the Region object for the full length of the sequence 
    #name queried.
    def region_for_entry(entry)
      @entries_map[entry]
    end
  end

  class Entry
    attr_reader :id, :length, :line_bases, :line_length, :offset
    alias_method :size, :length
    def initialize(id, length, offset = 0 , line_bases= 0 , line_length = 0 )
      @id=id
      @length=length.to_i
      @offset = offset.to_i
      @line_bases  = line_bases.to_i
      @line_length = line_length.to_i
    end

    def get_base_coordinate(coordinate)
      lines_for_offset = coordinate / line_bases
      line_offset = coordinate % line_bases
      #puts "get_base_coordinate"
      #puts "Coordinate: #{coordinate}"
      #puts "lines_for_offset: #{lines_for_offset}"
      #puts "line pffset: #{line_offset}"
      #puts self.inspect
      pointer = offset + (line_length * lines_for_offset) + line_offset - 1
      pointer
    end

    def get_full_region
      reg = Region.new
      reg.entry = id
      reg.start = 1
      reg.end = @length
      reg.orientation = :forward
      reg
    end

    alias_method  :to_region, :get_full_region

  end

  #Class to wrap a region of a chromosome
  class Region
    BASE_COUNT_ZERO =  {:A => 0, :C => 0, :G => 0,  :T => 0}
    attr_accessor :entry, :start, :end, :orientation

    attr_accessor :pileup, :average_coverage, :snps, :reference, :allele_freq, :consensus, :coverages, :bases, :total_cov, :called
  
    def initialize(args ={})
      @entry = args[:entry]
      @start = args[:start]
      @end = args[:end]
      @orientation = args[:orientation]
    end
  
    #TODO: Debug, as it hasnt been tested in the actual code. 
    def allele_freq_for_base(base)
      @all_ratios = Hash.new unless @all_ratios
      unless @all_ratios[base]
        ratios = Array.new
        for i in (0..region.size-1)
          ratios << @allele_freq[i][base]
        end
        @all_ratios[base] = ratios
      end
      @all_ratios[base]
    end

    alias_method :base_ratios_for_base, :allele_freq_for_base
    alias_method :base_ratios, :allele_freq

    #Calculates the concensus, base ratios, coverages and total coverages in the region
    #* min_cov minimum coverage to make a call (default 0)
    #* min_per minimum representation to make make a call. If more than one base 
    #  can be called, the IUAPC ambiguity code is returned
    def calculate_stats_from_pile(opts={})
      min_cov = opts[:min_cov] ? opts[:min_cov] : 0
      min_per =  opts[:min_per] ? opts[:min_per] : 0.20
      self.called = 0
      reference = self.reference.downcase

      self.allele_freq = Array.new(self.size, BASE_COUNT_ZERO) 
      self.bases = Array.new(self.size, BASE_COUNT_ZERO) 
      self.coverages = Array.new(self.size, 0)
      self.total_cov = 0

      self.pileup.each do | pile |

        if pile.coverage > min_cov
          self.allele_freq[pile.pos - self.start ] = pile.allele_freq
          reference[pile.pos - self.start   ] = pile.consensus_iuap(min_per).upcase
          self.coverages[pile.pos - self.start   ]  = pile.coverage.to_i
          self.bases[pile.pos - self.start       ]  = pile.bases
          self.called += 1 
        end
        #puts "#{pile.pos}\t#{bef}\t#{reference[pile.pos - region.start  - 1 ]} "
        self.total_cov += pile.coverage
      end

      self.consensus = Bio::Sequence.new(reference)
      self.consensus.na
      if self.orientation == :reverse
        self.consensus.reverse_complement!()
      end
      self.average_coverage = self.total_cov.to_f/self.size.to_f
      self
    end

    def to_s
      string = @entry + ":" + @start.to_s + "-" + @end.to_s 
      string
    end

    #Returns a region object from a string in form "name:start-end"
    def self.parse_region(reg_str)
      string = reg_str.delete("'")
      fields_1 = string.split(":")
      raise FastaDBException.new(), "Invalid region. #{string}" if fields_1.length != 2
      fields_2 = fields_1[1].split("-")
      raise FastaDBException.new(), "Invalid region. #{string}" if fields_2.length != 2

      reg = Region.new(:entry=> fields_1[0], :start=>fields_2[0].to_i, :end=>fields_2[1].to_i)

      if reg.end < reg.start 
        reg.orientation = :reverse
      else
        reg.orientation = :forward
      end
      reg
    end

    #Length of the region
    def size
      @end - @start
    end
    alias_method :length, :size

  end

  class FastaDBException < StandardError; end

  #Class that holds the fasta file. It is used as a database. 
  class FastaFile
    attr_reader :fasta_path

    #Initialize the fasta file. If the fai file doesn't exists, it is generated at startup
    #* fasta path to the fasta file
    #* samtools path to samtools, if it is not provided, use the bundled version
    def initialize(fasta: nil, samtools: false)
      #puts "The arguments are: '#{fasta}':'#{samtools}'"
      @fasta_path = fasta
      @samtools = samtools  
      @index = nil
      @fasta_file = nil
      @samtools = File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','samtools') if samtools == true
      raise FastaDBException.new(), "No path for the refernce fasta file. " if @fasta_path.nil?
      @fai_file = @fasta_path + ".fai" 
      unless File.file?(@fai_file) then
        command = "#{@samtools} faidx '#{@fasta_path}'"
        @last_command = command
        system(command)
      end
    end

    #Loads the fai entries 
    def load_fai_entries()
      return  @index.length if @index
      @index = Index.new
      fai_file = @fai_file
      File.open(fai_file).each do | line |
        fields = line.split("\t")
        @index << Entry.new(fields[0], fields[1], fields[2], fields[3], fields[4])
      end     
      @index.length
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
        command = "#{@samtools} faidx #{@fasta_path}"
        @last_command = command
        system(command)
      end
    end

    def index
      return @index if @index
      if @samtools
        faidx
      else
        samtools = File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','samtools')
        #TODO: make a ruby implementations 
        command = "#{samtools} faidx #{@fasta_path}"
        @last_command = command
        system(command)
      end
        load_fai_entries
      return @index
    end

    def fetch_sequence_samtools(region)
      query = region.to_s
      query = region.to_region.to_s if region.respond_to?(:to_region) 
      command = "#{@samtools} faidx #{@fasta_path} '#{query}'"
      puts "Running: #{command}"  if $DEBUG
      @last_command = command
      seq = ""
      yield_from_pipe(command, String, :text ) {|line| seq = seq + line unless line =~ /^>/}
      seq
    end

    def fetch_sequence_native(region)
      query = region
      query = Region.parse_region(region) unless region.is_a?(Region) 
      seq = ""
      #In order to make this reentrant, if we want to make a multithreaded
      #version of this function, we need to get a lock. Currently, only one thred
      #can be assosiated with eache fastadb object
      @fasta_file = File.open(@fasta_path) unless @fasta_file
      entry = index.region_for_entry(query.entry)
     
      start_pointer  =  entry.get_base_coordinate(query.start)
      @fasta_file.seek(start_pointer, IO::SEEK_SET)
      end_pointer  =  entry.get_base_coordinate(query.end)
      to_read = end_pointer - start_pointer + 1
      seq = @fasta_file.read(to_read)
      seq.gsub!(/\s+/, '')
      seq 
    end

    #The region needs to have a method to_region or a method to_s that ha the format "chromosome:start-end" as in samtools
    def fetch_sequence(region)
      load_fai_entries
      region = Region.parse_region(region.to_s) unless region.is_a?(Region) 
      entry = index.region_for_entry(region.entry)
      raise FastaDBException.new "Entry (#{region.entry})not found in reference" unless entry
      raise FastaDBException.new "Region in invalid range (#{region}): Valid range: #{entry.to_region.to_s} has a size of #{entry.size}." if region.end > entry.size or region.start < 1
      seq = @samtools ?  fetch_sequence_samtools(region): fetch_sequence_native(region)
      reference = Bio::Sequence::NA.new(seq)
      if region.respond_to? :orientation and region.orientation == :reverse
        reference.reverse_complement!()
      end
      reference
    end

    private
    #Returns Process::Status with the execution status. If run in a $DEBUG environment, stderr of the process
    #is forwarded to the default stdout
    def yield_from_pipe(command, klass, type=:text, skip_comments=true, comment_char="#", &block)
      stdin, pipe, stderr, wait_thr = Open3.popen3(command)
      #pid = wait_thr[:pid]  # pid of the started process.       
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
      puts stderr.read if $DEBUG 
      stdin.close
      pipe.close
      stderr.close
      return exit_status
    end
  end
end
