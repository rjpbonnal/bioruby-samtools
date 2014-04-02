#Module to hold the information about the fasta file

module Bio::DB::Fasta
  class Index
    include Enumerable
    attr_reader :entries
    
    def initialize
      @entries=[]
      @entries_map = Hash.new
    end
    
    #This doesnt validate if you are adding the same entry twice. I may add
    #a validation for that. 
    def << (entry)
      @entries << entry
      @entries_map[entry.id] = entry
    end
    
    def each(&block)
      @entries.entries(&block)
    end  
    
    def length
       @entries.length
    end
    
    #Returns a new Index just with the specified range, as if it was an Array. 
    #The return object is of type Index. 
    def [](args)
      tmp = @entries[args]
      new_index = Index.new
      tmp.each do | entry |
        @new_index << entry
      end       
    end
    
    def region_for_entry(entry)
       @entries_map[entry]
    end
  end
  
  class Entry
    attr_reader :id, :length
    
    def initialize(id, length)
      @id=id
      @length=length.to_i
    end
    
    def get_full_region
      reg = Region.new
      reg.entry = id
      reg.start = 0
      reg.end = @length
      reg.orientation = :forward
      reg
    end
    
    def to_region
      get_full_region
    end
  end
  
  #Class to wrap a region of a chromosome
  class Region
    attr_accessor :entry, :start, :end, :orientation
    
    def to_s
      string = @entry + ":" + @start.to_s + "-" + @end.to_s 
      string
    end
    
    def self.parse_region(reg_str)
        string = reg_str.delete("'")
        fields_1 = string.split(":")
        fields_2 = fields_1[1].split("-")
        raise FastaDBException.new(), "Invalid region. #{string}" if fields_1.length != 2 || fields_2.length != 2
        
        reg = Region.new
        reg.entry = fields_1[0]
        reg.start = fields_2[0].to_i
        reg.end = fields_2[1].to_i
        
        if reg.end < reg.start 
          reg.orientation = :reverse
        else
          reg.orientation = :forward
        end
        reg
    end
    
    def size
      @end - @start
    end
  
  end
  
  class FastaDBException < StandardError; end
  
  #Class that holds the fasta file. It is used as a database. It heavily relies ond samtools. 
  class FastaFile
    
    attr_reader :index, :fasta_path
    
    def FastaFile.finalize(id)
      #id.close()
      #puts "Finalizing #{id}  at #{Time.new}"       
    end
    
    def initialize(fasta_filename)
      @fasta_path = fasta_filename
      raise FastaDBException.new(), "No path for the refernce fasta file. " if @fasta_path.nil?
      @fasta_index = Bio::DB::SAM::Tools.fai_load(@fasta_path)
      if @fasta_index.null? then
        $stderr.puts "Generating index for: " + @fasta_path
        Bio::DB::SAM::Tools.fai_build(@fasta_path)
        @fasta_index =  Bio::DB::SAM::Tools.fai_load(@fasta_path)
        raise FastaDBException.new(), "Unable to generate fasta index for: " + @fasta_path if @fasta_index.nil? ||  @fasta_index.null?
      end
      ObjectSpace.define_finalizer(self,  self.class.method(:finalize).to_proc)
    end
    
    def load_fai_entries()
      return  @index.length if @index
      @index = Index.new
      fai_file = @fasta_path + ".fai"
      File.open(fai_file).each do | line |
        fields = line.split("\t")
        @index << Entry.new(fields[0], fields[1])
        
      end     
      @index.length
    end
    
    def close()
       Bio::DB::SAM::Tools.fai_destroy(@fasta_index) unless @fasta_index.nil? || @fasta_index.null?
       @fasta_index = nil
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
    
    
    #The region needs to have a method to_region or a method to_s that ha the format "chromosome:start-end" as in samtools
    def fetch_sequence(region)
      
        raise FastaDBException.new(), "No fasta index for " if @fasta_index.nil? || @fasta_index.null? 
        query = region.to_s
        query = region.to_region.to_s if region.respond_to?(:to_region) 
        
        len = FFI::MemoryPointer.new :int
        str = Bio::DB::SAM::Tools.fai_fetch(@fasta_index, query, len)
        raise FastaDBException.new(), "Unable to get sequence for reference: " + query if str.nil?
        reference = Bio::Sequence.auto(str)
        
        # 
        
        if region.orientation == :reverse
          #puts "reversing! #{reference.to_s}"
          reference.reverse_complement!()
        end
        reference
    end
    
     
  end
end