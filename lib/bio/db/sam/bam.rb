# require 'rubygems'
# require'ffi'
# require 'bio/db/sam/faidx'
# require 'bio/db/sam/sam'
module Bio
  module DB
    module SAM
      module Tools
        extend FFI::Library
        
        #ffi_lib "#{File.join(File.expand_path(File.dirname(__FILE__)),'external','libbam.dylib')}"
        ffi_lib Bio::DB::Sam::Library.filename

        BAM_FPAIRED = 1
        BAM_FPROPER_PAIR = 2
        BAM_FUNMAP = 4
        BAM_FMUNMAP = 8
        BAM_FREVERSE = 16
        BAM_FMREVERSE = 32
        BAM_FREAD1 = 64
        BAM_FREAD2 = 128
        BAM_FSECONDARY = 256
        BAM_FQCFAIL = 512
        BAM_FDUP = 1024
        BAM_OFDEC = 0
        BAM_OFHEX = 1
        BAM_OFSTR = 2
        BAM_DEF_MASK = (4|256|512|1024)
        BAM_CIGAR_SHIFT = 4
        BAM_CIGAR_MASK = ((1 << 4) -1)
        BAM_CMATCH = 0
        BAM_CINS = 1
        BAM_CDEL = 2
        BAM_CREF_SKIP = 3
        BAM_CSOFT_CLIP = 4
        BAM_CHARD_CLIP = 5
        BAM_CPAD = 6
        class Bam1CoreT < FFI::Struct
          #uint32_t bin:16, qual:8, l_qname:8;
        	#uint32_t flag:16, n_cigar:16;
          layout(
          :tid, :int32_t,
          :pos, :int32_t,
          :bin, :uint16,
          :qual, :uint8,
          :l_qname, :uint8,
          :flag, :uint16,
          :n_cigar, :uint16,
          :l_qseq, :int32_t,
          :mtid, :int32_t,
          :mpos, :int32_t,
          :isize, :int32_t
          )
        end
        class Bam1T < FFI::Struct
          layout(
          :core, Bam1CoreT,
          :l_aux, :int,
          :data_len, :int,
          :m_data, :int,
          :data, :pointer
          )
          def qname
            #bam1_qname(b) ((char*)((b)->data))
            data = self[:data]
            data.read_string() 
          end
          
    
          Bam_NT16_Rev_Table = "=ACMGRSVTWYHKDBN"
        end
        

        
        
        attach_function :sam_open, [ :string ], :pointer
        attach_function :sam_close, [ :pointer ], :void
        attach_function :sam_read1, [ :pointer, :pointer, :pointer ], :int
        attach_function :sam_header_read2, [ :string ], :pointer
        attach_function :sam_header_read, [ :pointer ], :pointer
        attach_function :sam_header_parse, [ :pointer ], :int
        #attach_function :sam_header_parse_rg, [ :pointer ], :int This is declared in the .h file, but is not implemented
        #attach_function :bam_strmap_put, [ :pointer, :string, :string ], :int
        #attach_function :bam_strmap_get, [ :pointer, :string ], :string
        #attach_function :bam_strmap_dup, [ :pointer ], :pointer
        #attach_function :bam_strmap_init, [  ], :pointer
        #attach_function :bam_strmap_destroy, [ :pointer ], :void
        attach_function :bam_header_init, [  ], :pointer
        attach_function :bam_header_destroy, [ :pointer ], :void
        attach_function :bam_header_read, [ :pointer ], :pointer
        attach_function :bam_header_write, [ :pointer, :pointer ], :int
        attach_function :bam_read1, [ :pointer, :pointer ], :int
        attach_function :bam_write1_core, [ :pointer, :pointer, :int, :pointer ], :int
        attach_function :bam_write1, [ :pointer, :pointer ], :int
        attach_function :bam_format1, [ :pointer, :pointer ], :string
        attach_function :bam_format1_core, [ :pointer, :pointer, :int ], :string
        attach_function :bam_get_library, [ :pointer, :pointer ], :string
        class BamPileup1T < FFI::Struct
          layout(
          :b, :pointer,
          :qpos, :int32_t,
          :indel, :int,
          :level, :int,
          :is_del, :uint32,
          :is_head, :uint32,
          :is_tail, :uint32
          )
        end
        attach_function :bam_plbuf_set_mask, [ :pointer, :int ], :void
        callback :bam_pileup_f, [ :uint32, :uint32, :int, :pointer, :pointer ], :int
        attach_function :bam_plbuf_reset, [ :pointer ], :void
        attach_function :bam_plbuf_init, [ :bam_pileup_f, :pointer ], :pointer
        attach_function :bam_plbuf_destroy, [ :pointer ], :void
        attach_function :bam_plbuf_push, [ :pointer, :pointer ], :int
        attach_function :bam_pileup_file, [ :pointer, :int, :bam_pileup_f, :pointer ], :int
        attach_function :bam_lplbuf_reset, [ :pointer ], :void
        attach_function :bam_lplbuf_init, [ :bam_pileup_f, :pointer ], :pointer
        attach_function :bam_lplbuf_destroy, [ :pointer ], :void
        attach_function :bam_lplbuf_push, [ :pointer, :pointer ], :int
        attach_function :bam_index_build, [ :string ], :int
        attach_function :bam_index_load, [ :string ], :pointer
        attach_function :bam_index_destroy, [ :pointer ], :void
        
        #The function for fetching stuff
        #typedef int ( *bam_fetch_f)(
        #    const bam1_t *b,
        #    void *data);
        callback(:bam_fetch_f, [ :pointer, :pointer ], :int)
        attach_function :bam_fetch, [ :pointer, :pointer, :int, :int, :int, :pointer, :bam_fetch_f ], :int
        attach_function :bam_parse_region, [ :pointer, :pointer, :pointer, :pointer, :pointer ], :int

        #The second parameter must be only 2 characters long
        attach_function :bam_aux_get, [ :pointer, :string], :pointer
        attach_function :bam_aux2i, [ :pointer ], :int32_t
        attach_function :bam_aux2f, [ :pointer ], :float
        attach_function :bam_aux2d, [ :pointer ], :double
        attach_function :bam_aux2A, [ :pointer ], :char
        attach_function :bam_aux2Z, [ :pointer ], :string
        attach_function :bam_aux_del, [ :pointer, :pointer ], :int
        #The second parameter must be only 2 characters long
        attach_function :bam_aux_append, [ :pointer,  :string, :char, :int, :pointer ], :void
        #The second parameter must be only 2 characters long
        attach_function :bam_aux_get_core, [ :pointer,:string ], :pointer
        attach_function :bam_calend, [ :pointer, :pointer ], :uint32
        attach_function :bam_cigar2qlen, [ :pointer, :pointer ], :int32_t

        #FIXME: if we see that we need this function, implement it on ruby, seems like FFI is having problems with
        #te static inline.
        #attach_function :bam_reg2bin, [ :uint32, :uint32 ], :int
        #FIXME: if we see that we need this function, implement it on ruby, seems like FFI is having problems with
        #te static inline.
        #attach_function :bam_copy1, [ :pointer, :pointer ], :pointer
        #FIXME: if we see that we need this function, implement it on ruby, seems like FFI is having problems with
        #te static inline.
        #attach_function :bam_dup1, [ :pointer ], :pointer
      end
    end
  end
end