#
#  sam.rb
#  
#
#  Created by Ricardo Ramirez on 3/25/10.
#
# require 'rubygems'
# require'ffi'
# require 'bio/db/sam/bam'
module Bio
  module DB
    module SAM
      module Tools
        extend FFI::Library
        
#        ffi_lib "#{File.join(File.expand_path(File.dirname(__FILE__)),'external','libbam.#{lib_os}')}"
         ffi_lib Bio::DB::Sam::Library.filename

#        typedef struct {
#        	int32_t n_targets;
#        	char **target_name;
#        	uint32_t *target_len;
#        	void *dict, *hash, *rg2lib;
#        	int l_text;
#        	char *text;
#        } bam_header_t;
        class BamHeaderT < FFI::Struct
          layout(
          :n_targets, :int32_t,
          :target_name, :pointer,
          :target_len, :pointer,
          :dict, :pointer,
          :hash, :pointer,
          :rg2lib, :pointer,
          :l_text, :int,
          :text, :pointer
          )
          def text=(str)
            @text = FFI::MemoryPointer.from_string(str)
            self[:text] = @text
          end
          def text
            @text.get_string(0)
          end

        end
        
        class SamfileTX < FFI::Union
          layout(
          :tamr, :pointer, #Text file, read.
          :bam,  :pointer, #bamFile,
          :tamw, :pointer #Text file, write. 
          )
        end
#        typedef struct {
#        	int type;
#        	union {
#        		tamFile tamr;
#        		bamFile bam;
#        		FILE *tamw;
#        	} x;
#        	bam_header_t *header;
#        } samfile_t;
        class SamfileT < FFI::Struct
          layout(
          :type, :int,          
          :x, SamfileTX,
          :header, :pointer
          #:header, BamHeaderT
          )
        end
        
        

        attach_function :samclose, [ :pointer ], :void
        attach_function :samread, [ :pointer, :pointer ], :int
        attach_function :samopen, [ :string, :string, :pointer ], :pointer
        attach_function :samwrite, [ :pointer, :pointer ], :int
        attach_function :sampileup, [ :pointer, :int, :bam_pileup_f, :pointer ], :int
        attach_function :samfaipath, [ :string ], :string
      end
    end
  end
end