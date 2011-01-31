#require 'rubygems'
#require'ffi'
#require 'bio/db/sam/bam'
module Bio
  module DB
    module SAM
      module Tools
        extend FFI::Library
        #ffi_lib "#{File.join(File.expand_path(File.dirname(__FILE__)),'external','libbam.dylib')}"
        ffi_lib Bio::DB::Sam::Library.filename

        attach_function :fai_build, [ :string ], :int
        attach_function :fai_destroy, [ :pointer ], :void
        attach_function :fai_load, [ :string ], :pointer
        attach_function :fai_fetch, [ :pointer, :string, :pointer ], :string
        attach_function :faidx_fetch_nseq, [ :pointer ], :int
        attach_function :faidx_fetch_seq, [ :pointer, :string, :int, :int, :pointer ], :string
      end
    end
  end
end