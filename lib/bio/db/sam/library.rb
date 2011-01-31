module Bio
  module DB
    module SAM
      module Library

        # Return the path with the file name of the library for the specific operating system
        def filename
          #TODO refactor this piece of code in all the files
          lib_os = case RUBY_PLATFORM
          when /linux/
            'a'
          when /darwin/
            'dylib'
          when /windows/
            'dll'  
          end

          File.join(File.expand_path(File.dirname(__FILE__)),'external',"libbam.#{lib_os}")
        end #filename
        module_function :filename
      end #Library
    end #Sam
  end #DB
end #Bio