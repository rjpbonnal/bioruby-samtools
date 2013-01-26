module Bio
  class DB
    module SAM
      module Library
        #IMPORTANT NOTE: Windows library is missing in this distribution

        # Return the path with the file name of the library for the specific operating system
        def filename
          #TODO refactor this piece of code in all the files
          lib_os = case RUBY_PLATFORM
          when /linux/
            'so.1'
          when /darwin/
            '1.dylib'
          when /windows/
            'dll'
          else
            case RUBY_DESCRIPTION
            when /jruby.*darwin/
              '1.dylib'
            when /jruby.*linux/
            'so.1'
            end
          end

          File.join(File.expand_path(File.dirname(__FILE__)),'external',"libbam.#{lib_os}")
        end #filename
        module_function :filename
      end #Library
    end #Sam
  end #DB
end #Bio