#(c) Copyright 2011 Raoul Bonnal. All Rights Reserved. 

# create Rakefile for shared library compilation

require 'rbconfig'
require 'open-uri'

def self.lib_extension
  case Config::CONFIG['host_os']
    when /linux/ then return 'so'
    when /darwin/ then return 'dylib'
    when /mswin|mingw/ then raise NotImplementedError, "BWA library is not available for Windows platform"  
  end
end

Version = "0.1.16"
URL = "http://sourceforge.net/projects/samtools/files/samtools/#{Version}/samtools-#{Version}.tar.bz2/download"
SamToolsFile="samtools-#{Version}.tar.bz2"

File.open(File.join(path,"Rakefile"),"w") do |rakefile|
rakefile.write <<-RAKE
require 'rake/clean'

task :download => URL do
  open(URL) do |uri|
    File.open(SamToolsFile,'wb') do |fout|
      fout.write(uri.read)
    end #fout 
  end #uri
end
    
task :compile do
  sh "tar xvfj #{SamToolsFile}"
  sh "make"
  case Config::CONFIG['host_os']
    when /linux/ then "sh make libbam.so.1-local"
    when /darwin/ then "sh make libbam.1.dylib-local"
    when /mswin|mingw/ then raise NotImplementedError, "BWA library is not available for Windows platform"  
  end#case
end
  
task :default => [:dowload, :compile, :clean]
  
RAKE
  
end