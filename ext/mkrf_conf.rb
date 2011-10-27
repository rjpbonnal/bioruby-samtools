#(c) Copyright 2011 Raoul Bonnal. All Rights Reserved. 

# create Rakefile for shared library compilation



path = File.expand_path(File.dirname(__FILE__))

path_external = File.join(path, "../lib/bio/db/sam/external")

version = File.open(File.join(path_external,"VERSION"),'r')
Version = version.read
version.close

url = "http://sourceforge.net/projects/samtools/files/samtools/#{Version}/samtools-#{Version}.tar.bz2/download"
SamToolsFile = "samtools-#{Version}.tar.bz2"

File.open(File.join(path,"Rakefile"),"w") do |rakefile|
rakefile.write <<-RAKE
require 'rbconfig'
require 'open-uri'
require 'fileutils'
include FileUtils::Verbose
require 'rake/clean'

URL = "#{url}"

task :download do
  open(URL) do |uri|
    File.open("#{SamToolsFile}",'wb') do |fout|
      fout.write(uri.read)
    end #fout 
  end #uri
end
    
task :compile do
  sh "tar xvfj #{SamToolsFile}"
  cd("samtools-#{Version}") do
    sh "patch < ../Makefile-bioruby.patch"
    case Config::CONFIG['host_os']
      when /linux/
        #sh "CFLAGS='-g -Wall -O2 -fPIC' make -e"
        sh "make"
        cp("libbam.a","#{path_external}")
        #sh "CFLAGS='-g -Wall -O2 -fPIC' make -e libbam.so.1-local"
        sh "make libbam.so.1-local"
        cp("libbam.so.1","#{path_external}")
      when /darwin/
        sh "make"
        cp("libbam.a","#{path_external}")
        sh "make libbam.1.dylib-local"
        cp("libbam.1.dylib","#{path_external}")      
      when /mswin|mingw/ then raise NotImplementedError, "BWA library is not available for Windows platform"  
    end #case
  end #cd
end
  
task :clean do
  cd("samtools-#{Version}") do
    sh "make clean"
  end
  rm("#{SamToolsFile}")
  rm_rf("samtools-#{Version}")
end

task :default => [:download, :compile, :clean]
  
RAKE
  
end