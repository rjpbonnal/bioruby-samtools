#(c) Copyright 2011 Raoul Bonnal. All Rights Reserved. 

# create Rakefile for shared library compilation



path = File.expand_path(File.dirname(__FILE__))

path_external = File.join(path, "../lib/bio/db/sam/external")

version = File.open(File.join(path_external,"VERSION"),'r')
Version = version.read
version.close

url = "http://sourceforge.net/projects/samtools/files/samtools/#{Version}/samtools-#{Version}.tar.bz2/download"
SamToolsFile = "samtools-#{Version}.tar.bz2"
url2 = "https://github.com/samtools/bcftools/releases/download/#{Version}/bcftools-#{Version}.tar.bz2"
VcfToolsFile = "bcftools-#{Version}.tar.bz2"

File.open(File.join(path,"Rakefile"),"w") do |rakefile|
rakefile.write <<-RAKE
require 'rbconfig'
require 'open-uri'
require 'fileutils'
include FileUtils::Verbose
require 'rake/clean'

URL = "#{url}"
URL2 = "#{url2}"

task :download do
  open(URL) do |uri|
    File.open("#{SamToolsFile}",'wb') do |fout|
      fout.write(uri.read)
    end #fout 
  end #uri

    open(URL2) do |uri|
    File.open("#{VcfToolsFile}",'wb') do |fout|
      fout.write(uri.read)
    end #fout 
  end #uri
end
    
task :compile do
  sh "tar xvfj #{SamToolsFile}"
  sh "tar xvfj #{VcfToolsFile}"
  cd("samtools-#{Version}") do
    sh "make"
    cp("samtools", "#{path_external}")
  end #cd
  cd("bcftools-#{Version}") do
    sh "make"
    cp('bcftools', "#{path_external}")
  end
end
  
task :clean do
  cd("samtools-#{Version}") do
    sh "make clean"
  end
  cd("bcftools-#{Version}") do
    sh "make clean"
  end
  rm("#{SamToolsFile}")
  rm_rf("samtools-#{Version}")
  rm("#{VcfToolsFile}")
  rm_rf("bcftools-#{Version}")
end

task :default => [:download, :compile, :clean]
  
RAKE
  
end