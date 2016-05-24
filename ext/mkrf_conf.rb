#(c) Copyright 2011 Raoul Bonnal. All Rights Reserved.

# create Rakefile for shared library compilation



path = File.expand_path(File.dirname(__FILE__))

path_external = File.join(path, "../lib/bio/db/sam/external")

version = File.open(File.join(path_external,"VERSION"),'r')
Version = version.read
version.close

#url = "http://sourceforge.net/projects/samtools/files/samtools/#{Version}/samtools-#{Version}.tar.bz2/download"
url="https://github.com/samtools/samtools/releases/download/#{Version}/samtools-#{Version}.tar.bz2"
SamToolsFile = "samtools-#{Version}.tar.bz2"
url_bcftools="https://github.com/samtools/bcftools/releases/download/#{Version}/bcftools-#{Version}.tar.bz2"
BcfToolsFile = "bcftools-#{Version}.tar.bz2"

File.open(File.join(path,"Rakefile"),"w") do |rakefile|
rakefile.write <<-RAKE
require 'rbconfig'
require 'open-uri'
#require 'open_uri_redirections'
require 'fileutils'
include FileUtils::Verbose
require 'rake/clean'

URL = "#{url}"
URL_bcf = "#{url_bcftools}"
task :download do
  open(URL) do |uri|
    File.open("#{SamToolsFile}",'wb') do |fout|
      fout.write(uri.read)
    end #fout
  end #uri

  open(URL_bcf) do |uri|
    File.open("#{BcfToolsFile}",'wb') do |fout|
      fout.write(uri.read)
    end #fout
  end #uri

end

task :compile do
  sh "tar xvfj #{SamToolsFile}"
  cd("samtools-#{Version}") do
    sh "make"
    cp('samtools', "#{path_external}")
  end #cd

  sh "tar xvfj #{BcfToolsFile}"
  cd("bcftools-#{Version}") do
    sh "make"
    cp('bcftools', "#{path_external}")
  end #cd
end

task :clean do
  cd("samtools-#{Version}") do
    sh "make clean"
  end
  rm("#{SamToolsFile}")
  rm_rf("samtools-#{Version}")
  rm("#{BcfToolsFile}")
  rm_rf("bcftools-#{Version}")
end

task :default => [:download, :compile, :clean]

RAKE

end
