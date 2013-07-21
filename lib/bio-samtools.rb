require 'ffi'
require 'bio/db/sam'
require 'bio/db/pileup'
require 'bio/db/vcf'
if (defined?(Rails) && Rails::VERSION::MAJOR >= 3)
  require 'bio-samtools-engine'
else
  #require the usual dependencies here
  
end
