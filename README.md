# bio-samtools

The original project samtools-ruby belongs to Ricardo H. Ramirez @ [https://github.com/homonecloco/samtools-ruby] (https://github.com/homonecloco/samtools-ruby)

## Introduction

Documentation and code come from that project and we'll adapt it for a better integration in BioRuby.

Binder of samtools for ruby, on the top of FFI. 

This project was born from the need to add support of BAM files to 
the [gee_fu genome browser] (http://github.com/danmaclean/gee_fu). 

## Installation

Add this line to your application's Gemfile:

    gem 'bio-samtools'
    
And then execute:
 
    bundle
    
Or install it yourself as:

    $ gem install bio-samtools

## Usage

### Creating a new SAM object

A SAM object represents the alignments in the BAM file, and is very straightforward to create, you will need a sorted BAM file, to access the alignments and a reference sequence in FASTA format to use the reference sequence. The object can be created and opened as follows:

	require 'bio-samtools'
	
	bam = Bio::DB::Sam.new(:bam=>"my_sorted.bam", :fasta=>'ref.fasta')
	bam.open

### Getting Reference Sequence

Retrieving the reference can only be done if the reference has been loaded, which isn't done automatically in order to save memory. Reference need only be loaded once, and is accessed using reference name, start, end in 1-based co-ordinates. A standard Ruby String object is returned.

	bam.load_reference 
	sequence_fragment = bam.fetch_reference("Chr1", 1, 500)
	
### Getting Alignments

Alignments can be obtained one at a time by looping over a specified region using the fetch() function.

	bam.load_reference 
	bam.fetch("1",3000,4000).each do |alignment|
		#do something with the alignment...
	end

See more detail on doc/tutorial.html or doc/tutorial.pdf for a walkthrough tutorial. 

## Dependencies

* BioRuby >= 1.4.3 [https://github.com/bioruby/bioruby](https://github.com/bioruby/bioruby)
* Ruby 1.9 and above. 
	
## FAQ
* I´m getting a **segmentation Fault**, what did I do wrong?

	[Answer] You are using an old version of bio-samtools, the current version doesn't link directly to the library. 

* I keep seeing this **Invalid gemspec in [some ruby gem path…]**, what is wrong?
	
	[Answer] This appears to be a bug in RubyGems that doesn't affect the running of the tools. It will keep happening until someone updates RubyGems. If it really bugs you, downgrade RubyGems.

## Contributing to bio-samtools
 
* Check out the latest master to make sure the feature hasn't been implemented or the bug hasn't been fixed yet
* Check out the issue tracker to make sure someone already hasn't requested it and/or contributed it
* Fork the project
* Start a feature/bugfix branch
* Commit and push until you are happy with your contribution
* Make sure to add tests for it. This is important so I don't break it in a future version unintentionally.
* Please try not to mess with the Rakefile, version, or history. If you want to have your own version, or is otherwise necessary, that is fine, but please isolate to its own commit so I can cherry-pick around it.

### TODO
1. Filter to the fetching algorithm (give a condition that has to be satisfied to add the alignment to the list)

### To whom do I complain?
Try [Ricardo.Ramirez-Gonzalez@tgac.ac.uk](Ricardo.Ramirez-Gonzalez@tgac.ac.uk)
 and [dan.maclean@tsl.ac.uk](dan.maclean@tsl.ac.uk)

### Important Notes
* Libraries (libbam) are downloaded, compiled and installed inside the gem at install time on the host system

    `openssl dgst libbam.so.1` MD5 is c45cfccfb41ffeb2730ee4b227d244c4

### Important Notes for developers

Remember that you must compile and install the right libbam library for you host system. In order to do that there are three possible solutions:

* download, compile and install the library in bioruby-samtools-your_clone/lib/bio/db/sam/external/libbam.xxxxx by yourself
* install the gem and then grab the compiled library `cp 'locate libbam.1.dylib' bioruby-samtools-your_clone/lib/bio/db/sam/external` (library name is an example)
* in your bioruby-samtools-your_clone create the Rakefile typing `cd ext; ruby mkrf_conf.rb; rake -f Rakefile`

The latest I think is the easiest way, cause you are replicating the automatic process.

For testing just run `rake test`. Tests must be improved.

## Copyright

Copyright (c) 2011 Raoul J.P. Bonnal. See LICENSE.txt for
further details.

