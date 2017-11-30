source "http://rubygems.org"
# Add dependencies required to use your gem here.
# Example:
gem "bio-svgenes", ">= 0.4.1"
gem "bio", ">= 1.4.2"
gem "rake" 
#gem 'open_uri_redirections'

# Add dependencies to develop your gem here.
# Include everything needed to run rake, tests, features, etc.
group :development do
  gem "shoulda", "= 2.10"
  gem 'test-unit'
if RUBY_VERSION.start_with?("2.1") or RUBY_VERSION.start_with?("2.2") or RUBY_VERSION.start_with?("2.0")
  gem "jeweler", "= 2.0.1"
else
	gem "juwelier" ,  :platforms => :ruby_23 #jeweler support is being dropped
end
  gem "rack", "1.6.4",  :platforms => :ruby_21


end
