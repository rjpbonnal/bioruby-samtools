Rails.application.routes.draw do

 #DEPRECATED 
 mount_at = BioSamtools::Engine.engine_name== "bio_samtools_engine" ? "" : BioSamtools::Engine.engine_name

#ROUTES from code below:
# YourPathToTheControllerFiles = bio/kb/gex
# OtherControllerName = slice
# 
#               gex     /gex(.:format)                        {:controller=>"bio/kb/gex/populations", :action=>"index"}
# filter_population GET /gex/populations/:id/filter(.:format) {:action=>"filter", :controller=>"bio/kb/gex/populations"}
#       populations GET /gex/populations(.:format)            {:action=>"index", :controller=>"bio/kb/gex/populations"}
#        population GET /gex/populations/:id(.:format)        {:action=>"show", :controller=>"bio/kb/gex/populations"}
#           samples GET /gex/samples(.:format)                {:action=>"index", :controller=>"bio/kb/gex/samples"}
#            sample GET /gex/samples/:id(.:format)            {:action=>"show", :controller=>"bio/kb/gex/samples"}
#            slices GET /gex/slices(.:format)                 {:action=>"index", :controller=>"bio/kb/gex/slices"}
#             slice GET /gex/slices/:id(.:format)             {:action=>"show", :controller=>"bio/kb/gex/slices"}
#  
# Uncomment to activate the resources
## or you can reconfigure mount_at 
## mount BioSamtools::Engine=>"your_mount_point", :as=>'your_personal_engine_name' # but is better to redefine the mount point in the Rails application.
yourPathToTheControllerFiles = 'samtools'
controllerName = :slices
# otherControllerName = :samples
## ---
scope mount_at, yourPathToTheControllerFiles do
  get '/slices/cut', to: 'slices#cut'	
  resources controllerName do #, :only => [ :index, :show, :new] do
    # member do
    # end
  end
  # resources otherControllerName, :only =>[:index, :show]
end #scope

 

end