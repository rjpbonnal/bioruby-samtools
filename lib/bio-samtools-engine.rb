require 'bio-samtools'
require 'rails'

module BioSamtools
  class Engine < Rails::Engine

    # Config defaults
# DEPRECATED
    config.widget_factory_name = "default factory name"
    config.mount_at = '/'
    
    # # Load rake tasks
    # rake_tasks do
    #   load File.join(File.dirname(__FILE__), 'rails/railties/tasks.rake')
    # end
    
    # Check the gem config
    initializer "check config" do |app|
      # make sure mount_at ends with trailing slash
      config.mount_at += '/'  unless config.mount_at.last == '/'
    end
    
    initializer "static assets" do |app|
      app.middleware.use ::ActionDispatch::Static, "#{root}/public"
    end
    
    # consider the possibility to keep the modules in the lib directory, which is more compatible
    # with a normal gem/package
    paths["app/models"]           << "lib"
  end
end

# #In your Rails application into 
# # config/initializers
# # create a file called bio-samtools.rb
# # you can change the mount point of you engine, using option
# # --with_engine=namespace like http://YourDomain:port/YourMountPoint
# # Otherwise, if you want your RailsEngine to be called at root level, leave everything as it now, without that file.
# module BioSamtools
#    class Engine < Rails::Engine
#      engine_name = '/Samtools'
#    end
# end
