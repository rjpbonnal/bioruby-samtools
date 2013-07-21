class SlicesController < ApplicationController

  def index
    @index = "Something relevant"
  end

  def show
    @item = "It's me"
  end

  def new
    # @example =  Bio::Samtools::Example.new
  end

  def create
    # @example = Bio::Samtools::Example.new(params[:example])
    if @example.save
      redirect_to example_url(@example)
    else
      # This line overrides the default rendering behavior, which
      # would have been to render the "create" view.
      render :action => "new"
    end
  end

  def cut
    @p=params
  end
end