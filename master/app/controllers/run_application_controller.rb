class RunApplicationController < ApplicationController
	def init_factor(factor_key=nil)
		@factor_colums = {}
    data_set_id = params[:data_set_id]||params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)
    if @data_set
			@data_set.samples.each do |sample|
				sample.to_hash.each do |header, value|
					if header.tag?('Factor') 
						key = header.split(/\[/).first.strip
						#@factor_colums[header] ||= []
						#@factor_colums[header] << value

						@factor_colums[key] ||= []
						@factor_colums[key] << value
					end
				end
			end
			@factor_colums.keys.each do |header|
				@factor_colums[header].uniq!
			end
    end
    unless @factor_colums.empty?
      factor_key = @factor_colums.keys.first unless factor_key
      @factors = @factor_colums[factor_key]
      params[:grouping] = factor_key
      params[:sampleGroup] = @factor_colums[params[:grouping]]
      params[:refGroup] = @factor_colums[params[:grouping]]
      @option_list = @factors.map{|v| [v,v]}
      @option_list.unshift(["please select", ''])
    end
	end
	def factor_select
		init_factor(params[:grouping])
	end
  def index
    @data_sets = if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
                   project.data_sets.reverse
                 else
                   []
                 end
  end
  def set_parameters
    class_name = params[:app]
    require class_name
    @sushi_app = eval(class_name).new
    @sushi_app.workflow_manager = @@workflow_manager
    data_set_id = params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)

    @sushi_app.dataset_sushi_id = data_set_id.to_i
    @sushi_app.set_input_dataset
    @sushi_app.set_default_parameters
    @nodes = @sushi_app.cluster_nodes
    if SushiFabric::Application.config.fgcz? and !session['employee']
      # Comment-out the next line if you want to activate all nodes in a course
      @nodes = @nodes.select{|node| node =~ /fgcz-h/ or node =~ /fgcz-c-065/}
      if current_user and FGCZ.get_user_projects(current_user.login).include?('p1535')
        @nodes['fgcz-c-047: cpu 32,mem   1 TB,scr  28T'] = 'fgcz-c-047'
      end
    end
		unless @factors
			#init_factor('Condition')
			init_factor
		end
  end
  def confirmation
    @params = params
    class_name = params[:sushi_app][:class]
    require class_name
    @sushi_app = eval(class_name).new
    data_set_id = params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)
    if next_dataset = params[:next_dataset] 
      if name = next_dataset[:name] and !name.to_s.strip.empty?
        @sushi_app.next_dataset_name = name.to_s.strip.gsub(/\s/,'_')
      end
      if comment = next_dataset[:comment] and !comment.to_s.strip.empty?
        @sushi_app.next_dataset_comment = comment.to_s.strip
      end
    end
    params[:parameters].each do |key, value|
      @sushi_app.params[key] = if key == 'node'
                                 if value.instance_of?(Array)
                                   value.map{|v| v.chomp}.join(',')
                                 else
                                   value
                                 end
                               elsif @sushi_app.params.data_type(key) == String
                                 value
                               else
                                 eval(value)
                               end
    end
    if @sushi_app.params['node'].empty?
      @sushi_app.workflow_manager = @@workflow_manager
      @sushi_app.params['node'] = @sushi_app.default_node
    end
    @sushi_app.params.each do |key, value|
      if @sushi_app.required_params and @sushi_app.required_params.include?(key) and value.to_s.empty? 
        @requires ||= {}
        @requires[key] = true 
      end
    end

  end
  def submit_jobs
    @params = params
    class_name = params[:sushi_app][:class]
    require class_name
    @sushi_app = eval(class_name).new
    @sushi_app.logger = logger
    @sushi_app.workflow_manager = @@workflow_manager
    @sushi_app.user = if current_user 
                        current_user.login
                      else
                        'sushi_lover'
                      end
    data_set_id = params[:data_set][:id]
    if next_dataset = params[:next_dataset] 
      if name = next_dataset[:name] and !name.to_s.strip.empty?
        @sushi_app.next_dataset_name = name.to_s.strip.gsub(/\s/,'_')
      end
      if comment = next_dataset[:comment] and !comment.to_s.strip.empty?
        @sushi_app.next_dataset_comment = comment.to_s.strip
      end
    end
    params[:parameters].each do |key, value|
      @sushi_app.params[key] = if @sushi_app.params.data_type(key) == String
                                       value
                                     else
                                       eval(value)
                                     end
    end
    if current_data_set = DataSet.find_by_id(data_set_id.to_i) and
       project = current_data_set.project and
       project_number = project.number
      @sushi_app.project = 'p' + project_number.to_s
    end
    @sushi_app.dataset_sushi_id = data_set_id.to_i
    @sushi_app.current_user = current_user
    @sushi_app.off_bfabric_registration = session[:off_bfabric_registration]
    @sushi_app.run
  end
end
