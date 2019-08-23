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
  def make_fgcz_node_list
    node2scr = {}
    command = "qhost -F scratch"
    keep = nil
    IO.popen(command) do |out|
      while line=out.gets
        hostname, arch, ncpu, loading, memtot, memuse, *others = line.split
        if hostname =~ /fgcz/
          keep = hostname
        elsif scratch_ = line.chomp.split.last and
              scratch = scratch_.split('=').last
          node2scr[keep] = scratch.to_i
          keep = nil
        end
      end
    end

    node_list = {}
    keep = nil
    command = 'qhost -q'
    IO.popen(command) do |out|
      while line=out.gets
        # HOSTNAME                ARCH         NCPU  LOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS
        hostname, arch, ncpu, loading, memtot, memuse, *others = line.split
        if hostname =~ /fgcz/
          #puts [hostname, ncpu, loading, memtot, memuse].join("\t")
          mem = memtot.gsub(/G/, '').to_i
          keep = [hostname, ncpu, "#{mem}G"]
        elsif hostname == "GT" and keep and cores = line.chomp.split.last and cores !~ /[du]/
          hostname = keep.shift
          keep[0] = cores
          if scr = node2scr[hostname] and scr >= 1000
            scr = "%.1f" % (scr.to_f / 1000)
            scr << "T"
          else
            scr = scr.to_s + "G"
          end
          keep << scr
          node_list[hostname] = keep 
          keep = nil
        end
      end
    end

    # reformat
    nodes = {}
    node_list.each do |hostname, specs|
      cores, ram, scr = specs
      key = "#{hostname}: cores #{cores}, ram #{ram}, scr #{scr}"
      value = hostname
      nodes[key] = value
    end
    nodes
  end
  def set_parameters
    class_name = params[:app]
    require class_name
    @sushi_app = eval(class_name).new
    @sushi_app.workflow_manager = @@workflow_manager
    resubmit_data_set_id = nil
    data_set_id = if data_set = params[:data_set] # usual case
                    data_set[:id]
                  elsif id = params[:data_set_id] and resubmit_data_set_id = params[:resubmit_data_set_id] # resubmit case
                    id
                  end
    @data_set = DataSet.find_by_id(data_set_id.to_i)
    @sushi_app.dataset_sushi_id = data_set_id.to_i
    @sushi_app.set_input_dataset
    @sushi_app.set_default_parameters
    @params_selected = {}
    if resubmit_data_set_id
      # load parameters
      parent_data_set_path = sample_path(@data_set)
      resubmit_data_set = DataSet.find_by_id(resubmit_data_set_id)
      resubmit_data_set_path = sample_path(resubmit_data_set)
      parameters_path = File.join(SushiFabric::GSTORE_DIR, resubmit_data_set_path)
      parameters_tsv = File.join(parameters_path, "parameters.tsv")
      parameterset_tsv = CSV.readlines(parameters_tsv, :col_sep=>"\t")
      parameterset_tsv.each do |row|
        header, value = row
        unless header == "sushi_app"
          @params_selected[header] = if @sushi_app.params.data_type(header) == String or value == nil or @sushi_app.params.data_type(header) == NilClass
                                       value
                                     else
                                       eval(value)
                                     end
          if !@sushi_app.params[header].instance_of?(Array) and
             !@sushi_app.params[header].instance_of?(Hash) and
             @sushi_app.params.data_type(header) != NilClass
            @sushi_app.params[header] = @params_selected[header]
          end
        end
      end
    end
    @nodes = if SushiFabric::Application.config.fgcz?
               make_fgcz_node_list
             else
               @sushi_app.cluster_nodes
             end
    if SushiFabric::Application.config.fgcz?
      if session['employee']
        @nodes = @nodes.select{|node| node =~ /fgcz-h-00[89]/} if SushiFabric::Application.config.course_mode
      elsif !current_user # demo sushi
        @nodes = @nodes.select{|node| node =~ /fgcz-h-00[89]/} unless SushiFabric::Application.config.course_mode
      else # normal user in prod sushi
        @nodes = if SushiFabric::Application.config.course_mode
                   @nodes.select{|node| node =~ /fgcz-h-00[89]/}
                 else
                   @nodes.select{|node| node =~ /fgcz-h/ or node =~ /fgcz-c-048/}
                 end
      end
      if current_user and FGCZ.get_user_projects(current_user.login).include?('p1535')
        @nodes['fgcz-c-047: cpu 32,mem   1 TB,scr  28T'] = 'fgcz-c-047'
      end

      # 20190823 masa tentatively in 2-3 months for Gwyneth
      @nodes.delete('fgcz-c-047: cores 8/16, ram 1009G, scr 900G')
      @nodes
    end
		unless @factors
			#init_factor('Condition')
			init_factor
		end

    @process_mode = @sushi_app.params['process_mode']
    @samples = @data_set.samples.map{|sample|
      name = sample.to_hash["Name"]
      [name, name]
    }
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
      @sushi_app.params[key] = if key == 'node' or key == 'samples' or @sushi_app.params[key, "multi_selection"]
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
    @sushi_app.params.each do |key, value|
      if @sushi_app.params[key, "multi_selection"] and !params[:parameters][key]
        @sushi_app.params[key] = ''
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
    @data_set_id = params[:data_set][:id]
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
    if current_data_set = DataSet.find_by_id(@data_set_id.to_i) and
       project = current_data_set.project and
       project_number = project.number
      @sushi_app.project = 'p' + project_number.to_s
    end
    @sushi_app.dataset_sushi_id = @data_set_id.to_i
    @sushi_app.current_user = current_user
    @sushi_app.off_bfabric_registration = session[:off_bfabric_registration]
    submit_type = params[:submit_type]
    if submit_type == "Submit"
      @sushi_app.run
    elsif submit_type == "MockRun"
      @sushi_app.mock_run
    end
  end
end
