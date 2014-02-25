class RunApplicationController < ApplicationController
  class Prefecture
    def initialize
      @cities = []
    end
    attr_accessor :id, :name
    attr_accessor :cities
  end
  class City
    attr_accessor :id, :name
  end
  def init(pref_id)
    p1 = Prefecture.new
    p1.id = 1
    p1.name = 'Tokyo'
    p2 = Prefecture.new
    p2.id = 2
    p2.name = 'Saitama'
    p3 = Prefecture.new
    p3.id = 3
    p3.name = 'Kanagawa'
    @prefectures = [p1, p2, p3]

    c11 = City.new
    c11.id = 11
    c11.name = 'Asakusa'
    c12 = City.new
    c12.id = 12
    c12.name = 'Hachioji'
    p1.cities = [c11, c12]

    c21 = City.new
    c21.id = 21
    c21.name = 'Kawaguchi'
    c22 = City.new
    c22.id = 22
    c22.name = 'Tokorozawa'
    p2.cities = [c21, c22]

    c31 = City.new
    c31.id = 31
    c31.name = 'Yokohama'
    p3.cities = [c31]

    pref = @prefectures.find{|pref| pref.id == pref_id}
    @cities = pref.cities

    params[:pref_id] = pref_id
    params[:city_id] = pref.cities[0].id
  end
  def city_select
    init(params[:pref_id].to_i)
  end
  def result
    @pref_id = params[:pref_id]
    @city_id = params[:city_id]
  end
  def index
    unless @prefectures
      init(1)
    end
    @data_sets = if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
                   project.data_sets.reverse
                 else
                   []
                 end
  end
  def set_parameters
    unless @prefectures
      init(1)
    end
    class_name = params[:app]
    #require class_name
    @sushi_app = eval(class_name).new
    data_set_id = params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)
    @nodes = {
      'fgcz-c-045: cpu 64,mem 504 GB,scr  15T' => 'fgcz-c-045',
      'fgcz-c-046: cpu 64,mem 504 GB,scr  11T' => 'fgcz-c-046',
      'fgcz-c-047: cpu 32,mem   1 TB,scr  28T' => 'fgcz-c-047',
      'fgcz-c-048: cpu 48,mem 252 GB,scr 3.5T' => 'fgcz-c-048',
      'fgcz-c-049: cpu  8,mem  63 GB,scr 1.7T' => 'fgcz-c-049',
      'fgcz-c-051: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-051',
      'fgcz-c-052: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-052',
      'fgcz-c-053: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-053',
      'fgcz-c-054: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-054',
      'fgcz-c-055: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-055',
      'fgcz-c-057: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-057',
      'fgcz-c-058: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-058',
      'fgcz-c-059: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-059',
      'fgcz-c-061: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-061',
      'fgcz-c-063: cpu 12,mem  70 GB,scr 450G' => 'fgcz-c-063',
      'fgcz-c-065: cpu 24,mem  70 GB,scr 197G' => 'fgcz-c-065',
    }
  end
  def confirmation
    @params = params
    class_name = params[:sushi_app][:class]
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
                                 value.map{|v| v.chomp}.join(',')
                               elsif @sushi_app.params.data_type(key) == String
                                 value
                               else
                                 eval(value)
                               end
    end
  end
  def submit_jobs
    @params = params
    class_name = params[:sushi_app][:class]
    @sushi_app = eval(class_name).new
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
    if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
      @sushi_app.project = 'p' + project_number.to_s
    end
    @sushi_app.dataset_sushi_id = data_set_id.to_i
    @sushi_app.run
    @sushi_app.job_ids.each do |job_id|
      new_job = Job.new
      new_job.submit_job_id = job_id.to_i
      new_job.next_dataset_id = @sushi_app.next_dataset_id
      new_job.save
      new_job.data_set.jobs << new_job
      new_job.data_set.save
    end
  end
end
