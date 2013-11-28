#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'

$global_binding = Kernel.binding
class SushiWrap
  def initialize(job_script)
    @template =<<-EOS
class CLASS_NAME < SushiFabric::SushiApp
  def initialize
    super
INIT
  end
  def next_dataset
NEXT_DATASET
  end
  def commands
COMMANDS
  end
end
    EOS
    @job_script = job_script
    @class_name = File.basename(job_script).gsub(/\.sh/,'')
    @init = ''
    @next_dataset = ''
    @commands = ''
  end
  def generate_class
    File.readlines(@job_script).each do |line|
      if line =~ /^\#\@/ 
        if line =~ /next_dataset/
          dataset = eval(line.gsub(/\#\@/,''))
          dataset = Hash[*dataset]
          dataset.keys.each do |key|
            value = dataset[key]
            value1, value2 = value.split('.')
            if @required_columns.include?(value1)
              dataset[key] = if key.include?('File') 
                               "File.join(@result_dir, @dataset['" + value1 + "'].to_s + '." + value2.to_s + "')"
                             else
                               "@dataset['" + value1 + "']" 
                             end
            end
          end
          @next_dataset = " "*4 + "{" + dataset.map{|key, value| "'" + key + "'=>" + value}.join(',') + "}"
        else
          eval(line.gsub(/\#/,''))
          @init << " "*4 + line.gsub(/\#/,'') 
        end
      elsif line !~ /\#/ and !line.chomp.empty?
        values = line.scan(/\$\w+/)
        replaces = {}
        values.each do |val|
          if @required_columns.include?(val.gsub(/\$/,''))
            replaces[val] = '#{@dataset[\'' + val.gsub(/\$/,'') + "']}"
          else
            replaces[val] = val
          end
        end
        @commands << " "*4 + '"' + line.gsub(/\$\w+/, replaces).chomp + '"' + "\n"
      end
    end

    @template.gsub(/CLASS_NAME/, @class_name).gsub(/INIT/,@init).gsub(/NEXT_DATASET/, @next_dataset).gsub(/COMMANDS/,@commands)
  end
  def define_class
    eval(generate_class, $global_binding)
  end
  def debug
    print <<-EOS
#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'

    EOS
    print generate_class
    print <<-EOS
if __FILE__ == $0
  usecase = #{@job_script.gsub(/\.sh/,'')}.new

  usecase.project = "p1001"
  usecase.user = 'sushi_lover'
  usecase.parameterset_tsv_file = 'sample_parameterset.tsv'
  usecase.dataset_tsv_file = 'sample_dataset.tsv'

  # run (submit to workflow_manager)
  #usecase.run
  usecase.test_run
end
    EOS
  end
end

if __FILE__ == $0
  unless job_script=ARGV[0] and run_or_debug=ARGV[1]
    puts "Usage:\n ruby #{__FILE__} [job_script.sh] [run|test|debug]"
    exit
  end
  sushi_wrap = SushiWrap.new(job_script)

  if run_or_debug == 'run' or run_or_debug == 'test'
    sushi_wrap.define_class
    class_name = job_script.gsub(/\.sh/,'')
    usecase = eval(class_name).new

    usecase.project = "p1001"
    usecase.user = 'sushi_lover'
    usecase.parameterset_tsv_file = 'sample_parameterset.tsv'
    usecase.dataset_tsv_file = 'sample_dataset.tsv'

    # run (submit to workflow_manager)
    if run_or_debug == 'run'
      usecase.run
    else
      usecase.test_run
    end
  else
    sushi_wrap.debug 
  end

end
