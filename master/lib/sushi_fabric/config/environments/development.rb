SushiFabric::Application.configure do
  # sushi_fabric
  config.workflow_manager = 'druby://localhost:12345'
  config.gstore_dir = File.join(Dir.pwd, 'public/gstore/projects')
  config.sushi_app_dir = Dir.pwd
  config.scratch_dir = '/tmp/scratch'
end
