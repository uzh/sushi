# SushiFabric

This library includes a core process called from [SUSHI](https://github.com/uzh/sushi). sushi_fabric command becomes available after installation. It requires [WorkflowManager](https://github.com/uzh/workflow_manager) to execute sushi_fabric command. SUSHI application must inherit the SushiFabric:SushiApp class and overwrite #next_dataset and #commands methods. Please refer to [SUSHI](https://github.com/uzh/sushi) for more details about how to make SUSHI application. 

## Installation

Add this line to your application's Gemfile:

    gem 'sushi_fabric'

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install sushi_fabric

## Usage

~~~~
  $ sushi_fabric -h                                                                                       16-05-10 15:37
  Usage: 
    sushi_fabric --class [SushiApp class] [options]
      (Either --dataset_id or --dataset_tsv is required)
      -c, --class class_name           SushiApp class name (required)
      -i, --dataset_id dataset_id      DataSet ID in Sushi DB
      -d, --dataset dataset_tsv        DataSet file (.tsv) (This option is prior to dataset_id option)
      -s, --dataset_name dataset_name  DataSet name in Sushi (This will be used with --dataset option, default: tsv file base name)
      -m parameterset_tsv,             Parameterset file (.tsv)
          --parameterset
      -r, --run                        Real run mode. without this option, it runs with test run mode which checks only DataSet and Parameters and no submittion
      -p, --project project            Project Number (default: 1001)
      -u, --user user                  Submit user (default: sushi_lover)
      -I, --load_path load_path        Add path where SushiApp class is located (default: ./lib)
      -n next_dataset_name,            Next DataSet Name (default: Analysis_Category+ID+Date )
          --next_dataset_name
~~~~

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request
