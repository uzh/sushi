#!/usr/bin/env ruby
# encoding: utf-8

module GlobalVariables
  SUSHI = 'Supercalifragilisticexpialidocious!!'


  def builder_selector(base_dir, shown_pattern=nil, value_pattern=nil)
    selector = {}
    Dir[base_dir].sort.select{|dir| File.directory?(dir)}.each do |dir|
      key = if shown_pattern
              dir.gsub(shown_pattern.keys.first, shown_pattern)
            else
              dir
            end
      value = if value_pattern
                dir.gsub(value_pattern.keys.first, value_pattern)
              else
                File.basename(dir)
              end
      selector[key] = value
    end
    selector
  end
  def ref_selector
    selector = {'select'=>''}

    base_pattern = "/srv/GT/reference/*/*/*"
    shown_replace_pattern = {/\/srv\/GT\/reference\//=>''}
    builds = builder_selector(base_pattern, shown_replace_pattern)

    base_pattern = "/srv/GT/reference/*/*/*/Annotation/Version*"
    value_replace_pattern = {/\/srv\/GT\/reference\/.+?\/.+?\//=>''}
    versions = builder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

    builds.keys.each do |build_key|
      if versions.keys.find{|version_key| version_key=~/#{build_key}/}
        builds.delete(build_key)
      end
    end
    builds.merge!(versions)
    builds = builds.sort.to_h

    selector.merge!(builds)
    selector
  end
  def factor_dataset
    factors = get_columns_with_tag 'Factor'
    dataset = {}
    factors.first.keys.each do |colname|
      dataset[colname+" [Factor]"] = @dataset[colname]
    end
    dataset
  rescue
    {}
  end
end
