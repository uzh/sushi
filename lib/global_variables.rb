#!/usr/bin/env ruby
# encoding: utf-8

module GlobalVariables
  SUSHI = 'Supercalifragilisticexpialidocious!!'

  def build_selector
    selector = {'select'=>''}

    base_pattern = "/srv/GT/reference/*/*/*"
    shown_replace_pattern = {/\/srv\/GT\/reference\//=>''}
    builds = builder_selector(base_pattern, shown_replace_pattern)

    base_pattern = "/srv/GT/reference/*/*/*/Annotation/Version*"
    versions = builder_selector(base_pattern, shown_replace_pattern)

    builds.merge!(versions)
    builds = builds.sort.to_h

    selector.merge!(builds)
    selector
  end
end
