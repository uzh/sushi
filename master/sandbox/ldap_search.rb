#!/usr/local/ngseq/bin/ruby

$:.unshift File.dirname(__FILE__)

require 'fgcz'

groups = FGCZ.get_user_groups "trxcopy", "fXSC7o8g"
p groups
