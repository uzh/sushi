#!/usr/bin/env ruby
# encoding: utf-8

require 'optparse'
# Version = '20130517-120455'

class OptionParser
  attr_reader :p
  alias :_on :on
  def on(attr, *args, &block)
    if attr.is_a?(Symbol)
      self.class.class_eval do
        unless method_defined?(attr)
          attr_accessor attr 
        else
          raise "Method #{attr.to_s} is already defined in OptionParser class"
        end
      end
      unless args[0] =~ /\-/
        default = args.shift
        self.send(attr.to_s+"=", default)
      end
      _on(*args) do |i|
        self.send(attr.to_s+"=", i)
        block.call(i) if block
      end
    else
      args.unshift attr
      _on(*args, block)
    end
  end
end

if __FILE__ == $0
  opt = OptionParser.new do |o|
    o.banner = "Last update: #{o.version}\nUsage: ruby #{__FILE__} [options]"
    o.on(:size, 100, '-N N', '--pop_size', Integer, 'population size (default: 100)'){|i| p i}
    o.on(:seed, '-R R', '--rseed', 'random seed')
    o.on(:flag, '-f', 'flag')
    o.on('-M M', '--hoge', Integer, 'hoge'){|i| p i}
    o.parse!(ARGV)
  end

  print "opt.flag = "
  p opt.flag
  print "opt.size = "
  p opt.size
  print "opt.seed = "
  p opt.seed
  puts
  puts opt.help
end
