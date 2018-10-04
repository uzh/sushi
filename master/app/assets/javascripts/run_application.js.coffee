# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/

$ ->
    $("[with='alphanum']").on "keyup.inputcontrol.alphanum", ->
            match_part = $(this).val().match(/[^0-9a-zA-Z_]/)
            if match_part
              alert "WARNING: #{match_part} is not allowed"
              $(this).addClass('attention')
              $(this).val($(this).val().replace(/[^0-9a-zA-Z_]/g,""))
    $("[with='numeric']").on "keyup.inputcontrol.numeric", ->
            match_part = $(this).val().match(/[^0-9]/)
            if match_part
              alert "WARNING: #{match_part} is not allowed"
              $(this).addClass('attention')
              $(this).val($(this).val().replace(/[^0-9]/g,""))
    $("[with='ascii']").on "keyup.inputcontrol.ascii", ->
            match_part = $(this).val().match(/[^ -~]|[\\{}$%#!*]/)
            if match_part
              alert "WARNING: #{match_part} is not allowed"
              $(this).addClass('attention')
              $(this).val($(this).val().replace(/[^ -~]|[\\{}$%#!*]/g,""))
