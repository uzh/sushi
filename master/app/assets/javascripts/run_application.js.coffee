# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/

$ ->
    $("[with='alphanum']").on "keyup.inputcontrol.alphanum", ->
            $(this).val($(this).val().replace(/[^0-9a-zA-Z_]/g,""))
    $("[with='numeric']").on "keyup.inputcontrol.numeric", ->
            $(this).val($(this).val().replace(/[^0-9]/g,""))
    $("[with='ascii']").on "keyup.inputcontrol.ascii", ->
            $(this).val($(this).val().replace(/[^ -~]/g,""))
