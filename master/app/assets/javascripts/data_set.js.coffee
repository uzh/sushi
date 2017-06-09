# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/

$ ->
    $(document).tooltip()
    $("#refresh_sushi_list").click (event) ->
        event.preventDefault()
        data_set_id = $("#data_set_id").val()
        dispLoading("updating...")
        $.ajax "/data_set/"+data_set_id+"/refresh_apps",
           type:"GET"
           complete:(data) ->
               removeLoading()

