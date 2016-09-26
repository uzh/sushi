# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/

$ ->
    $("#refresh_sushi_list").click ->
        data_set_id = $("#data_set_id").val()
        dispLoading("updating...")
        $.ajax "/data_set/"+data_set_id+"/refresh_apps",
           type:"GET"
           complete:(data) ->
               removeLoading()

dispLoading = (msg) ->
    dispMsg = ""
    if( msg != "" )
        dispMsg = "<div class='loadingMsg'>" + msg + "</div>"
    if($("#loading").size() == 0)
        $("body").append("<div id='loading'>" + dispMsg + "</div>")

removeLoading = () ->
    $("#loading").remove()

