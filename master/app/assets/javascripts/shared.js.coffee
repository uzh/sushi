# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/

@dispLoading = (msg) ->
    dispMsg = ""
    if( msg != "" )
        dispMsg = "<div class='loadingMsg'>" + msg + "</div>"
    if($("#loading").size() == 0)
        $("body").append("<div id='loading'>" + dispMsg + "</div>")

@removeLoading = () ->
    $("#loading").remove()

