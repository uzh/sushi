# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/


$ ->
  $("#refresh_sushi_table").click ->
    alert("Updating table will take some seconds. Please do not click refresh button many times quickly.")
    $.get("/sushi_application/refresh_table")
