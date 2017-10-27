# Place all the behaviors and hooks related to the matching controller here.
# All this logic will automatically be available in application.js.
# You can use CoffeeScript in this file: http://jashkenas.github.com/coffee-script/

$ ->
    $('#page_unit').change (event) ->
      idx = this.selectedIndex
      page_number = this.options[idx].text
      $('#page_unit').val(page_number)
      this.form.submit()
