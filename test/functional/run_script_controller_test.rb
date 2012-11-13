require 'test_helper'

class RunScriptControllerTest < ActionController::TestCase
  test "should get index" do
    get :index
    assert_response :success
  end

end
