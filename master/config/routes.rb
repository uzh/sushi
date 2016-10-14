SushiFabric::Application.routes.draw do
  get "sushi_application/index"

  root :to => "home#index"
  post "/" => "home#index"
  
  devise_for :users
  
  resources :job_monitoring do
    member do
      post :kill_job
      post :change_status
    end
    collection do
      post :index
      get :print_log
      get :print_script
    end
  end
  
  #resources :data_set do
  resources :data_set, :only => [:index, :show, :destroy] do
    member do
      post :edit
      post :save_as_tsv
      get :refresh_apps
      post :confirm_delete_only_data_files
      post :run_delete_only_data_files
    end
    collection do
      get :whole_treeviews
      get :partial_treeviews
      post :import
      post :delete
      get :script_log
      get :job_parameter
      post :multi_delete
      post :multi_destroy
      get :index_full
      post :report
      post :bfabric
    end
  end
  # aliases
  #get "/data_set/:project_id/index" => "data_set#index"
  #get "/data_set/:project_id/index_full" => "data_set#index_full"
  get "/data_set/:project_id/:id(.:format)" => "data_set#show"

  resources :sample, :only => [] do
    member do
      post :show
      post :edit
      post :multiedit
    end
  end

  resources :sushi_application, :only => [:index] do
    collection do
      post :refresh
      get :refresh_table
    end
  end

  get "/api/:method" => "api#index"
  post "/api/:method" => "api#index"
  
  resources :run_application, :only => [:index] do
    collection do 
      post :set_parameters
      post :confirmation
      post :submit_jobs
			get  :factor_select
			post  :factor_result
    end
  end

  #get "/projects/:project_id" => redirect("http://fgcz-gstore.uzh.ch/projects/%{project_id}")
  get "/projects/:project_id(/*dirs)" => "home#gstore"
  get "/check_sushi_constants" => "home#sushi_constants"
  get "/import/*dataset" => "data_set#import_from_gstore"
  get "/sushi_rank" => "home#sushi_rank"


#	match "/city_select" => "run_application#city_select"
#	match "/result" => "run_application#result"

  # The priority is based upon order of creation:
  # first created -> highest priority.

  # Sample of regular route:
  #   match 'products/:id' => 'catalog#view'
  # Keep in mind you can assign values other than :controller and :action

  # Sample of named route:
  #   match 'products/:id/purchase' => 'catalog#purchase', :as => :purchase
  # This route can be invoked with purchase_url(:id => product.id)

  # Sample resource route (maps HTTP verbs to controller actions automatically):
  #   resources :products

  # Sample resource route with options:
  #   resources :products do
  #     member do
  #       get 'short'
  #       post 'toggle'
  #     end
  #
  #     collection do
  #       get 'sold'
  #     end
  #   end

  # Sample resource route with sub-resources:
  #   resources :products do
  #     resources :comments, :sales
  #     resource :seller
  #   end

  # Sample resource route with more complex sub-resources
  #   resources :products do
  #     resources :comments
  #     resources :sales do
  #       get 'recent', :on => :collection
  #     end
  #   end

  # Sample resource route within a namespace:
  #   namespace :admin do
  #     # Directs /admin/products/* to Admin::ProductsController
  #     # (app/controllers/admin/products_controller.rb)
  #     resources :products
  #   end

  # You can have the root of your site routed with "root"
  # just remember to delete public/index.html.
  # root :to => 'welcome#index'

  # See how all your routes lay out with "rake routes"

  # This is a legacy wild controller route that's not recommended for RESTful applications.
  # Note: This route will make all actions in every controller accessible via GET requests.
  # match ':controller(/:action(/:id))(.:format)'
end
