module ApplicationHelper
  def linebreak_to_br(text)
    text.gsub(/\r\n|\r|\n/, "<br />")
  end
  def project_init
    if !session[:projects] or params[:select_project]
      @fgcz = SushiFabric::Application.config.fgcz?
      session[:employee] = (@fgcz and FGCZ.get_user_groups(current_user.login).include?('Employees'))
      session[:projects] = if @fgcz 
                             FGCZ.get_user_projects(current_user.login).map{|project| project.gsub(/p/,'').to_i}.sort
                           else
                             [1001]
                           end
      session[:project] = if @fgcz
                            if project=params[:select_project] and number=project[:number] and number.to_i!=0 or
                               project=params[:project] and number=project[:number] and number.to_i!=0 and 
                               (session[:employee] or session[:projects].include?(number.to_i))
                              current_user.selected_project = number
                              current_user.save
                              number.to_i
                            elsif project_id = params[:project_id]
                              number = project_id.gsub(/p/,'')
                              current_user.selected_project = number
                              current_user.save
                              number.to_i
                            elsif current_user.selected_project != -1
                              current_user.selected_project
                            else
                              session[:projects].first
                            end
                          else
                            session[:projects].first
                          end
      if @fgcz and current_user.selected_project == -1
        current_user.selected_project = session[:project]
        current_user.save
      end
    end
  end

end
