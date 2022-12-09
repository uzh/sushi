module ApplicationHelper
  def linebreak_to_br(text)
    text.gsub(/\r\n|\r|\n/, "<br />")
  end
  def project_init
    if !session[:projects] or params[:select_project]
      @fgcz = SushiFabric::Application.config.fgcz?
      session[:employee] = (@fgcz and current_user and FGCZ.get_user_groups(current_user.login).include?('Employees'))
      session[:projects] = if @fgcz and current_user
                             FGCZ.get_user_projects2(current_user.login).map{|project| project.gsub(/p/,'').to_i}.sort
                           elsif defined?(SushiFabric::Application.config.course_users) and user_projects = SushiFabric::Application.config.course_users
                             user_projects
                           else
                             [1001]
                           end
      session[:project] = if @fgcz and current_user
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
                            if project=params[:select_project] and number=project[:number] and number.to_i!=0 or
                               project=params[:project] and number=project[:number] and number.to_i!=0 and 
                               session[:projects].flatten.include?(number.to_i)
                               number.to_i
                            elsif session[:project]
                               session[:project]
                            else
                               session[:projects].first
                            end
                          end

      session[:partition] = if SushiFabric::Application.config.course_mode
                              "course"
                            elsif session[:employee]
                              "employee"
                            else
                              "user"
                            end

      if @fgcz and current_user and current_user.selected_project == -1
        current_user.selected_project = session[:project]
        current_user.save
      end
    end
  end
  def td(str)
    if str.to_s.length > 16
      str="<span title='"+str+"'>"+str.to_s.split(//)[0,16].join+"...</span>"
    end
    str.to_s.html_safe
  end
  def tooltip(str, help)
    str="<span title='"+help+"'>"+str+"</span>"
    str.to_s.html_safe
  end
  def remove_tag(header)
    header.gsub(/\s*\[.+\]/, '')
  end
end
