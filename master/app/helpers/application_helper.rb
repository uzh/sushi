module ApplicationHelper
  def linebreak_to_br(text)
    text.gsub(/\r\n|\r|\n/, "<br />")
  end
  def employee?
    SushiFabric::Application.config.fgcz? and current_user and FGCZ.employee?(current_user.login)
  end
  def user_projects
    if SushiFabric::Application.config.fgcz? and current_user
      FGCZ.get_user_projects2(current_user.login).map{|project| project.gsub(/p/,'').to_i}.sort
    elsif SushiFabric::Application.config.course_mode and user_projects_ = SushiFabric::Application.config.course_users
      user_projects_.flatten.uniq.sort
    else
      [1001]
    end
  end
  def project_init
    if !session[:projects] or params[:select_project] or params[:project_id] or params[:id]
      @fgcz = SushiFabric::Application.config.fgcz?
      session[:employee] = employee?
      session[:projects] = user_projects
      session[:project] = if @fgcz and current_user
                            if params[:project_id].nil? and params[:id].nil? and project=params[:select_project] and number=project[:number] and number.to_i!=0 or
                               params[:project_id].nil? and params[:id].nil? and project=params[:project] and number=project[:number] and number.to_i!=0 and
                               (session[:employee] or session[:projects].include?(number.to_i))
                              # project text field or selection list event
                              current_user.selected_project = number
                              current_user.save
                              number.to_i
                            elsif project_id = params[:project_id] and number = project_id.gsub(/p/,'')
                              # direct link case with pXXX
                              if !session[:employee] and !session[:projects].include?(number.to_i)
                                number = session[:projects].first
                              end
                              current_user.selected_project = number
                              current_user.save
                              number.to_i
                            elsif id = params[:id] and data_set = DataSet.find_by_id(id) and number = data_set.project.number
                              # direct link case without pXXX
                              if !session[:employee] and !session[:projects].include?(number.to_i)
                                number = session[:projects].first
                              end
                              current_user.selected_project = number
                              current_user.save
                              number.to_i
                            elsif current_user.selected_project != -1 and session[:projects].include?(current_user.selected_project.to_i)
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
  def partitions
    command = "sinfo --format=%R"
    list = `#{command}`.split(/\n/)
    list.delete("PARTITION")
    if i = list.index("employee")
      list.delete("employee")
      list.unshift("employee")
    end
    list
  end
  
  # System Maintenance Announcement helper methods
  def maintenance_announcement_enabled?
    Rails.application.config.respond_to?(:maintenance_announcement_enabled) &&
      Rails.application.config.maintenance_announcement_enabled
  end
  
  def maintenance_announcement_type
    Rails.application.config.maintenance_announcement_type || 'info'
  end
  
  def maintenance_announcement_message
    Rails.application.config.maintenance_announcement_message || ''
  end
  
  def maintenance_announcement_alert_class
    case maintenance_announcement_type
    when 'info'
      'alert-info'
    when 'warning'
      'alert-warning'
    when 'danger'
      'alert-danger'
    else
      'alert-info'
    end
  end
  
  def maintenance_announcement_icon
    case maintenance_announcement_type
    when 'info'
      'fa-info-circle'
    when 'warning'
      'fa-exclamation-triangle'
    when 'danger'
      'fa-exclamation-circle'
    else
      'fa-info-circle'
    end
  end
  
  def maintenance_announcement_label
    case maintenance_announcement_type
    when 'info'
      'Notice'
    when 'warning'
      'Warning'
    when 'danger'
      'Important'
    else
      'Notice'
    end
  end
end
