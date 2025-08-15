module ApplicationHelper
  # Lightweight Markdown-like renderer without external gems.
  # Supports:
  # - Bold (**text**), Italic (*text*), Inline code (`code`)
  # - Headings (1. , 2. or # , ## are converted to <hN>)
  # - Unordered list (- , * ) and ordered list (1. 2. ...)
  # - Paragraphs and line breaks
  def render_markdown(text)
    return ''.html_safe if text.to_s.strip.empty?

    safe = h(text)

    # Inline formatting
    safe.gsub!(/\*\*(.+?)\*\*/, '<strong>\1</strong>')
    safe.gsub!(/(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)/, '<em>\1</em>')
    safe.gsub!(/`([^`]+)`/, '<code>\1</code>')

    lines = safe.split(/\r?\n/)
    html_lines = []
    in_ul = false
    in_ol = false

    lines.each do |line|
      if line.strip.empty?
        if in_ul
          html_lines << '</ul>'
          in_ul = false
        end
        if in_ol
          html_lines << '</ol>'
          in_ol = false
        end
        html_lines << '<br />'
        next
      end

      # Headings like "1. Title:" or markdown-style '# Title'
      if line =~ /^\s*#+\s+(.*)$/
        level = [line[/^\s*(#+)/, 1].length, 6].min
        html_lines << "<h#{level}>#{Regexp.last_match(1).strip}</h#{level}>"
        next
      elsif line =~ /^\s*(\d+)\.\s+(.*)$/
        # Ordered list detection will handle list items, but for a single leading number followed by bold title like '1. **Job Summary:**'
        # treat as paragraph if it contains bold colon title
        content = Regexp.last_match(2)
        if content =~ /<strong>.*:<\/strong>/
          html_lines << "<p>#{Regexp.last_match(1)}. #{content}</p>"
          next
        end
      end

      # List items
      if line =~ /^\s*[-\*]\s+(.*)$/
        html_lines << '<ul>' unless in_ul
        in_ul = true
        html_lines << "<li>#{Regexp.last_match(1)}</li>"
        next
      elsif line =~ /^\s*\d+\.\s+(.*)$/
        html_lines << '<ol>' unless in_ol
        in_ol = true
        html_lines << "<li>#{Regexp.last_match(1)}</li>"
        next
      end

      # Paragraph
      html_lines << "<p>#{line.strip}</p>"
    end

    html_lines << '</ul>' if in_ul
    html_lines << '</ol>' if in_ol

    allowed_tags = %w(p br b i em strong a ul ol li code pre h1 h2 h3 h4 h5 h6 blockquote hr)
    allowed_attrs = %w(href title)
    sanitize(html_lines.join("\n"), tags: allowed_tags, attributes: allowed_attrs)
  end
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
end
