module DataSetHelper
  # Minimal renderer tailored to methods.md's known structure (## / ### headings,
  # blank-line-separated paragraphs). Not a general Markdown parser.
  def render_methods_md(content)
    return "".html_safe if content.nil?

    html = []
    paragraph = []

    flush_paragraph = -> {
      unless paragraph.empty?
        html << "<p>#{paragraph.join('<br>')}</p>"
        paragraph = []
      end
    }

    content.each_line do |raw_line|
      line = raw_line.chomp
      if line =~ /^##\s+(.*)/
        flush_paragraph.call
        html << "<h2>#{ERB::Util.html_escape($1)}</h2>"
      elsif line =~ /^###\s+(.*)/
        flush_paragraph.call
        html << "<h3>#{ERB::Util.html_escape($1)}</h3>"
      elsif line.strip.empty?
        flush_paragraph.call
      else
        paragraph << ERB::Util.html_escape(line)
      end
    end
    flush_paragraph.call

    html.join.html_safe
  end

  def dataset_tree(indent, root)
    html = ""
    root.each do |node|
      html << "<tr><td>"
      html << "&nbsp;"*(indent*2)
      if indent > 0
        if node == root.last
          html << "&#9492;"
        else
          html << "&#9500;"
        end
      end
      html << node["text"]
      html << dataset_tree(indent+1, node["children"]) unless node["children"].empty?
      html << "</td></tr>"
    end
    html
  end
end
