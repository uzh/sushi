module DataSetHelper
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
