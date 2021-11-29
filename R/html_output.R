#' @import xml2

html_output <- function(html_dir, config, replace = 3, prefix = "./", label = "") {
    html <- xml2::read_html(html_dir)
    xml2::xml_set_text(xml2::xml_child(html, replace), config)
    xml2::write_html(html, file.path(prefix, paste(label, "summary.html", sep = "")))
}