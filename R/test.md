---
  title: "R Notebook"
output:
  html_document
---

  ## Example

  ```{r, results='asis'}
require(ggplot2)
for(i in 1:5){
  cat("### Section ", i, "\n")

  df <- mtcars[sample(5),]
  tb <- knitr::kable(df, caption = paste0("Table",i))
  g1 <- ggplot2::ggplot(df, aes(x = mpg, y = disp, fill = gear)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = paste0("Figure ", i))
  cat("\n")
  print(g1)
  print(tb)
  cat("\n")
}
