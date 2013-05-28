make_interactive <- function(filename, script, toggle="toggle") {
  require(gridSVG)
  grobs <- grid.ls()
  
  idx <- grep("panel-", grobs$name)
  for (i in idx) { 
    grid.garnish(grobs$name[i],
                 onmouseover=paste("frame('",grobs$name[i+2], ".1')", sep=""),
                 onmouseout=paste("deframe('",grobs$name[i+2], ".1')", sep=""), 
                 onmousedown=paste(sprintf("%shigh(evt, '", toggle),grobs$name[i+2], ".1')", sep=""))
  }

  # use script on server to get locally executable javascript code
  # or use inline option
  grid.script(filename=script)
  gridToSVG(filename)
}
