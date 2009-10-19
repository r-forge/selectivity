#
# sound stuff
#

play.sound <-
function (sound.file)
{
  if (missing(sound.file))
  {
    if (!exists(".__default.sound.file")) default.sound()
    sound.file <- .__default.sound.file
  }
  
  if (file.exists(sound.file))
  {
    Command <- paste("sndrec32 /play /close \"", sound.file ,"\"", sep="")
    try( system(Command, invisible=TRUE, wait=FALSE, show.output.on.console=FALSE) )
  } else
  {
    warning(paste("file: \"",sound.file,"\" does not exist", sep=""))
  }

  return ( invisible(NULL) )
}

set.sound <-
function ()
{
  if (exists(".__default.sound.file") && file.exists(.__default.sound.file))
  {
    default <- .__default.sound.file
  } else
  if (file.exists("c:/windows/media"))
  {
    default <- "c:/windows/media/*.*"
  } else
  {
    default <- ""
  }
  default <-
      choose.files(default = default, caption = "Select a sound file", multi = FALSE,
                   filters = matrix(c("windows wave files","*.wav"),1), index = 1)
  assign(".__default.sound.file", default, "library:ColinsLib")
  return ( TRUE )
}

default.sound <-
function ()
{
  default <- "c:/windows/media/tada.wav"
  if (file.exists(default))
  {
    assign(".__default.sound.file", default, "library:ColinsLib")
  } else
  {
    set.sound()
  }
  return ( TRUE )
}

