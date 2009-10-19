
#
#  error function
#
error <-
function(code, call. = match.call()) # think about a better default value ...
{
  error.message <-
  if (!missing(code) && valid.error(code))
  {
    error.message(code)
  } else
  {
    "unknown or unprovided error."
  }
  call. <- paste(as.expression(call.))
  error.message <- paste(error.message, " in function: ", call., sep="")
  stop(error.message, call.=FALSE)
}

#
# check if error code is valid
#
valid.error <-
function (code)
{
  if (code >= 1L && code <= 1L)
    return ( TRUE )
  else
    return ( FALSE )
}

#
# get error message
#
error.message <-
function (code)
{
  switch(code,
    "Graph bounds exceeded"
  )
}

#
# a list of error codes
#
GRAPH_ERROR <- 1L

