module Humbert

# imported dependencies
using Oscar

# files defining functions
include("main.jl")
include("utils.jl")

# public functions
export primeList, polyForm, fileParser, classNumbers, getEll, polz, polzSQI, RHI, degForm, allRHI, minDeg

end # module Humbert
