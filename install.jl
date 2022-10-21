using Pkg

Pkg.develop(path="./ChimericAnalysis")
Pkg.build("ChimericAnalysis")
Pkg.develop(path="./ChimericBrowser")
Pkg.build("ChimericBrowser")

println("Install finished.")
