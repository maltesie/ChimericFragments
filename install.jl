using Pkg

Pkg.activate("ChimericAnalysis")
Pkg.instantiate()
Pkg.build()
Pkg.activate("ChimericBrowser")
Pkg.instantiate()
Pkg.build()

println("Install finished.")
