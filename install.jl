using Pkg

Pkg.activate("ChimericAnalysis")
Pkg.add(url="https://github.com/maltesie/RNASeqTools#master")
Pkg.instantiate()
Pkg.build()
Pkg.activate("ChimericBrowser")
Pkg.instantiate()
Pkg.build()
#Pkg.add(url="https://github.com/maltesie/RNASeqTools#master")
#Pkg.build("RNASeqTools")
#Pkg.develop(path="./ChimericAnalysis")
#Pkg.build("ChimericAnalysis")
#Pkg.develop(path="./ChimericBrowser")
#Pkg.build("ChimericBrowser")

println("Install finished.")
