push!(LOAD_PATH,"../src/")
using Documenter, replication_Stock_Watson

makedocs(modules = [replication_Stock_Watson], sitename = "replication_Stock_Watson.jl")

deploydocs(repo = "github.com/GeraudDM/replication_Stock_Watson.jl.git", devbranch = "main")
