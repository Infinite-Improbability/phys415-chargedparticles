# phys415-chargedparticles

Simulation of charged particle motion in arbitrary electric and magnetic fields.
Currently supports a single particle modelled using the Boris method.

## Requirements and Setup
To run this code you will need Julia. The [latest stable release](https://julialang.org/downloads/#current_stable_release) (1.7.2 at time of writing) is recommended. You will also need some additional packages, listed in `Project.toml`. The easiest way to install them is to use Julia's [inbuilt package manager](https://pkgdocs.julialang.org/v1/), as follows:
1. Clone the repo and `cd` inside
2. Start a Julia REPL session with `julia`
3. Switch to package management mode by typing `]`
4. Activate the project with `activate .`
5. Run `instantiate` to install the dependencies defined in the project manifest.

To run a Julia script inside this environment use `julia --project=. scriptName.jl`.