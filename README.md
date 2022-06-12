# phys415-chargedparticles

Simulation of charged particle motion in arbitrary electric and magnetic fields.
Currently supports an arbitrary number of particles modelled using the Boris method. Interpolation of fields, which may be loaded from HDF5 data, is also supported.

## Requirements and Setup
To run this code you will need Julia. The [latest stable release](https://julialang.org/downloads/#current_stable_release) (1.7.3 at time of writing) is recommended. You will also need some additional packages, listed in `Project.toml`. The easiest way to install them is to use Julia's [inbuilt package manager](https://pkgdocs.julialang.org/v1/), as follows:
1. Clone the repo and `cd` inside
2. Start a Julia REPL session with `julia`
3. Switch to package management mode by typing `]`.

4. Activate the project with `activate .`
5. Run `instantiate` to install the dependencies defined in the project manifest.
 
To run a Julia script inside this environment use `julia --project=. BorisSim.jl` from your terminal. If you prefer to launch from the REPL, just run `include("BorisSim.jl")` with the project environment active, as described above. You'll need to exit package management first, by pressing backspace on an empty line in the REPL. `BorisSim.jl` is the main file. Running the others will achieve little.

### Updating
If you have previously run this code run through the setup steps again to get the latest dependencies.

## Configuration
To define the parameters of the scenario, see the code comments. Options of note are the number of particles (line 125 at present) and the field definitions (third block). Be aware that combining interpolation and large numbers of particles will make things rather slow.

There is an *experimental* animated output option using Javis.jl. Just uncomment the makeRender call to try it.

The `interpolations.jl` includes some test functions. These are currently commented out but do work. The 2D function has graphical output but the 3D function does not, given the difficulty in clearly plotting a three dimensional scalar field.

The default configuration simulates 100 particles with randomly generated initial properties in uniform z-directed electric and magnetic fields. Other included examples are interpolation over data from an HDF5 file and interpolation over a uniform magnetic field. These examples have been designed to show clear gyration in the output. Working with, for example, a smaller magnetic field may make this effect unnoticeable over the simulation period.