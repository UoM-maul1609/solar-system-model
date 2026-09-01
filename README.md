# Solar System Model (SSM)

The Solar System Model (SSM) is a Fortran model for integrating the 3-D
trajectories of ten Solar System bodies. It solves the Cartesian equations of
motion using Newtonian gravity and DVODE, with an optional approximate
relativistic correction referenced to the Sun.

The source code contains a longer Doxygen description of the governing
equations.

## Bodies

The model uses the following fixed ordering:

1. Sun
2. Mercury
3. Venus
4. Earth
5. Mars
6. Jupiter
7. Saturn
8. Uranus
9. Neptune
10. Pluto

## Requirements

- Fortran compiler (`gfortran` for the default Linux build)
- NetCDF Fortran
- the bundled `osnf/` numerical library

## Building

The normal build remains:

```bash
make
```

This retains the existing optimised heliocentric build.

Two optional convenience targets are also provided:

```bash
make release
make debug
```

The original `PLATFORM=LINUX` / `PLATFORM=WIN64` mechanism is retained.

### NetCDF configuration

Existing local configurations remain supported. For example:

```make
NETCDF_FOR=/path/to/netcdf-fortran
NETCDF_C=/path/to/netcdf-c
NETCDF_LIB=-lnetcdff
```

You may also continue to provide `NETCDFLIB` and `NETCDFMOD` directly.

If no legacy NetCDF paths are supplied on Linux, the Makefile tries:

```bash
nf-config --fflags
nf-config --flibs
```

The assembled flags can also be overridden directly with
`NETCDF_FFLAGS` and `NETCDF_FLIBS`.

## Running

Run directly with:

```bash
./main.exe namelist.in
```

or use the supplied helper script:

```bash
./run.sh
```

## Namelist

### `&initial_state`

The input units are:

| Variable | Units |
| --- | --- |
| `x`, `y`, `z` | km |
| `ux`, `uy`, `uz` | km s-1 |
| `gm` | km3 s-2 |

The defaults in `solar_system.f90` are now stored in these same units and are
converted to SI once after the namelist is read. This prevents default values
from being converted twice in a partially specified namelist.

### `&run_vars`

| Variable | Meaning |
| --- | --- |
| `outputfile01` | output NetCDF filename |
| `run_forward_in_time` | reverse the initial velocities when false |
| `general_relativity` | enable the existing approximate solar GR correction |
| `dt` | output interval in days; also DVODE maximum step |
| `tfinal` | integration duration in years |
| `interval_io` | interval between progress messages in years |
| `interact(1:10)` | bodies enabled as Newtonian gravity sources |

No new solver controls have been added to the namelist. The existing DVODE
configuration and tolerances remain in the source.

## Interaction mask

`interact(j)=.true.` means body `j` contributes to the Newtonian acceleration.

With only the Sun enabled, for example:

```fortran
interact(1:10)=.true.,.false.,.false.,.false.,.false.,
               .false.,.false.,.false.,.false.,.false.
```

all planets feel the Sun but do not gravitationally interact with one another.

A body with no active Newtonian gravity sources still evolves according to
its velocity (`dx/dt=v`) rather than being frozen.

## Relativistic correction

The original relativistic approximation has been retained. The implementation
now explicitly calculates the correction relative to body 1 (the Sun), rather
than relying on the first enabled entry in the interaction list.

The `general_relativity` switch remains independent of the Newtonian
`interact` mask.

## Time integration

The existing DVODE setup, work arrays, tolerances and `ISTATE` handling are
unchanged.

The only integration-loop change is that the next requested output time is
limited to the final time:

```fortran
tout=min(tt+dt,tfinal)
```

This prevents the final requested integration time from overshooting `tfinal`
when the run duration is not an exact multiple of `dt`.

## NetCDF output

The trajectory output remains single precision to keep file size down:

| Variable | Units |
| --- | --- |
| `time` | s |
| `pos` | m |
| `vel` | m s-1 |

The output also includes:

| Variable | Description |
| --- | --- |
| `gm` | gravitational parameter of each body |
| `gravity_source` | integer form of the `interact` mask |

Global metadata records the initial epoch, body ordering, namelist filename,
integration direction, relativistic switch, heliocentric build switch,
`dt`, `tfinal`, and `interval_io`.

The NetCDF file remains open during the integration. Buffered output is
written and `nf90_sync` is called when the buffer is flushed, avoiding repeated
open/close cycles.

## Notes on this update

This update is intentionally conservative. It does **not**:

- change DVODE behaviour or error handling;
- add solver tolerances to the namelist;
- change trajectory output to double precision;
- remove currently unused state or diagnostic variables;
- restructure the original scientific/Doxygen comments.

The changes are limited to the selected unit, final-time, source-free motion,
solar-GR reference, NetCDF-output/metadata, and build portability fixes.
