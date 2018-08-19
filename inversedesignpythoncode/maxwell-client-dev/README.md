Matlab client for the new maxwell.

# Universal server-client interface

`maxwell-client` uses a universal interface to communicate
 with backend Maxwell simulation servers. 
This interface is briefly described below.


## Communication model

Maxwell's client-server communication model is basically just HTTP
 with a few wrinkles here and there.
As such, the basics are simply:

*   `POST` requests are used to upload files.
*   `GET` requests are used to download files.

The wrinkles are:

*   Uploading a `request-*` file (via `POST`) indicates that a simulation 
    is ready to be solved.
*   An empty `GET` request returns information on the simulation queue
    (i.e. how many problems are currently in the queue).


## File model

Now that you know how communication works, here's what the server expects 
 to be communicated. In terms of input files:

*   A `*.grid` file which contains datasets `/omega_r` and `/omega_i`, 
    as well as 12 additional datasets for the `s_prim` and `s_dual` values:
    `/sp_xr`, `/sp_xi`, ..., `/sd_xr`, ..., `/sd_zi`.

*   24 files total for the components of `eps`, `mu`, `J`, and `E0`,
    which use prefixes `e`, `m`, `J`, and `E` respectively.
    In other words, 24 files from `*.e_xr` to `*.E_zi`.
    Each of these are broken down into 3 directional components 
    which each have a real and imaginary component.
    Each component has it's own file.
    Data is stored in the `/data' dataset of the file.
    Lastly, note that compressed (lossless) hdf5 files are almost always
    preferred for these input files since almost all the values 
    are heavily repeated.

*   Lastly, a `request-*` text file with arbitrary contents, 
    although it may be used to hold login data in the future.
    Note that this should be the *last* file uploaded
    since it signals to the server that it should now have all the files it 
    needs to run the simulation.


And output files:

*   12 files total for the components of `E` and `H`:
    `*.E_xr` to `*.H_zi`.
*   A `*.status` file which contains the residual data for every iteration
    of the solve.
*   A `*.log` file which contains the contents of stdout and stderr 
    of the solver.
*   A `*.finished` file which signals that the solve is complete 
    (and thus ready to be downloaded).

The questions now remains, what should the `*` be?
Well, pretty much anything, but some kind of quasi-unique identifier 
 would be best.
On the server-side, all files are sorted by IP so there's no need
 to worry about name collisions with other users, only yourself.










