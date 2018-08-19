%% maxwell_help
% Welcome to Maxwell's documentation!
% If anything is unclear or needs to be updated please don't hesitate to
% open an issue at <https://github.com/JesseLu/maxwell-client/issues>, 
% or to contact me at <jesselu@stanford.edu>.
% Enjoy!

%%% Introduction to Maxwell
% Along with the functionality of Maxwell, this library includes various
% examples which are useful in familiarizing oneself with Maxwell.

%%% Maxwell functions
%
% * <maxwell_axb.html maxwell_axb> - Matrices and vectors associated with the electromagnetic wave equation.
%
% * <maxwell_box.html maxwell_box> - Box of constant epsilon/mu within the simulation grid.
%
% * <maxwell_box_smooth.html maxwell_box_smooth> - Box of constant epsilon/mu within the simulation grid with smoothed boundaries.
%
% * <maxwell_cyl.html maxwell_cyl> - Cylinder of constant epsilon/mu within the simulation grid.
%
% * <maxwell_cyl_smooth.html maxwell_cyl_smooth> - Cylinder of constant epsilon/mu within the simulation grid with smoothed boundaries.
%
% * <maxwell_flux.html maxwell_flux> - Electromagnetic power passing through a finite plane of the simulation.
%
% * <maxwell_fsmode.html maxwell_fsmode> - Excitation source for arbitrary free-space modes.
%
% * <maxwell_gaussian.html maxwell_gaussian> - Excitation source for free-space Gaussian modes.
%
% * <maxwell_grid.html maxwell_grid> - Initialize a simulation domain.
%
% * <maxwell_pos2ind.html maxwell_pos2ind> - Translate from position to array indices.
%
% * <maxwell_shape.html maxwell_shape> - Add shapes of constant material to the simulation domain.
%
% * <maxwell_solve.html maxwell_solve> - Solve the electromagnetic wave equation.
%
% * <maxwell_solve_async.html maxwell_solve_async> - Electromagnetic solve without waiting for completion.
%
% * <maxwell_solve_eigenmode.html maxwell_solve_eigenmode> - Eigenmode based on an initial guess for E.
%
% * <maxwell_view.html maxwell_view> - View a slice of the field and/or structure.
%
% * <maxwell_wgmode.html maxwell_wgmode> - Excite a waveguide mode. Can also used to filter for a particular mode.


%%% Examples 麦克斯韦方程组各种结构下求解例子
%
% * <example0_waveguide.html example0_waveguide> - Excite the fundamental mode of a waveguide and measure output power.
%
% * <example1_ring.html example1_ring> - Excite the clock-wise顺时针方向的 mode of a ring resonator.
%
% * <example2_metalmode.html example2_metalmode> - Solve for the resonant mode of a metallic金属 resonator.
%
% * <example3_cavitymode.html example3_cavitymode> - Solve for the eigenmodes of dielectric resonators.
%
% * <example4_focusedbeam.html example4_focusedbeam> - Excite Gaussian or donut环状线圈 free-space modes.
%
% * <example5_hiresmetal.html example5_hiresmetal> - Solve for the resonant mode of a metallic resonator.


%%% Introduction to Maxopt
% The |maxopt| family of functions, examples, and cases案例 lives on top of 
% the |maxwell| functions 
% and implements basic structural optimization functions and examples.

%%% Maxopt functions
% These functions enable启用 gradient-based optimization schemes.
%
% * <maxopt_field_gradient.html maxopt_field_gradient> - Calculate structural gradients for the E-field of a simulation.
%
% * <maxopt_freq_gradient.html maxopt_freq_gradient> - Calculate structural gradients for frequency of the eigenmode.
%
% * <maxopt_gradient_descent.html maxopt_gradient_descent> - Simple gradient descent optimization algorithm.
%
% * <maxopt_solve_adjoint.html maxopt_solve_adjoint> - Asynchronously异步 solve for the adjoint伴随 of the operator.


%%% Maxopt examples 优化方法
% Various optimization examples. 
% The various optimization methods showcased展示的 include 
% brute-force parameter sweep, derivative-free search, and gradient-descent.
%
% * <maxopt_example0_sweep.html maxopt_example0_sweep> - Sweep the parameters of a waveguide-coupled disk resonator structure.
%
% * <maxopt_example1_search.html maxopt_example1_search> - Derivative-free optimization of a nanophotonic grating coupler.
%
% * <maxopt_example2_adjoint.html maxopt_example2_adjoint> - Form a cavity out of a square lattice using gradient optimization.
%
% * <maxopt_example3_eigenmode.html maxopt_example3_eigenmode> - Derivative-based optimization of an L3 cavity mode.


%%% Maxopt cases 结构构建，调用不同优化方法
% These cases are used in the Maxopt examples and show how structures can be
% set up for optimization.
%
% * <maxopt_case_2wbeam.html maxopt_case_2wbeam> - Sets up a frequency-doubling cavity optimization.
%
% * <maxopt_case_L3.html maxopt_case_L3> - Sets up the optimization for an L3 photonic crystal resonator.
%
% * <maxopt_case_grating.html maxopt_case_grating> - Sets up a grating coupler optimization problem.
%
% * <maxopt_case_metalfocus.html maxopt_case_metalfocus> - Used to optimize a metal focusing structure.
%
% * <maxopt_case_squarepc.html maxopt_case_squarepc> - Sets up optimization problem for a square photonic crystal array.
%
% * <maxopt_case_wdmgrating.html maxopt_case_wdmgrating> - Used to optimize a wavelength-splitting grating coupler.


% Allows the user to get here just by typing 'maxwell_help' in the command line.
showdemo maxwell_help
