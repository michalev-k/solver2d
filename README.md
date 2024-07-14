![Solver2D Logo](extras/solver2d.png)

This is a fork of:
# Solver2D 
A project for testing various rigid body solvers.
Read more about this project [here](https://box2d.org/posts/2024/02/solver2d/).
Here is a video showing test results: [video](https://www.youtube.com/watch?v=sKHf_o_UCzI)

## License
Solver2D is developed by Erin Catto, and uses the [MIT license](https://en.wikipedia.org/wiki/MIT_License).

This fork adds the following solvers:

### Grouped Mass Splitting
From the paper "Mass splitting for jitter-free parallel rigid body simulation"
(Not implemented for joints)

### Mass Split Temporal Jacobi with Soft Contacts
Thinking about how to improve on the above I came up with this solver.
Basically it is the TGS Soft equivalent of Mass Splitting.

Solving multiple times per substep seems to increase the quality of the solution, so I made it configurable.
Main Iters -> Substepping iterations
Extra Iters -> Inner solution iterations
Relax Iters -> inner relaxation iterations

It does not use grouping like in the mass-splitting paper, but still provides surprisingly excellent quality.
With mass-splitting you would need to do contact reordering inside groups to solve them in parallel.
As this solver does not use grouping, it is order-independent
As such it lends itself to a simple parallel implementation, such as Multithreading, SIMD and GPU-Physics.

That being said, you can also decide to add grouping to this for improved quality.
I just personally did not like the additional complexity of grouping and contact-reordering. 
Maybe its worth it though. 

Any changes are licensed under MIT as well.