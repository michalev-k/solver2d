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
Trying to improve on the above I came up with this solver.
Basically it is the TGS Soft equivalent of Mass Splitting.

You can see it here: [video](https://youtu.be/FG2N6LqkdD0)

Solving multiple times per substep seems to increase the quality of the solution, so I made it configurable.
- Main Iters -> Substepping iterations
- Extra Iters -> Inner solution iterations
- Relax Iters -> inner relaxation iterations

Unlike the Mass Splitting Solver it does not use grouping.
Still it provides surprisingly excellent quality.
It is jitter free with few iterations, and can deal with many difficult situations (such as high mass-ratios) very well.

With mass-splitting you would need to do contact reordering inside groups to solve them in parallel.
In comparison this solver is order-independent, and as such lends itself to parallel implementation. (Multithreading, SIMD, GPU-Physics)

That being said, you can also decide to add grouping to this for further improved quality.
I just personally did not like the additional complexity of grouping and contact-reordering.
Maybe it is worth it though. 

Like the original Solver2D code, any changes here are licensed under MIT as well.
