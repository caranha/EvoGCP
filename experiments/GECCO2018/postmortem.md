# Comparison of SI approaches for 3-GCP -- A postmortem

Submitted for GECCO2018, Accepted as a poster

## Outline
This paper compare 4 swarm intelligence approaches (ABC, CS, FFA, PSO) 
for solving the 3-GCP problem. The four approaches try to minimize 
the number of violations in the 3-coloring of a graph.

The algorithms are tested using the parameters suggested by the literature,
as well as parameters suggested by the application of the iRACE script.
The results highlight the differences in performances of the two 
parameter sets, and allow some discussion about how sensitive each algorithms
is to parameter setting.

This paper was a preparation for a wider study on a component-based 
approach for Bio-inspired algorithms.

## Files

- EvoGCP_0.3.1.tar.gz: Version of the code used for the experiments

### Data
- graphs.90.to.180.Rda contains the graph sets used in the experiments
- reproduceV3.Rda contains data frames with the results of the experiments
- irace_tuning_X contains the result data for the irace parameter optimization script

### images
- Contains images summarizing the results of the experiments

### scripts
- irace_analysis.R used to generate the images/data from the irace data results
- process_GCP_all_litpars.R used to generate the summarizing images from the experiments
- reproduceV3.R used to run the experiments

## Reviewer Comments
- Differentiation between decision and optimization problem

Two reviewers argued that messed up the problem definition. The Graph Coloring Problem 
can be seen as an optimization problem (find a minimum number of colors to colorize a 
certain graphs) or a decision problem (decide whether a graph is colorable or not). In this 
paper, it was not clear whether we were going for A or for B. 

The reviewers suggested that meta-heuristics are better suited for treating the GCP as a 
minimization problem, and using it this way allows/requires us to use both solvable and insolvable 
graph instances.

- Use of DIMACs benchmarks

It was suggested that we use DIMACs benchmarks for graphs problems in addition to the 
generated benchmarks that we used. 

- Use of other methods

One reviewer was really mad that we didn't add ACO

## TODO

- Add DIMACs benchmarks
- Consider whether to change the definition of the problem (decision problem, optimization problem,
  or keep as a minimization of violation (design) problem).
- Move on to decomposition of SI approaches
