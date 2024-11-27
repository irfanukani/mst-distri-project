
# Minimum Spanning Tree | MPI

This projects aims to finding Minimum Spanning Tree from a given weighted graph. 

## Run Locally

Clone the project

```bash
  git clone https://github.com/irfanukani/mst-distri-project
```

Go to the project directory and compile the mst.cpp file.

```bash
  mpic++ mst.cpp -std=c++20
```

For generating input you can provide your own txt file or you can create one using:

```bash     
  python generator.py > input.txt
```

Now, we can run the algorithm and it'll print the output:

```bash
  mpiexec -n <num_process> a.out <input_file_path>
```
## Authors

- [@irfanukani](https://www.github.com/irfanukani)
- [@ram-ek](https://www.github.com/ram-3k)



