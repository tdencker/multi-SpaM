# About

`Multi-SpaM` is a program to infer phylogenies for a set of genomes. It is based on 4-sets of spaced words that are highly likely to represent true homologies. These blocks are evaluated with the Maximum-Likelihood method `RAxML` and the resulting trees of 4 sequences, or quartet trees, are amalgamated into a supertree (using the `Quartet MaxCut` tool). Since the number of blocks used is limited, it is suitable even for large datasets.

All information around installation and standard usage can be found in this readme.

# Installation and Usage

Currently, `multi-SpaM` is only available on 64-bit linux distributions (due to the limitation of the `Quartet MaxCut` tool).

In order to install `multi-SpaM` simply use in the base directory:

	$ make

The program itself can then be used via the python script `multispam.py`. Example usage:

	$ python multispam.py -t <number of threads> -i <input file> -o <output file>

where the input file is a FASTA file containing multiple genomes. The output file will be a tree file in newick format.

# Options

-i              | Input file in FASTA format

-o              | Output file in newick format

-k / -w         | Weight of the pattern (i.e. the number of matching positions) [ can't be larger than 16 ]

-d              | Number of don't care positions (i.e. the number of positions that don't have to match)

-t              | Number of threads used

-n              | Number of sampled blocks

--mem-save / -m | Memory save mode (higher runtime, but much less RAM usage for larger files)


**Tips:**
In general, the parameters don't have to be changed. Only the number of threads, input and output need to be specified. If the resulting trees seem unreasonable, you can try lowering the number of don't care positions to 50. In case of large input files, it is recommended to increase the weight to 12 or even higher. Also, if you have rather limited RAM, you can use the memory save mode. For input files larger than 200 mb or so, the required RAM will exceed 8 gb. With the memory saving mode, the RAM requirement could be reduced to 10.5 gb for a 4.8 gb dataset (doubling the runtime). The number of sampled blocks doesn't have to be increased unless (potentially) for very large datasets.

## Additional Resources

Additional information can be found in our [paper](https://arxiv.org/abs/1803.09222).

## License

Copyright © 2018 - Thomas Dencker
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

Some files may be licensed differently.

## Contact

In case of bugs or unexpected errors don't hesitate to send me a mail: thomas.dencker@stud.uni-goettingen.de
