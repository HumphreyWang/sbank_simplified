
# What is `sbank` ?

`sbank` is a PYTHON package for stochastic gravitational waveform template bank generation.
If you are interested in the original code, please check [this link](https://github.com/gwastro/sbank) for details.

![sbank_structure](./docs/sbank_structure.png?raw=true)
This is an illustration of the `sbank` structure drawn by me.

# What is `sbank_simplified` ?

This was originally a coursework for *Computational Astrophysics(2023)* at SYSU, and then I decided to make it public.

---

True to its name, this is a(n over-)simplified version of `sbank`. I removed all stuff about gravitational waves
and only kept the key algorithms like how we regard one template as near the other one,
and how we filter those proposals that have already been covered in generating a bank.

---

To make it even simpler (so that can be illustrated easily),
here we do not calculate the inner product of two signals (no need for a detector's PSD),
but just calculate the proper distance between two points in parameter space,
and we only consider Euclidean-like 2D space. You can see the simplified structure below:

![sbank_structure](./docs/sbank_simplified_structure.png?raw=true)

# Get started

There are some examples in `sbank_tests.sh`, and their results are shown in the ` \examples\ ` folder.

We use Parser to parse input parameters, you can type
```shell
python3 sbank.py --help
```
in command line for help docs. You can also check files in `\docs\` for further introductions.

Enjoy and hope it lets you have a better understanding of what a bank generation actually does
and what a template bank should *look* like :)

(This PYTHON code does not directly generate videos, it generates a figure set that is later combined into a video.)

# One more thing

There an additional file `matched_filtering.py`,
which is used for illustrating the matched filtering method (GW150914 as an example).
This file requires [PyCBC](https://pycbc.org/), you can install it by 
```shell
pip install pycbc
```
