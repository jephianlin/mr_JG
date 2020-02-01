**This code is a copy of  [Minimum Rank Library](https://github.com/jasongrout/minimum_rank) from Jason Grout's GitHub.  
**Changes are made mainly to make the code compatible to current version of Sage.  

Minimum Rank Sage Library
=========================

One goal of this library is to be able to be used both by loading the necessary files in the Sage notebook via load statements, as well as installation as a python module in Sage.  Thus, it is important that the names in the separate submodules do not trample on each other (since the load command imports everything at the top level)

See http://sage.cs.drake.edu/home/pub/69/ for an example of how to load and use this library in Sage.  In particular, the following code in a Sage notebook cell will load this library:

```python
URL='https://raw.githubusercontent.com/jephianlin/mr_JG/master/'
files=['Zq_c.pyx','Zq.py','zero_forcing_64.pyx','zero_forcing_wavefront.pyx','minrank.py', 'inertia.py']
for f in files:
    print("Loading %s..."%f);
    load(URL+f)
```
  
CoCalc user
-----------

For CoCalc (previously SageMathCloud) user, if you are using a free acount, then the code above does not work for you.  It is because the connection from CoCalc to outside is forbidden.  You may do the following:

1. Click the green button "Clone or download" and choose "Download Zip".
2. Go to your project in CoCalc and upload this zip file.
3. Create a terminal and open it.
4. Type in the terminal::
```bash
unzip mr_JG-master.zip
```
5. Open your Sage worksheet and execute::
```python
URL='mr_JG-master/'
files=['Zq_c.pyx','Zq.py','zero_forcing_64.pyx','zero_forcing_wavefront.pyx','minrank.py', 'inertia.py']
for f in files:
    print("Loading %s..."%f);
    load(URL+f)
```
If you are tired of downloading and uploading all the time, you may consider to subscribe CoCalc, which give you internet access from CoCalc to outside (in particular, GitHub).


List of major functions
-----------------------

#### `Zq_c.pyx`:
* `push_zeros`
* `push_zeros_looped`
* `neighbors_connected_components`

#### `Zq.py`:
* `subsets`
* `zero_forcing_sets`
* `Z_pythonBitsetold`
* `Zq_inertia_lower_bound`
* `plot_inertia_lower_bound`
* `Zq_graph_info`
* `Zq_bitset`
* `Zqhat_recurse`
* `Zqhat`
* `Zq_compute`
* `Zplus`
* `check_trees`

#### `zero_forcing_64.pyx`:
* `zero_forcing_set_bruteforce_cython_connected`
* `zero_forcing_set_bruteforce_cython`

#### `zero_forcing_wavefront.pyx`:
* `zero_forcing_set_wavefront`

#### `minrank.py`:
* `get_mr_from_list`
* `zerosgame`
* `zero_forcing_set_bruteforce`
* `find_Z`
* `has_forbidden_induced_subgraph`
* `min_rank_by_bounds`
* `find_rank_spread`
* `cut_vertex_connected_graph_mr`
* `minrank_bounds`
* `is_outerplanar`
* `cut_vertex_balanced`
* `cliques_containing_edge`
* `edge_clique_cover_minimum`

#### `inertia.py`:
* `InertiaSet` (class)
* `basic_inertia_set`
* `inertia_set`
* `find_Mq_list`
