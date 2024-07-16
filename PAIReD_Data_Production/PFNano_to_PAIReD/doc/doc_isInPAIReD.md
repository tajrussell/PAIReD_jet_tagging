# Documentation for the `isInPAIReD` files

Below you can find the documentation for the file `src/tools/PAIReD_geometries/*.py`.

Those modules contain a function which is used for defining the PAIReD jet geometry/shape. A detailed example of the functioning, step by step, can be found in the jupyter notebook [here](../src/notebooks/makeNtuplesPAIReDjointMC.ipynb).

The module requires the following python packages:
* `numpy`
* `awkward`


## Function documentation

### `isInPAIReD` <a name="isInPAIReD">  </a>

**`isInPAIReD(eta_j1, eta_j2, phi_j1, phi_j2, eta_p, phi_p)`**

Determines which of the input particles are inside the PAIReD jet defined by the given pair of AK4 jets. The output array has the shape $(N_e, N_{jp}, N_p)$, where $N_e$ is the number of events, $N_{jp}$ the number of jet pairs (length of the jet input arrays) and $N_p$ the number of particles (length of the particle input arrays). Note that $N_{jp}$ and $N_p$ are event dependend which is why the lengths of these axes vary for the different entries along the event axis.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>eta_j1 : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the first jet, or jets if array is given.</p>
    </dd>
    <dt>eta_j2 : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the second jet, or jets if array is given.</p>
    </dd>
    <dt>phi_j1 : float or array of floats</dt>
    <dd>
      <p>Phi angle of the first jet, or jets if array is given.</p>
    </dd>
    <dt>phi_j2 : float or array of floats</dt>
    <dd>
      <p>Phi angle of the second jet, or jets if array is given.</p>
    </dd>
    <dt>eta_p : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the particle, or particles if array is given.</p>
    </dd>
    <dt>phi_p : float or array of floats</dt>
    <dd>
      <p>Phi angle of the particle, or particles if array is given.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>isInPAIReD : bool or array of bools</dt>
    <dd>
      <p>Boolean holding the information whether the particle is inside the PAIReD jet or not. The with index pair $(i,j)$ holds the information whether particle $j$ is in jet pair $i$.</p>
    </dd>
    </dl>
  </dd>
</dl>

#### Examples

```python
eta_j1 = ak.Array([[0, 1, 2],
                   [],
                   [-2]])
phi_j1 = ak.Array([[1, 2,-3],
                   [],
                   [-1]])

eta_j2 = ak.Array([[1, 2,-6],
                   [],
                   [1]])
phi_j2 = ak.Array([[2,-2, 0],
                   [],
                   [0]])

eta_p = ak.Array([[2, 0,-5, 8],
                  [2,1,3],
                  [7,2,-2]])
phi_p = ak.Array([[2, 2,-1, 3],
                  [-2,-1,3],
                  [1,1,-3]])

isInPAIReD(eta_j1, eta_j2, phi_j1, phi_j2, eta_p, phi_p)
```
```
Out[]:   [[[True, True, False, False], [True, ...], [False, False, True, False]],
         [],
         [[False, False, False]]]
        ------------------------------------------------------------------------
        type: 3 * var * var * bool
```
