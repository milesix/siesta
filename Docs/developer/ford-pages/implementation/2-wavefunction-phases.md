title: Wavefunction phase convention

When dealing with periodic systems and Brillouin Zone sampling, one inevitably
arrives to the definition of Bloch-like basis functions from a given basis set.
The definition of said functions can follow one of two conventions, depending
whether one prefers to use something similar to a regular Bloch function
\(\chi_{j}\) (Convention I) or a *cell periodic* Bloch function
\(\tilde\chi_{j}\) (Convention II).

$$ \left |       \chi_{j} \right > = \sum_{ \mathbf R}
   e^{i \mathbf k·(\mathbf R + \tau_{j})}
   \left | \phi_{j}(\mathbf R) \right >  $$

$$ \left | \tilde\chi_{j} \right > = \sum_{ \mathbf R}
   e^{i \mathbf k· \mathbf R}
   \left | \phi_{j}(\mathbf R) \right >  $$

The only difference is the addition of the phase term \( e^{ i \mathbf k
· \tau_{j} } \). In SIESTA we follow the first convention.
Of course, under both conventions matrices such as the Hamiltonian matrix are
expressed differently:

$$        H_{ij} = \sum_{ \mathbf R} e^{i \mathbf k·(\mathbf R
                                         + \tau_{j} - \tau_{i})} h_{ij}$$
$$ \tilde H_{ij} = \sum_{ \mathbf R} e^{i \mathbf k· \mathbf R } h_{ij} $$

With \( h_{ij} = \left < \phi_{i}  | \hat H | \phi_{j} \right > \). When solving
the standard eigenvalue problem, this difference gets carried to the eigenvector
coefficients \( C_{ij} \), so that:

$$ \tilde C_{nj} = C_{nj} e^{i \mathbf k · \tau_{j} } $$

Nevertheless, unless we need to explicitly take into account the coefficients,
this becomes transparent when dealing with wavefunctions or density matrices.
For example, given a Bloch function \( \psi_{n} \):

$$ \left | \psi_{n} \right > =
   \sum_{R}\sum_{j} \tilde C_{nj} e^{i \mathbf k·(\mathbf R)}
   \left | \phi_{j}(\mathbf R) \right >  =
   \sum_{R}\sum_{j} C_{nj} e^{i \mathbf k·(\mathbf R + \tau_{j})}
   \left | \phi_{j}(\mathbf R) \right > $$