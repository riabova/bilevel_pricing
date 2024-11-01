# Mathematical model
Consider an online retail company offering a set of products and a delivery location choice to their customers. The delivery options include home delivery or delivery to one of the retailer’s pick-up locations,
which we will further call stores. Since the retailer is interested in some customers choosing particular store locations, they provide individual discounts for the corresponding delivery options. Each customer has a subset of offered products in their cart, and for each product, they choose a delivery option: either home or a store location. After the customers make their delivery choices, the retailer solves a routing problem on the chosen locations (homes and stores) as illustrated below. 

<p align="center">
<img src="https://github.com/riabova/bilevel_pricing/blob/main/img/game_scheme1a.png" width="400">
</p>

## Bilevel formulation
This decision-making process can be described as a single-leader multi-follower bilevel optimization problem defined over the set of *K* customers, *S* stores, *P* products and *R* vehicles.

| Symbol                                      | Description                                                                        |
|---------------------------------------------|------------------------------------------------------------------------------------|
| **Sets**                                    |                                                                                    |
| $\mathcal{K}$                               | Set of customers                                                                   |
| $\mathcal{S}$                               | Set of stores                                                                      |
| $\mathcal{P}$                               | Set of products                                                                    |
| $\mathcal{R}$                               | Set of vehicles                                                                    |
| $\mathcal{N} = \{0\}\cup K \cup S$          | Set of all nodes (depot + homes + stores)                                          |
| $\mathcal{P}_k \subseteq \mathcal{P}$       | Set of products considered by customer $k\in \mathcal{K}$                          |
| $\mathcal{S}_k \subseteq \mathcal{S}$       | Set of stores available to customer $k\in \mathcal{K}$                             |
| $\mathcal{S}'_k = \mathcal{S}_k \cup \{k\}$ | Set of delivery options for customer $k\in \mathcal{K}$                            |
| **Parameters**                              |                                                                                    |
| $c_{ij}$                                    | Cost of traversing arc $(i, j) \in \mathcal{N}\times\mathcal{N}$                   |
| $q_r$                                       | Capacity of vehicle $r \in \mathcal{R}$                                            |
| $w_p$                                       | Weight of product $p \in \mathcal{P}$                                              |
| $h_p$                                       | Intrinsic inconvenience of product $p\ in \mathcal{P}$                             |
| $u_p^k$                                     | Utility of product $p\in \mathcal{P}_k$ for customer $k\in\mathcal{K}$             |
| $I_s^k$                                     | Inconvenience of visiting store $s\in\mathcal{S}_k$ for customer $k\in\mathcal{K}$ |
| $z_p^{base}$                                | Base (market) price of product $p\in \mathcal{P}$                                  |
| **Leader's variables**                      |                                                                                    |
| Discount                                    |                                                                                    |
| $z_{ip} \geq 0$                             | discount for delivery option $(i, p)$                                              |
| Routing                                     |                                                                                    |
| $x_{ij}^r \in \{0, 1\}$                     | 1 if vehicle $r$ uses edge $(i, j)$, 0 otherwise                                   |
| $x_{i0}^r \in \{0, 1, 2\}$                  | for edges connecting depot                                                         |
| Auxiliary                                   |                                                                                    |
| $g_i^r \in \{0, 1\}$                        | 1 if vehicle $r$ needs to visit location $i$, 0 otherwise                          |
| $v_{ip}^r \in \{0, 1\}$                     | 1 if vehicle $r$ carries product $p$ to location $i$, 0 otherwise                  |
| **Follower $k$'s variables**                |                                                                                    |
| $y_{ip} \in \{0, 1\}$                       | if customer chooses buying option $(i, p)$, 0 otherwise                            |
| $\pi_s^k \in \{0, 1\}$                      | if follower $k$ needs to visit location    $s$, 0 otherwise                        |

The upper-level problem is given by

The bilevel model for the described problem is then given by:

**Maximize:**

```math
\max_{x, y, z} \quad \sum_{i \in \mathcal{S}_k'} \sum_{p \in \mathcal{P}_k} (z_p^{base} - z_{ip}) y_{ip} - \sum_{r \in \mathcal{R}} \sum_{i \in \mathcal{N}} \sum_{j \in \mathcal{N} \setminus \{i\}} c_{ij} x_{ij}^r \quad (1a)
```

**Subject to:**

```math
\sum_{j \in \mathcal{N} \setminus \{i\}} x_{ij}^r = 2g_i^r \quad \forall i \in \mathcal{N}, ~ r \in \mathcal{R} \quad (1b)
```
```math
\sum_{(i, j) \in E(S)} x_{ij}^r \leq |S| - 1 \quad \forall S \subset \mathcal{N} \quad (1c)
```
```math
\sum_{p \in \mathcal{P}} w_p \sum_{i \in \mathcal{N} \setminus \{0\}} v_{ip}^r \leq q_r \quad \forall r \in \mathcal{R} \quad (1d)
```
```math
g_i^r \geq v_{ip}^r \quad \forall i \in \mathcal{N} \setminus \{0\}, ~ p \in \mathcal{P}, ~ r \in \mathcal{R} \quad (1e)
```
```math
g_0^r \geq v_{ip}^r \quad \forall i \in \mathcal{N} \setminus \{0\}, ~ p \in \mathcal{P}, ~ r \in \mathcal{R} \quad (1f)
```
```math
\sum_{r \in \mathcal{R}} v_{ip}^r = y_{ip} \quad \forall i \in \mathcal{N} \setminus \{0\}, ~ p \in \mathcal{P} \quad (1g)
```
```math
z_{ip} \leq z^{base}_p \quad \forall i \in \mathcal{S}, ~ p \in \mathcal{P} \quad (1h)
```
```math
y_{ip} \in argmax_y \{f_k(y) \text{ s.t. } y \in \mathcal{Y}_k\} \quad \forall k \in \mathcal{K}, ~ i \in \mathcal{S}'_k, ~ p \in \mathcal{P}_k \quad (1i)
```

**Variables:**
```math
x_{ij}^r, ~ g_i, ~ v_{ip}^r \in \{0, 1\}, ~ z_{ip} > 0
```

The leader’s objective (1a) maximizes the revenue calculated as net sales minus delivery (routing) costs. Constraints (1b) are the degree constraints for each node. Exponential constraint set
(1c) contains subtour elimination constraints for which we use standard delayed constraint generation in a callback loop. Constraints (1d) restrict the capacities of the vehicles. Constraints (1e)
and (1f) link the need to visit the location with the fact that a vehicle carries an item to that location. Constraints (1g) make sure that each item is carried by exactly one vehicle. Constraints (1h)
guarantee the optimality of the K lower-level problems. Yk represents the feasible set of follower’s
k problem.

Each follower’s discrete choice problem is defined as:

**Maximize:**

```math
\max_{y, \pi} \quad f_k(y) = \sum_{i \in S'_k} \sum_{p \in \mathcal{P}_k} \left(u_p - (z_p^{base} - z_{ip})\right) y_{ip} - \sum_{s \in \mathcal{S}_k} I_s^k \pi_s^k + \sum_{p \in \mathcal{P}_k} h_p \sum_{s \in \mathcal{S}_k} y_{sp} \quad (2a)
```

**Subject to:**

```math
   \sum_{i \in \mathcal{S}'_k} y_{ip} = 1 \quad \forall p \in \mathcal{P}_k \quad (2b)
```
```math
   y_{sp} \leq \pi_s^k \quad \forall s \in \mathcal{S}_k, ~ p \in \mathcal{P}_k \quad (2c)
```
```math
   y_{ip} + y_{jp'} \leq 1 \quad \forall (i, j) \in \mathcal{D}, \quad p, p' \in \mathcal{P}_k \quad (2d)
```

**Variables:**
```math
   y_{ip} \in \{0, 1\}
```
Objective function (2a) maximizes the utility customer k gets from their purchases. It is given
by the difference between the monetary equivalent of product utility for that customer and the price
they pay for this product minus the monetary measure of inconvenience they perceive from visiting
locations other than their home plus the sum of parameters intrinsic to the products that represent
how convenient it is to pick them up. These parameters are unrestricted in sign meaning that a
positive value of hp stands for a preference toward pick up (for example a very expensive item for
which a more secure delivery is preferred) while a negative value of hp stands for the inconvenience
associated with carrying the product (large or heavy items). Constraints (2b) ensure that exactly
one buying option is chosen for each product. Constraint set (2c) determines if the customer pays
the inconvenience of going to each particular store location. The packing constraints in (2d) restrict
the pairs of conflict locations.

## Relaxed single-level reformulation

Taking a closer look at the lower-level problem (2), it is
not hard to see that without packing constraints (2d), the problem is total dual integral, which
means that the optimal solution to that problem can be obtained by solving its linear relaxation.

We can obtain the following single-level reformulation for our problem using the relaxed lower-level problem:

**Maximize**
```math
\max_{x, y, z}\quad\sum_{i \in \mathcal{S}_k'}\sum_{p\in \mathcal{P}_k} (z_p^{base} - z_{ip})y_{ip} - \sum_{r \in \mathcal{R}}\sum_{i \in \mathcal{N}}\sum_{j \in \mathcal{N}\setminus\{i\}}c_{ij}x_{ij}^r \quad (3a)
```

**Subject to**
```math
\sum_{j\in \mathcal{N}\setminus\{i\}} x_{ij}^r = 2g_i^r \quad \forall i \in \mathcal{N}, ~r \in \mathcal{R} \quad (3b)
```
```math
\sum_{(i, j)\in E(S)}x_{ij}^r\leq |S| - 1 \quad \forall S\subset \mathcal{N} \quad (3c)
```
```math
\sum_{p\in \mathcal{P}}w_p\sum_{i\in N\setminus\{0\}}v_{ip}^r \leq q_r \quad\forall r \in \mathcal{R} \quad (3d)
```
```math
g_i^r \geq v_{ip}^r \quad \forall i \in \mathcal{N}\setminus\{0\}, ~p\in \mathcal{P}, ~r \in \mathcal{R} \quad (3e)
```
```math 
g_0^r \geq v_{ip}^r \quad \forall i\in \mathcal{N}\setminus\{0\}, ~p\in \mathcal{P}, ~r \in \mathcal{R} \quad (3f)
```
```math
\sum_{r\in \mathcal{R}}v_{ip}^r = y_{ip} \quad \forall i\in \mathcal{N}\setminus\{0\}, ~p\in \mathcal{P} \quad (3g)
```
```math
z_{ip} \leq z^{base}_p \quad \forall i \in \mathcal{S}, ~ p \in \mathcal{P} \quad (3h)
```
```math
\sum_{p\in \mathcal{P}_k}\alpha_p = \sum_{i \in \mathcal{S}'_k}\sum_{p\in \mathcal{P}_k}(u_p -(z_p^{base} - z_{ip})))y_{ip} - \sum_{s \in \mathcal{S}_k} I_s^k\pi_s^k + \sum_{p\in \mathcal{P}_k} h_p\sum_{s\in \mathcal{S}_k}y_{sp} \quad \forall k\in \mathcal{K} \quad (3i)
```
```math
\sum_{i \in \mathcal{S}'_k}y_{ip} = 1 \quad \forall p\in \mathcal{P}_k \quad (3j)
```
```math
y_{sp} \leq\pi_s^k \quad \forall s \in \mathcal{S}_k, ~p\in \mathcal{P}_k \quad (3k)
```
```math
\alpha_p + \delta_{kp} \geq u_p - z_{p}^{base} \quad\forall p\in \mathcal{P}_k \quad (3l)
```
```math
\alpha_p + \beta_{sp} + \delta_{sp} \geq u_p - (z_p^{base} - z_{sp}) + h_p \quad\forall s\in \mathcal{S}_k, p\in \mathcal{P}_k \quad (3m)
```
```math
\sum_{p\in \mathcal{P}_k}\beta_{sp} - \theta_s^k \leq I_s^k \quad\forall s\in \mathcal{S}_k \quad (3n)
```
**Variables**
```math
x_{ij}^r, ~g_i^r, ~v_{ip}^r, ~y_{ip} \in \{0, 1\}, ~\beta_{sp}, ~\delta_{ip}, ~\theta_s^k, ~z_{ip} > 0
```

Here (3a) - (3h) is the upper-level problem described above, constraints (3i) represent strong duality for the lower-level problem, constraints (3j) and (3k) are primal feasibility constraints from (2) 
and constraints (3l) - (3n) present dual feasibility.

## Restoring integrality

The relaxed packing constraints essentially say that for each conflicting pair of
locations, either one or both corresponding variables should be zero. Therefore, when we encounter
a violation of such a constraint, we can create two branches, where we fix one of the two variables
to be equal to zero and solve the two resulting problems with fixed variables. The Branch-and-Bound procedure pictured below starts by solving the relaxed single-level problem at the root node. Then it checks for the packing constraints violation and if a conflict pair of locations is
found, creates two branches, each one of them forces one of the conflicting delivery options to be
equal to zero.

<p align="center">
<img src="https://github.com/riabova/bilevel_pricing/blob/main/img/bnb.png" width="600">
</p>
