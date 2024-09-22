# 2D Medical Torus - System of two hormones

<img src="https://github.com/auan0001/AFEM-Medical-Torus/blob/main/medtorus.gif" width="400" height="400" />

## Running the script
The FEM computation depends on (the now legacy) [FEniCS](https://fenicsproject.org/download/archive/). The files generated can be visualized in [ParaView](https://www.paraview.org/). 
## Numerical Solution
Consider the system of equations modelling a medical torus containing
two hormones

$$
\begin{aligned}
   \partial_t u - \alpha_1 \Delta u&= -uv^{2} + c_f(1-u), &\quad (\boldsymbol x,t) \in \mathcal{B}\times (0,T],\\ 
   \partial_t v - \alpha_1 \Delta v&= uv^{2} - (c_f+c_k)v, &\quad (\boldsymbol x,t) \in \mathcal{B}\times (0,T],\\
   u(\boldsymbol x,0) &= u_0(\boldsymbol x), &\quad \boldsymbol x\in \mathcal{B}, \\
   v(\boldsymbol x,0) &= v_0(\boldsymbol x), &\quad\boldsymbol x\in \mathcal{B}.
\end{aligned}
$$

The objective is to obtain the weak form to find

$$
u, v \in W = \lbrace\boldsymbol w(\boldsymbol x,t), \quad \boldsymbol w(\boldsymbol x,t) \in \mathcal{H}^{1}(\mathcal{B}), \quad \forall t \in [0,T)\rbrace
$$

such that

$$
\begin{aligned}
  \begin{split}
    \int_{\mathcal{B}}\partial_t uw_1 d\boldsymbol x+ \alpha_1\int_{\mathcal{B}}\nabla u \cdot \nabla w_1 d\boldsymbol x+\int_{\mathcal{B}} uv^{2}w_1d\boldsymbol x+ c_f\int_{\mathcal{B}} uw_1d\boldsymbol x&=c_f\int_{\mathcal{B}} w_1d\boldsymbol x, \\
    \int_{\mathcal{B}}\partial_t vw_2 d\boldsymbol x+ \alpha_2\int_{\mathcal{B}}\nabla v \cdot \nabla w_2 d\boldsymbol x-\int_{\mathcal{B}} uv^{2}w_2d\boldsymbol x+(c_f+c_k)\int_{\mathcal{B}} vw_2d\boldsymbol x&=0, \\
    \qquad \forall \boldsymbol w \in W.
  \end{split}
\end{aligned}
$$

by integrating the SE and multiplying with
$\boldsymbol w = (w_1,w_2)^{\top}$. When integrating the Laplacian term,
Green's formula is used, which in this case is integration by parts in
two dimensions with applied homogeneous Neumann boundary conditions in
the outward normal direction of $\partial \mathcal{B}$. From the weak
form, it is possible to derive the Galerkin FEM formulation

$$
W_{h,0} = \lbrace\boldsymbol w(\boldsymbol x,t): \quad \boldsymbol w \in \mathcal{C}^{0}(\mathcal{B})\forall t \in [0,T), \quad w\mid_k \in \mathcal{P}\_1(k), \quad \forall k \in \tau_h\rbrace
$$

such that

$$
\begin{aligned}
  \begin{split}
    \int_{\mathcal{B}}\partial_t u_hw_1 d\boldsymbol x+ \alpha_1\int_{\mathcal{B}}\nabla u_h \cdot \nabla w_1 d\boldsymbol x+\int_{\mathcal{B}} u_hv_h^{2}w_1d\boldsymbol x+ c_f\int_{\mathcal{B}} u_hw_1d\boldsymbol x&=c_f\int_{\mathcal{B}} w_1d\boldsymbol x, \\
    \int_{\mathcal{B}}\partial_t v_hw_2 d\boldsymbol x+ \alpha_2\int_{\mathcal{B}}\nabla v_h \cdot \nabla w_2 d\boldsymbol x-\int_{\mathcal{B}} u_hv_h^{2}w_2d\boldsymbol x+(c_f+c_k)\int_{\mathcal{B}} v_hw_2d\boldsymbol x&=0, \\
    \qquad \forall \boldsymbol w \in W.
  \end{split}
\end{aligned}
$$

The non-linear terms $u_hv_h^2$ from the RHS can be
substituted and by using the Lagrange interpolant $\pi_h \in  W_{h,0}$

$$
\begin{split}
    S &= uv^{2}, \quad S \approx \pi_h S=\sum_{N_j \in  \mathcal{N}\_h}S_{j}\varphi_j\\
  \end{split}
$$

which can be prepared for the Galerkin FEM

$$
\begin{split}
    \pi_h S&=\sum_{N_j \in \mathcal{N}\_h} \xi_{1,j}\xi_{2,j}^{2}\\
  \end{split}
$$

and inserted along with the linear terms consisting of $u_h, v_h$ 

$$
\begin{aligned}
  \begin{split}
    \sum_{N_j\in \mathcal{N}\_h}\partial_t\xi_{1,j} \int_{\mathcal{B}}\varphi_i\varphi_i d\boldsymbol x+\sum_{N_j\in \mathcal{N}\_h}\alpha_1 \xi_{1,j}\int_{\mathcal{B}}\nabla \varphi_j \cdot \nabla \varphi_i d\boldsymbol x \\
    +\sum_{N_j\in \mathcal{N}\_h}\xi_{1,j}\xi_{2,j}^{2}\int_{\mathcal{B}}\varphi_j \varphi_id\boldsymbol x + \sum_{N_j\in \mathcal{N}\_h} c_f\xi_{1,j}\int_{\mathcal{B}}\varphi_j \varphi_id\boldsymbol x
    = c_f\int_{\mathcal{B}}\varphi_j d\boldsymbol x, \\
    \sum_{N_j\in \mathcal{N}\_h}\partial_t\xi_{2,j} \int_{\mathcal{B}}\varphi_i\varphi_i d\boldsymbol x+ \sum_{N_j\in \mathcal{N}\_h}\alpha_2 \xi_{2,j}\int_{\mathcal{B}}\nabla \varphi_j \cdot \nabla \varphi_i d\boldsymbol x\ \\
    -\sum_{N_j\in \mathcal{N}\_h}\xi_{1,j}\xi_{2,j}^{2}\int_{\mathcal{B}}\varphi_j \varphi_id\boldsymbol x + \sum_{N_j\in \mathcal{N}\_h}(c_f+c_k)\xi_{2,j}\int_{\mathcal{B}}\varphi_j \varphi_id\boldsymbol x=0, \\
    \qquad \forall \boldsymbol w \in W.
  \end{split}
\end{aligned}
$$

Identify the matrices

$$
\begin{aligned}
  \begin{split}
    A = \int_{\mathcal{B}}\nabla \varphi_j \cdot \nabla \varphi_i dx, \quad M = \int_{\mathcal{B}}\varphi_j \varphi_i d\boldsymbol x\\  
  \end{split}
\end{aligned}
$$

and see that one of the source terms consists of a load vector

$$
\begin{aligned}
  \begin{split}
    b=c_f\int_{\mathcal{B}}\varphi_j d\boldsymbol x.
  \end{split}
\end{aligned}
$$

In order to discretize the system, it is convenient to
change into matrix ODE notation

$$
\begin{aligned}
  \begin{split}
    M \frac{d\xi_{j,1}}{dt}+\alpha_1 A\xi_{j,1}+MS+c_fM\xi_{j,1} = b\\
    M \frac{d\xi_{j,2}}{dt}+\alpha_2 A\xi_{j,2}-MS+(c_f+c_k)M\xi_{j,2}= 0.
  \end{split}
\end{aligned}
$$

The Crank-Nicolson formulation

$$
\begin{split}
    &M\frac{\xi^{n+1}\_{j,1}-\xi^{n}\_{j,1}}{k_n} + \alpha_1 A\frac{\xi^{n+1}\_{j,1}+\xi^{n}\_{j,1}}{2} + M\xi^{n}\_{j,1}(\xi^{n})^{2}+c_fM\frac{\xi_{j,1}^{n+1}+\xi^{n}\_{j,1}}{2} = b \\
    &M\frac{\xi^{n+1}\_{j,2}-\xi^{n}\_{j,2}}{k_n} + \alpha_2 A\frac{\xi_{j,2}^{n+1}+\xi^{n}\_{j,2}}{2} - M\xi^{n}\_{j,1}(\xi_{j,2}^{n})^{2}+(c_f+c_k)M\frac{\xi^{n+1}\_{j,2}+\xi^{n}\_{j,2}}{2} = 0
  \end{split}
$$

is solved for the next timestep and ready to be implemented as

$$
\begin{split}
    \xi_{j,1}^{n+1}&=\Big(M+\frac{\alpha_1 k_n}{2}A+\frac{c_fk_n}{2}M\Big)^{-1} \
    \Big(\Big(M-\frac{\alpha_1 k_n}{2}A-\frac{c_f k_n}{2}M-k_n M (\xi_{j,2}^{n})^{2}\Big)\xi_{j,1}^{n}+k_nb\Big)\\
    \xi_{j,2}^{n+1}&=\Big(M+\frac{\alpha_2 k_n}{2}A+\frac{(c_f+c_k) k_n}{2}M\Big)^{-1} \
    \Big(M-\frac{\alpha_2 k_n}{2}A-\frac{(c_f+c_k) k_n}{2}M+k_n M \xi_{j,1}^{n}\xi_{j,2}^{n}\Big)\xi_{j,2}^{n}. \\
  \end{split}
  $$
