# Detect Module

>   Stochastic Simulation Version

This is the stochastic simulation version of the **Detect Module** designed by our team BNU-China 2022. It's, in respect to biology, a Two-component System (TCS) modified for detecting Δ9-THC in human body, which is the principal psychoactive constituent of cannabis. There are several other modules included in the project as well but this is the one I'm assigned to (at least currently :D). Feel free to contact me at 202011081001@mail.bnu.edu.cn if any question occurs.

Basically, it's based on the **Gillespie/SSA algorithm**, which is usually adopted in the fields of computational biology and chemistry, etc. But as the situation in our project is rather complicated (or exactly, "real"), the open-source versions online are mostly not ideal when simulation. So I came up with this **MODIFIED** version of Gillespie algorithm

## Biological Design

As I'm a student majoring in Computer Science, this is the part that I fail to understand very clearly. And I suppose that most of you come here for the code, not the elaborated, boring explanation of how the system works biologically (which could be wrong actually lol). So I just post the pathway of the Detect Module.

```mermaid
flowchart LR
PpmrB.->RNApmrB.->PmrB
PpmrA.->RNApmrA.->PmrA
RNApmrB.-> A[empty]
PmrB.-> B[empty]
RNApmrA.-> C[empty]
PmrA.-> D[empty]
THC--bind-->PmrB--ATP-->PmrA-->E[PmrA-P]-->PmrC-->F[Target Protein]
PmrC-->G[Target Enzyme]


```

However, I will still leave some references here that support most of my understanding of TCS and Δ9-THC detecting mechanism throughout the developing process. Hopefully it can help you.

| Reference                                                    | URL                                                          |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| Stochastic kinetic model of two component system signalling reveals all-or-none, graded and mixed mode stochastic switching responses | [Stochastic kinetic model of two component system signalling reveals all-or-none, graded and mixed mode stochastic switching responses - PubMed (nih.gov)](https://pubmed.ncbi.nlm.nih.gov/20174681/) |
| Equation-free analysis of two-component system signalling model reveals the emergence of co-existing phenotypes in the absence of multistationarity | https://pubmed.ncbi.nlm.nih.gov/22761552/#:~:text=Equation-free%20analysis%20of%20two-component%20system%20signalling%20model%20reveals,attributed%20to%20the%20inherent%20stochasticity%20of%20biochemical%20processes. |
| The biology of the PmrA/PmrB two-component system: the major regulator of lipopolysaccharide modifications | [The biology of the PmrA/PmrB two-component system: the major regulator of lipopolysaccharide modifications - PubMed (nih.gov)](https://pubmed.ncbi.nlm.nih.gov/23799815/) |

## Gillespie Algorithm

Gillespie/SSA algorithm is a well-known stochastic simulating algorithm for complicated chemical reaction system. 

>   It generates a statistically correct trajectory (possible solution) of a stochastic equation system for which the reaction rates are known. It was created by [Joseph L. Doob](https://en.wikipedia.org/wiki/Joseph_L._Doob) and others (circa 1945), presented by [Dan Gillespie](https://en.wikipedia.org/wiki/Dan_Gillespie) in 1976, and popularized in 1977 in a paper where he uses it to simulate chemical or biochemical systems of reactions efficiently and accurately using limited computational power.
>
>   Wikipedia

To understand how the algorithm works, let's consider a biochemical system, which contains $N$ types of molecules.
$$
\{S_1,S_2,\cdots, S_N\}
$$
We assume that the temperature of the system is stable, that all the molecules are fully mixed and that all reactions can happen immediately. Now, suppose that there are $M$ possible reactions.
$$
\{R_1, R_2,\cdots, R_M\}
$$
And let 
$$
X(t)=(X_1(t), X_2(t), \cdots, X_N(t))
$$
 denote the system's status at $t$, with $X_i(t)$ denoting the molecules of $S_i$ at $t$. Let $v_{ji},\ j=1,2,\cdots,M;\ i=1,2,\cdots,N$ denote the variation of the molecules of $S_i$ due to reaction $R_j$. Hence, after the reaction $R_j$ happens, molecules of $S_i$ will be changed from $X_i(t)$ to $X_i(t)+v_{ji}$.

Unlike the deterministic system of ODE, Gillespie/SSA algorithm defines a **propensity function** for each reaction $R_j$, let's say, $a_j$. And $a_j(x)dt$ describes the probability that $R_j$ happens between the time interval $[t, t+dt)$.

Since it's difficult to accurately measure or calculate $a_j$ in any possible manner, we consider the 2 important conditions in thermodynamics to approximate it: **molecule collision**, and **high enough activation energy**. And $a_j$ is the product of the two probablities.

For the former, let's take the simple reaction $A+2B\rightarrow C$ for an example: if there are $x_1$ molecules of $A$ and $x_2$ molecules of $B$ in the system, then the probability of 1 $A$ and 2 $B$ colliding is $\frac{C_{x_1}^1\times X_{x_2}^2}{C_{x_1+x_2}^3}$, which is in direct proportion to $\frac{x_1x_2(x_2-1)}{2}$ .

For the latter, we all know that the activation energy of our reactants should be higher than a "barrier". According to some references, the probability of climbing over this "barrier" is in direct proportion to $e^{\frac{-\Delta \mu}{k_BT}}$, where $\Delta \mu$ refers to the energy barrier, $k_B$ is Boltzmann constant, and $T$ refers to the temperature.

Therefore, we can write $a_j(x)$ like this:
$$
a_j(x)=c_jh_j(x)
$$
where $c_j$ refers to the reaction constant of reaction $R_j$, and $h_j$ refers to the the number of possible combinations of the reactants in $R_j$.

The work above is merely on how to describe a reaction in a system. Ans now we introduce the **Chemical Master Equation**, CME for short, into this topic. For a complex biochemical system, reactions won't happen in any fixed order or any fixed time interval; as long as there are certain amount of some reactants, a reaction may happen, no matter whether it's generated from the previous reactions or pre-set as the initial condition for simulation. It's a stochastic Markov Process.

Supposing that the status of the system at $t_0$ is $X(t_0)=x_0$, we define a probability density function (PDF) $P(x, t|x_0, t_0)$ as the probability of the system's status being $x_0$ at $t_0$. i.e., 
$$
P(x, t|x_0, t_0)=Prob\{X(t)=x,\ if\ X(t_0)=x_0\}
$$
This could also be understood as the percentage of the number of systems with status $x_0$ at $t_0$among all the possible systems. 

With the status of our system at $t_0$ and its PDF at any time, it's safe to conclude the statistical properties of our system at any time. For instance, **the expectation of our system's status** at any time is
$$
<X>=\sum_{x\in \chi}x\times P(x, t|x_0, t_0)
$$
where $\chi$ is the set of all possible status of our system.

Now, the variation of the PDF by time can be defined as
$$
P(x,t+dt|x_0, t_0)-P(x, t|x_0,t_0)
$$
And the transition from $x-v_j$ to $x$ and $x$ to $x+v_j$ in $dt$ can be described respectively as
$$
P(x,t+dt|x_0, t_0)=\sum\limits_{j=1}^MP(x-v_j, t|x_0,t_0)a_j(x-v_j)dt\\
P(x,t|x_0, t_0)=\sum\limits_{j=1}^MP(x, t|x_0, t_0)a_j(x)dt
$$
Combining (9) and (10) and we can get
$$
P(x,t+dt|x_0, t_0)-P(x, t|x_0,t_0)=\sum\limits_{j=1}^MP(x-v_j, t|x_0,t_0)a_j(x-v_j)dt-\sum\limits_{j=1}^MP(x, t|x_0, t_0)a_j(x)dt
$$
Now divide both sides by $dt$ and let $dt \rightarrow 0$, we can arrive at the Chemical Master Equation:
$$
\frac{\partial}{\partial t}P(x, t|x_0,t_0)=\sum\limits_{j=1}^M[P(x-v_j, t|x_0,t_0)a_j(x-v_j)-P(x, t|x_0, t_0)a_j(x)]
$$
It's nearly impossible to obtain an analytical solution for a quite complex biochemical reaction system. But if there are enough **trajectories**, then the evolution process of the system, i.e., the PDF, can be obtained by simulation.

Now let's go back to Gillespie algorithm. It generates trajectories by calculating the propensity of every possible reaction, and choosing one randomly to let it happen in a random time interval, i.e., changing the status of the whole system accordingly. So the next question is, how we decide the random time interval and choose the reaction that will happen.

Define the time interval of a reaction $\tau$ and the reaction we choose $R_{\mu}$, then $(\tau, \mu)$ has to satisfy a given PDF as follows:
$$
P_{0}(\tau, \mu ; x)=\left\{\begin{array}{l}
a_{\mu}(x) e^{-a_{0}(x) \tau}, \text { if } 0 \leq \tau<\infty \text { and } \mu=1, \ldots, M \\
0, \text { otherwise }
\end{array}\right.
$$
where $P(\tau, \mu;x)$ denotes the probability under conditions that the system's status is $x$ at $t$, that the next reaction will be $R_{\mu}$, and will take place at time $t+\tau$. Next we try to prove the PDF $(11)$. We write $P(\tau, \mu;x)$  in this way:
$$
P(\tau, \mu;x) = P_0(\tau, x)\times a_{\mu}(x)d\tau
$$
where $P_0(\tau, x)$ refers to the probability of no reaction happening during $(t, t+\tau)$ when $X(t)=x$, while $a_{\mu}(x)d\tau$ denotes the probability of reaction $R_{\mu}$ happening in $d\tau$, an infinitesimal time interval.

In order to derive the analytical expression of $P_0(\tau, x)$, we consider it as a function of $\tau$, and regard $x$ as a constant. Since there is no reaction taking plat at time $t$, it's apparent that $P_0(\tau, x)=1$.

Let $P_0(\tau^{'}, x)$ be the probability of no reaction happening during $(t, t+\tau^{'})$, and $[1-\sum_{v=1}^Ma_v(x)d\tau^{'}]$ is the probability of no reaction taking place in $d\tau^{'}$, then the probability that no reaction takes place in interval $(t, t+\tau^{'}+d\tau^{'})$ can be written as
$$
P(\tau^{'}+d\tau^{'}, x)=P_0(\tau^{'}, x)\times[1-\sum_{v=1}^Ma_v(x)d\tau^{'}]
$$
Then we can obtain
$$
\frac{\partial}{\partial \tau^{'}}P_0(\tau^{'}, x)=-\sum_{v=1}^Ma_v(x)P_0(\tau^{'}, x)P_0(0,x)=1
$$
Integrate the equation above and we get
$$
P(\tau, x)=e^{-\sum_{v=1}^Ma_v(x)\tau}
$$
Let $a_0=\sum_{v=1}^Ma_v(x)$, and then we have the PDF stated in $(11)$ that $(\tau, \mu)$ has to satisfy.

Eventually, we describe the process of the algorithm as follows:

1.   Initialize the time $t=t_0$ and the status of the system $x=x_0$
2.   With status of the system $x$ at time $t$, calculate all the propensity functions $a_j(x)$, and let $a_0=\sum\limits_{v=1}^Ma_v$
3.   Generate 2 random real numbers $\tau,\ \mu$ that satisfy the PDF $(11)$. More specifically
     1.   we can first choose 2 numbers $r_1, r_2$ that are distributed evenly over $[0, 1]$
     2.   Let $\tau = -\frac{\ln(r_1)}{a_0}$
     3.   Let $\mu$ be an integer that satisfies $\sum_{v=1}^{\mu-1} a_{v}<r_{2} a_{0} \leq \sum_{v=1}^{\mu} a_{v}$
4.   Let $t = t+\tau$, and update the system's status according to $R_{\mu}$, i.e., $X_i\rightarrow X_i+v_{\mu i}$
5.   Repeat Step 2~5 until $t$ reaches the preset simulation time limit

## Simulation Result

We set simulation time limit to 86400s, which is 24hr, the time where THC is taken in is 21600s (6hr) and the time it is fully degraded is set to 43200s (12hr). The simulation result is shown below in a MATLAB figure[^1].

[^1]: When trying simulation, please note that you'd better add breakpoints to line 142 and line 154 with the condition of "~isempty(find(state < 0, 1))". This is because even though we have debugged and improved the robustness of our program, it's still possible (but not much) that a reaction will happen when there are no reactants in the system, and this will cause negative value in the trajectories.

<img src="https://github.com/HalveLuve/Images/blob/master/uPic/image.png?raw=true" style="zoom:50%;" ></img>

## Parameters

Most of the parameters used in our Gillespie-based simulation were found in the literature given above, and we will update the parameter table here in the near future. The rest are either assumed or gained from experiment.