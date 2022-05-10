# Abstract

To be implemented...


# Introduction



# Physics

# Building up the Algorithm
To create the algorithm, we begin with diffusion limited aggregation (DLA). A statistical mechanics growth model initially introduced by Witten and Sander \cite{wittensander}, random walkers are sent out from infinity to a central cluster; upon being adjacent to a particle in the cluster, the random walker sticks and becomes part of the cluster, thus causing the nucleus to grow. The name comes from the concept that the particles' movement is defined by diffusion; it is diffusion-limited because the particles that are growing the aggregate are assumed to be in a low enough concentration that they never come in contact with each other; and the central site where the growth is is considered to be the aggregate. As the number of random walkers grows, the aggregate begins to take on a fractal-like shape, as shown in Figure \ref{fig:DLACluster}. DLA is applicable to a plethora of situations: Hele-Shaw flow, dendritic growth, and dielectric breakdown, to name a few \cite{martineau}.

<!-- 
    Insert figure here.
    Filepath: sections/Images/Extraneous/DLA.png
    Caption: An example of a cluster that is formed through Diffusion Limited Aggregation (DLA), courtesy of \cite{paulbourke}. Random walkers congregate against an aggregate and grow outwards in a fractal-like pattern.
-->

## Connection of DLA to Physics
The governing physics begins with the continuity equation, which tells us that water is conserved throughout the process. Explicitly, this states that:
\begin{equation}
    \dv{c}{t} = -\nabla \cdot J,
\end{equation}
where \(c\) is the concentration of water at a node and \(J\) is the flux, which can be expressed as:
\begin{equation}
    J = D \nabla \cdot \hat{n},
\end{equation}
where \(D\) is some diffusion constant, \(\nabla c\) is the gradient of water on the nucleus's surface, and \(\hat{n} = \nabla{S}/\norm{\nabla{S}}\) is the surface normal vector \cite{PREFrosting}. We can put these two relations together to end up with the diffusion equation: 
\begin{equation}
    \dv{c}{t} = - D \nabla^2 c.
\end{equation}
These variables can be shown in Figure \ref{fig:fluxGraphic}. 

<!-- 
    Insert figure here.
    Filepath: sections/Images/Extraneous/Flux.png
    Caption: The flux felt from an ice particle to the water droplet, as shown in \cite{PREFrosting}. These refer to the Equations derived in \ref{section:dlatophysics}
 -->

Following the method outlined in \cite{richmond}, we can connect DLA to the diffusion equation. Specifically, let there be a two-dimensional array with lattice spacing \(\Delta x\). Furthermore, define \(P_{ij}(n)\) as the probability of finding the walker at the position \(ij\) after \(n\) steps. Then, at the previous step there is an equal probability of finding the walker at any of the neighboring sites, granted we only consider the "cardinal" directions, and not any movement along diagonals. This can be written as:
\begin{multline}
    P_{ij}(n) = \frac{1}{4}\left[P_{i-1,j}(n-1) + P_{i+1,j}(n-1) + \right. \\ \left. + P_{i,j-1}(n-1) + P_{i,j+1}(n-1)\right],
\end{multline}
which can be rewritten as: 
\begin{align}
    &P_{ij}(n) - P_{ij}(n-1) = \\ \nonumber
    &= \frac{1}{4}\left[ P_{i-1,j}(n-1) - 2P_{i,j}(n-1)+ P_{i+1,j}(n-1) \right. \\ \nonumber
    &+ \left. P_{i,j-1}(n-1) - 2 P_{i,j}(n-1) + P_{i,j+1}(n-1) \right].
\end{align}
Then, we let \(t = n \Delta t\) and can express:
\begin{align}
    P(n) - P(n-1) &= P\left(\frac{t}{\Delta t}\right) - P\left(\frac{t - \Delta t}{\Delta t}\right) \\
    &\approx \Delta t \dv{P(t)}{t}.
\end{align}
We also let \(x = i \Delta x\), and concern ourselves, without loss of generality to a one dimensional case. Then:
\begin{align}
    P_{i-1} - 2P_i + P_{i+1} \approx (\Delta x)^2 \dv[2]{P(x)}{x}.
\end{align}
Combining the results up to the dimension of our problem, our equation becomes: 
\begin{equation}
    \dv{P(t)}{t} = \frac{(\Delta x)^2}{2N \Delta t} \nabla^2 P(t),
\end{equation}
where \(N\) is the dimension of the lattice. Thus, random walkers can be used to solve a diffusion equation, and the model that we are working with is governed by a diffusion equation. 


## Creating the Algorithm
The underlying algorithm is very similar to the DLA algorithm, however we introduce a probability of sticking such that if a random walker is adjacent to the aggregate it has a probability \(p\) of sticking instead of immediately ensuring that it will stick. 

## Playing Field
The algorithm begins with creating the grid on which the random walkers exist, which we called the "playing field." This was done by defining a length and a width, which we labeled as {\fontfamily{qcr}\selectfont rx} and {\fontfamily{qcr}\selectfont ry}. Then, a domain that is {\fontfamily{qcr}\selectfont 2 rx + 1} by {\fontfamily{qcr}\selectfont 2 ry + 1} is created, such that the origin is at the center and the four corners are  {\fontfamily{qcr}\selectfont (rx,ry)},  {\fontfamily{qcr}\selectfont (rx,-ry)},
{\fontfamily{qcr}\selectfont (-rx,ry)}, and {\fontfamily{qcr}\selectfont (-rx,-ry)}. To keep track of all of the data that we are storing throughout the process, we use {\fontfamily{qcr}\selectfont MATLAB} cell arrays. We initialize the cell array to store five pieces of data: the ordered pair of each lattice point, the current state of the lattice point, the neighbors (in any of the four major directions), the neighbors that are ice, and the neighbors that are dry. All of these concepts will become clearer as this progresses, however, it feels important to mention this is the origin of all of the data. 

Each lattice point is initialized such that their current state is liquid. However, to begin the frosting process, we manually change the state of the origin to be ice. 

## Sending out Random Walkers
Then, we begin to send out random walkers. Random walkers are spawned on a circle of radius {\fontfamily{qcr}\selectfont r} such that angle angle has equal probability of being chosen. Then, however, we define a piecewise function that maps the chosen point to a grid point, such: 
\begin{equation}\label{eq:gridPoint}
    (x,y) \rightarrow
    \begin{cases}
        (\lceil x \rceil, \lceil y \rceil) & \text{if} \hspace{2mm} 0 \leq arg(x,y) < \pi/2 \\
        (\lfloor x \rfloor, \lceil y \rceil) & \text{if} \hspace{2mm} \pi/2 \leq arg(x,y) < \pi \\
        (\lfloor x \rfloor, \lfloor y \rfloor) & \text{if} \hspace{2mm} \pi \leq arg(x,y) < 3\pi/2 \\
        (\lceil x \rceil, \lfloor y \rfloor) & \text{if} \hspace{2mm} 3\pi/2 \leq arg(x,y) < 2\pi
    \end{cases}
\end{equation}
which ensures that the walkers are not on the interior of the circle to begin their random walk, as is shown in Figure \ref{fig:outsideCircle}. This is done such that all walkers are given enough time to do a random walk such that their behavior is still governed by diffusion. 

We also ensure this by providing a "buffer region" between the aggregate and the radius at which the random walkers are being sent out. Let {\fontfamily{qcr}\selectfont buffer} be equal to an arbitrary value; as the aggregate grows, we calculate the further radial point to the highest whole number - that is:
\begin{equation}
    \text{{\fontfamily{qcr}\selectfont furthest}} = \max\{\left\lceil\sqrt{x_i^2 + y_i^2}\right\rceil\}, 
\end{equation}
for all \(i\) in the number of particles. We can succinctly do this by tracking, over each iteration, the radius of the newly introduced particle against the previous iteration's maximum value. Then we can update the value of {\fontfamily{qcr}\selectfont furthest} according to any changes that occur. 

Then, we send walkers out at r = {\fontfamily{qcr}\selectfont buffer + furthest}. This ensures that on each newly sent out random walker, the radius at which they are being released relative to the aggregate is constant. For the runs that we show in the upcoming sections, {\fontfamily{qcr}\selectfont buffer} = 30. 

On each iteration, we also check whether this new radius is equal to the size of the grid that we are working on; if it is, we stop the algorithm. Thus, we are not growing the domain at all and only have it such that we limit the number of particles that are sent out or we abort once the radius is equal to the size of the grid. 

\begin{figure}
    \centering
    \includegraphics[width = \textwidth/2]{sections/Images/Algorithm/circleStarting.png}
    \caption{Random walkers are spawned on a circle of an arbitrary radius (plotted in the figure is r = 5). All angles are equally likely, as is shown in the spawned walkers, which are the black dots. To ensure that these walkers end up on grid points, we choose, according to Equation \ref{eq:gridPoint}, where to send the walkers. These are shown in red. All walkers are pushed outside of the circle. For more information as to why this should be refer to Section \ref{section:spawningWalkers}.}
    \label{fig:outsideCircle}
\end{figure}

## Probability of Freezing
As explored before, the difference between the algorithm developed and DLA is that there is a probability that the random walker freezes instead of, as in DLA, the particle immediately freezing if it is adjacent to the aggregate. This allows different types of aggregates to form, as is explored in Section \ref{section:results}. We let the \(p_{\mathrm{freeze}} \in [0,1]\). 

## Pseudocode
The pseudocode for the algorithm, then, is quite simple.
\begin{lstlisting}
create_playing_field(); 
seed = freeze(0,0);
while (length(stuck) <= particles)
    random_walker();
    while (~stick)
        walk();
        if(on_aggregate)
            if (rand < p_freeze)
                stick();
                stuck += [walker.x, walker.y];
                update_radius();
            else
                walk();
            end
        end
    end
end
\end{lstlisting}

## Further Details and Fine Lines
Some fine points are outlined in this section for the explicit details of the algorithm. For further clarity, refer to the MATLAB code included at the end of this report. 

The walker will only move in the four cardinal directions - up, down, left, or right. 

If the walker is against the aggregate and does not freeze, it does not immediately walk to its next point. Instead, the algorithm checks to see if the walker's next walk will cause it to walk into ice. If it does, the walker does not walk; it will continue to choose a new location to walk to until it finds one that is available. Then, it will walk. 

We limit the number of iterations that a random walker can have during its lifetime. Each time the random walker walks, this number is increased by one. Once the number of iterations exceeds the limit, the particle is discarded and a new random walker is sent out. This ensures that the algorithm does not run for infinite time if the walker somehow avoids the aggregate for long enough. 

If the walker exits the domain, it is sent to the opposite side. For example, if the current \(x\) position exceeds {\fontfamily{qcr}\selectfont nx}, it will be sent to \(x\) = {\fontfamily{qcr}\selectfont -nx}.

\subsection{Algorithm Results}
## Algorithm Results
In the following section, the results of the algorithm with different values of \(p_{\mathrm{freeze}}\) are investigated. In all of the runs, the number of particles is limited to \(N = 3000\), each walker is limited to \(n = 5000\) iterations, the domain size if {\fontfamily{qcr}\selectfont nx = ny = 200}, and the probability of sticking will be explicitly mentioned.

<!-- 
    Insert figure here.
    Filepath: sections/Images/Results/Normal/DLAwithRightColors.png
    Caption: Algorithm run with \(p_{\mathrm{freeze}} = 1\). This exhibits the same properties as DLA, which is as expected.
-->

<!-- 
    Insert figure here.
    Filepath: sections/Images/Results/Normal/noDry_Ice0.5.png
    Caption: Algorithm run with \(p_{\mathrm{freeze}} = 0.5\). As the probability begins to decrease, the branches begin to get thicker.
-->

<!-- 
    Insert figure here.
    Filepath: sections/Images/Results/Normal/noDry_ice0.1.png
    Caption: Algorithm run with \(p_{\mathrm{freeze}} = 0.1\). Decreasing the probability even further makes the arms even larger.
-->

<!-- 
    Insert figure here.
    Filepath: sections/Images/Results/Normal/noDry_ice0.05.png
    Caption: Algorithm run with \(p_{\mathrm{freeze}} = 0.05\). The lower probability really begins to encourage circular growth.
-->

<!-- 
    Insert figure here.
    Filepath: sections/Images/Results/Normal/noDry_ice0.01.png
    Caption: Algorithm run with \(p_{\mathrm{freeze}} = 0.01\). For this low of a probability, the output is much more circular. Although the scale is not shown, this image was much more compact than the previous images.
 -->


All of the results are shown in Figures \ref{fig:noah_dla}, \ref{fig:noah_05}, \ref{fig:noah_01}, \ref{fig:noah_005}, and \ref{fig:noah_001}. For high probabilities, the patterns begin to resemble DLA. For low probabilities, the aggregate is much more circular. This is because high probabilities of sticking encourage the aggregate to grow outwards; once a branch begins to form, it is much more likely to continue to grow. On the other hand, if the probability of sticking is low, then the particles will remain around the aggregate for longer. For example, if the particle is adjacent to the aggregate at some point, and does not stick, it will continue to walk; since the probability of sticking is low, it will most likely continue to walk for a long period of time and ultimately end up causing the aggregate to remain more circular. 

\subsection{Including Drying Regions}
Recall from \ref{section:physics} that the fourth stage of frosting can result in either ice bridging or a dry region forming. If the distance between two neighboring sites is too large, or if the liquid droplet is too small, the liquid will form a bridge, but the particle itself will never freeze. Thus, it will become a dry region. In the algorithm that we have, the domain has liquid particles every unit apart; therefore, this distance will never be an issue. Furthermore, there is no size restriction that we have on the droplets, thus the length will also never be an issue. Therefore, we have to find a way to add it to the algorithm that we already have in order to include dry regions on the model. Recall that the purpose of this project is to simulate results from previous physical studies. Since dry regions are included inherently in the physics, we need to implement it into our algorithm. 

The simplest way of doing this was by simply adding a probability of drying on the {\fontfamily{qcr}\selectfont if} statement that was already implemented in the algorithm from Section \ref{section:psuedocode_noah}.



# Vectorization / Optimization
<!-- TODO: Next section 1 -->

# Results
<!-- TODO: Next section 2 -->


# Discussion


# Conclusion