Reviewer 1 comments
================================================================================

This paper considers timely questions concerning the stability of complex systems, mainly using random matrix theory. The central problem consists in studying the stability of random matrices which were originally seen as random community matrices from ecology. In a seminal paper, R May challenged the consensus that complexity begets stability in complex ecosystems. These issues are of current concern in ecology. Most of the works dealing with such problems (either for stability or feasibility) consider networks of high species richness S, and develop or use techniques and results from random matrix theory in the large S limit, like e.g. Girko circular law that describes the asymptotic behavior of the related spectra measures (a system is stable when all of its eigenvalues have negative real parts).

The paper of Duthie is interesting since it deals with finite size systems where S is not too large. Mathematically, the author focus on random matrices of the form D A-id where D is a random diagonal matrix and A is a random matrix having i.i.d. centered entries. The author suggests numerically that stability is enhanced by multiplying A by D. This idea seems to be new. 

Recent works dealing with similar models in the large S limit have been published recently (Gibbs et al (2018) and Ahmadian (2018). Such works are not dealing with finite S and focus on limiting spectral measures. On the other hand, real empirical ecological networks have finite sizes S, so that the paper might be of practical interest for ecologists dealing with empirical data. In this setting, most interaction weights that define Lotka-Volterra dynamics are unknown at present time, and a natural idea consists in randomizing such weights to arrive at random community matrices of finite size, see, e.g. Ensemble ecosystem modeling for predicting ecosystem response to predator reintroduction, Baker et al. (2016) Conservation Biology, Volume 31, No. 2, 376–384.

In summary I think that the paper is interesting and timely, but I strongly recommend that the author
considers also empirical webs (like for example those given in Dougoud et al. (2018)), the feasibility
of equilibria in large ecosystems: A primary but neglected concept in the complexity-stability debate,
PLoS Comput. Biol, and check numerically if the observation that variations of component responses
rates enhance stability and feasibility is valid in experimentally observed webs. The author should
also test this hypothesis on standard webs models for structured food web models like the niche, or
the nested-hierarchy random web models.

Comment:

The manuscript is well written. I have just a remark concerning Figures 3 and 4. The black line is
intended to show the proportion of systems that are stable when the variance is positive but would
be unstable for zero variance. Is the author considering the conditional probability that A-id is
unstable given that DA-id is stable ? I think that this point should be explained in more details.


Reviewer 2 comments
================================================================================

The manuscript explores the stability of matrices of the form

A = diag(gamma) * M,

where gamma is an S-component vector, diag(gamma) is the diagonal matrix with gamma along its diagonal, * denotes matrix multiplication, and M is a circular random matrix (independently and randomly sampled offdiagonal entries with mean 0 and variance V=C*sigma^2 where C is the connectance and sigma^2 the variance of nonzero entries; and with all diagonal entries equal to -1). The main claim is that whenever the random matrix M has parameters putting it near the boundary of stability and instability, the likelihood of A being stable is increased by having a larger variance in gamma.

The manuscript is very clear and its results were easily reproducible. My concerns are threefold: first, there is no attempt at an analytical understanding, making the breadth of the results difficult to judge; second, the numerical exploration is restricted to a very narrow class of random matrices, further hampering a good appreciation of the results' breadth; and third, it feels as if the effect of response rate variation was a second- or even third-order effect, and thus only of minor importance in explaining the stability of large complex systems. Below I elaborate on each of these concerns in turn.

The results crucially build on the matrices being dangerously close to the stability boundary. At that point, minor fluctuations in the position of the leading eigenvalue can make or break a system. In fact, there are known results on the uncertainty in the leading eigenvalue's position - but, to my knowledge, those are only available for circular matrices (here: constant gamma). I am wondering whether combining those with the methods of Gibbs et al. (2018) could perhaps be used to gain further understanding. In either case, the manuscript would of course be much stronger with analytical results, but I do realize this may be uncharted territory.

Since numerical explorations are used in lieu of analytical ones, I would look at a broader class of matrices M. For example, the average pairwise correlation rho between the (i,j)th and (j,i) th matrix entries is known to strongly affect the local stability of a random matrix - in the large S limit, the leading eigenvalue converges to 

lambda = sqrt(S*V)*(+rho).

I did some preliminary exploration of having different rho values, and this did not change the Author's results. So adding this new aspect to the manuscript would strengthen its results. I would also maybe consider matrices with even more structure - e.g., varying diagonal entries in M, or having different means/variances in the top and bottom triangles (corresponding to food webs with cascade structures; Allesina et al. 2015 Nature Communications).

Finally, the general relevance of the results is hampered by being confined to a narrow parameter space, and even there the fraction of occurrences where var(gamma) makes a difference is minuscule. Looking at Figure 3, the only cases in which there is an appreciable effect of a variance in the gammas is when only a small handful of the matrices end up being stable in the first place. As such, var(gamma) only makes a difference when there is not much difference to be made anyway, because the vast majority of matrices end up unstable regardless of the magnitude of var(gamma). In light of this, I would be much more careful with the manuscript's framing - in particular, I would change the title to express better that we are dealing with second-order effects, and would not say in the Abstract that system stability is "markedly increased" by variation in the gammas. 

Reviewer 3 commments
================================================================================

This manuscript describe the effect of variation in the response rate of individual systems to the stability of a complex system by A. Bradley Duthie. The manuscript is well-written and the calculations, codes and so on sound technically correct, but I have problem to see the importance of the response rates to understand the stability of the systems.

In my opinion 3 main questions are more important in all these discussions:

1. The diagonal of the matrix, thinking in species interactions, represents the intraspecific interactions in the same species, and in the study of the stability is the main point because it gives us the center of the cloud of eigenvalues. In these study, as many all of them, is taken as -1. All infraspecific interactions in every species are negatives and equal? It is difficult to agree with this idea.

2. In these discussions the author use random matrix to implement complex systems. Again it is difficult to understand how a random matrix can represent a complex system. If we think in species again, different ecological interactions have different distributions, has can be observed in the paper of Tang and I Sundstrand that the random matrices can be a perfect null model, but nothing more.

3. I can understand the importance of the study of the responds rates to understand the stability of the matrices. Maybe I am wrong but sound like the author is comparing apples and oranges.


Remarks to the author:

This manuscript describe the effect of variation in the response rate of individual systems to the stability of a complex system by A. Bradley Duthie. The manuscript is well-written and the calculations, codes and so on sound technically correct, but I have problem to see the importance of the response rates to understand the stability of the systems.

In my opinion 3 main questions are more important in all these discussions:

1. The diagonal of the matrix, thinking in species interactions, represents the intraspecific interactions in the same species, and in the study of the stability is the main point because it gives us the center of the cloud of eigenvalues. In these study, as many all of them, is taken as -1. All infraspecific interactions in every species are negatives and equal? It is difficult to agree with this idea.

2. In these discussions the author use random matrix to implement complex systems. Again it is difficult to understand how a random matrix can represent a complex system. If we think in species again, different ecological interactions have different distributions, has can be observed in the paper of Tang and I Sundstrand that the random matrices can be a perfect null model, but nothing more.

3. I can understand the importance of the study of the responds rates to understand the stability of the matrices. Maybe I am wrong but sound like the author is comparing apples and oranges.

Other minor questions.

I've missed some classic articles from random-matrix literature as:

H.J. Sommers, A. Crisanti, H. Sompolinsky, and Y. Stein, Spectrum of Large Random Asymmetric Variables, Phys. Rev. Lett. 60, 1895 (1988).

S.E. Townsend, D.T. Haydon, and L. Matthews, On the generality of stability-complexity relationships in Lotka-Volterra ecosystems, J. Theor. Biol. 267, 243 (2010).

A. Mougi, and M. Kondoh, Diversity of Interaction Types and Ecolog- ical Community Stability, Science 337, 349 (2012).

K.Z. Coyte, J. Schluter, and K.R. Foster, The ecology of the micro- biome: Networks, competition, and stability, Science 350, 663 (2015).

