Reviewer 1
================================================================================

Reviewer 1 General comments
--------------------------------------------------------------------------------

This paper takes on a timely question: given a random matrix A and diagonal non-negative matrix X, when will M= XA be stable.  It’s important b/c this more accurately reflects the dynamics near equilibrium in a system where species abundances are variable (and then are equal to the diagonal elements of X).  The current paper also claims a different interpretation of the stability of a matrix of this form, in connection with variable response rates across species.

Even though this question probably did lie dormant for many years, it has been tackled more recently, and some of the work cited here seems to tackle very similar questions (e.g. Gibbs et al). But those connections/differences aren’t really explained here. In short, I’d like to understand if what the author has identified goes beyond what is known from other work about the effect of a diagonal matrix X on the stability of M.

I found the presentation and framing of these results overall very confusing. Here are some comments: 

Response to Reviewer 1 General comments
--------------------------------------------------------------------------------





Reviewer 1 specific comment 1
--------------------------------------------------------------------------------

(a) The most obvious interpretation (to me) for gamma_i is as the set of equilibrium abundances for the set of populations, N_i. But  this intepretation is not explained clearly, and indeed would presumably not lead one to restrict the support of gamma to (0,2), or (given what we know about realistic distributions of abundances) to consider cases with half the species abundant, and half very rare.

Response to Reviewer 1 specific comment 1
--------------------------------------------------------------------------------





Reviewer 1 specific comment 2
--------------------------------------------------------------------------------

(b) If there is another interpretation of gamma_i in terms of evolution of traits, then this is not adequately explained.  Patel et al is cited in support of this, where stability of an eco-evolutionary system is explored using a block matrix of S species and a set of traits that can evolve. But it isn’t clear to me in what limit the stability of such a system will reduce to the stability of a matrix of the form M=XA for community matrix A.  Whereas the effect of having a distribution of equilibrium abundances is very clearly of this form.

In summary I would recommend making the connection with species equilibrium abundances clearer, and  in parallel including the derivation of gamma_i by considering (presumably extremely fast) trait evolution.

Response to Reviewer 1 specific comment 2
--------------------------------------------------------------------------------






Reviewer 2
================================================================================

Reviewer 2 General comments
--------------------------------------------------------------------------------

This work it is proposes that accounting for rate variation of nodes in responding to perturbation increases (it is incorrectly stated "drives") the stability in large complex systems.

The paper is well written and presented, and the results are reproducible. However, as I will try to explain, I do not find it meets the criteria for publication in PLOS COMPUTATIONAL BIOLOGY (i.e., Originality, Innovation, High importance to researchers in the field, Significant biological and/or methodological insight, Substantial evidence for its conclusions)

Response to Reviewer 2 General comments
--------------------------------------------------------------------------------





Reviewer 2 specific comment 1
--------------------------------------------------------------------------------

1) The Jacobian matrix A (in the notation of the author -  see pag. 8) already incorporates heterogeneity of the response rate, as s_i=\sum_j A_ij is not a priori constant, i.e. different nodes do have different response rates based on their connections. Adding an external heterogeneity therefore, should be justified, having in mind a non-linear model that when linearized, gives a Jacobian that is M_j=\gamma_i A_ij.

Response to Reviewer 2 specific comment 1
--------------------------------------------------------------------------------








Reviewer 2 specific comment 2
--------------------------------------------------------------------------------

2) The latter case can be realised for a GLV, with \gamma_i=x^*_i. This case corresponds to the case when \gamma_i and A are correlated, and it has been recently fully analysed [ref].

Response to Reviewer 2 specific comment 2
--------------------------------------------------------------------------------







Reviewer 2 specific comment 3
--------------------------------------------------------------------------------

3) There is no analytical analysis of the results, and yet it is known that stability of the linearized system dynamics only depends on the mean, the variance and the correlation of the matrix M (C, or S or sigma are not the right control parameters) [see Grilli et al, Nature Comm 2017]. As the variance of M increases with the heterogeneity in \gamma, it means that when \gamma is different from zero, the induced correlations <M_ij \gammi_i M_ki \gamma_i> have sometimes a stabilize effect. Nevertheless, as it is now, this result is only shown with numerical simulations and there is not an analytical understanding of the result, preventing precise quantitative claims on the effect of the heterogeneity of response rates to the system stability. No mention to diagonal stability theory.

Response to Reviewer 2 specific comment 3
--------------------------------------------------------------------------------







Reviewer 2 specific comment 4
--------------------------------------------------------------------------------

4) The numerical results do no lead any substantial evidence for the conclusion of the manuscript ("rate variation in nodes drive the stability of large complex systems"). In fact, by replicating the analysis (see attached Figures) and by the tables provided by the authors, it is evident that the increase of stability is at best around 2-3% on the overall analysed matrix. Nevertheless, the unstable matrices are typically more unstable for the \gamma \neq 1 case. Too me these results shows that the variation of the stability between gamma_i=1 and heterogeneous gamma_i is weak.

Response to Reviewer 2 specific comment 4
--------------------------------------------------------------------------------

Reviewer 2 states that the numerical results both do not produce substantial evidence for the manuscript's conclusions (sentence 1), but also acknowledges an increase in stability ("at best around 2-3%"; sentence 2). Technically this is a contradiction, so I interpret the reviewer to mean that, while stability does in fact increase, the actual magnitude of this increase is not large enough to be interesting. 




Reviewer 3
================================================================================

Reviewer 3 General comments
--------------------------------------------------------------------------------

Reviewer #3: In this contribution, the author studies the variability of component response rates (different time scales) in a general model of population dynamics. He finds that the (asymptotic) stability of the system is largely increased if variability in rates is introduced, compared to the stability of a system where the vector of rates is constant. Although the paper is well written, I have serious concerns about its merits for publication in PLoS Computational Biology.

Response to Reviewer 3 General comments
--------------------------------------------------------------------------------







Reviewer 3 specific comment 1
--------------------------------------------------------------------------------

First, although the motivation of the author was (in principle) different from that of Gibbs et al., PRE (2018), at the end there is substantial overlap with the latter paper. Gibbs et al. conducted a thorough study of the effect of random population abundances on the stability of community matrices for generalized Lotka-Volterra systems, which reduces to analyzing the stability of the matrix M=XA, where X is a diagonal matrix with random (positive) abundances and A is a (random) interaction matrix. This turns out to be the same problem discussed in the present manuscript, so novelty is substantially reduced (in part because the analysis conducted by Gibbs et al. is deeper and much more extensive).

Response to Reviewer 3 specific comment 1
--------------------------------------------------------------------------------







Reviewer 3 specific comment 2
--------------------------------------------------------------------------------

Second, the main point of the present manuscript is that stability is largely increased by introducing variability in response rates. Although the author cites the work by Gibbs et al., is surprising that he didn't acknowledge that his results are contradictory with those provided by Gibbs et al. The key point of the latter paper is that population abundances have no effect on stability, as May assumed in his seminal paper. In fact, Gibbs et al. state that 'the community matrix is stable if an only if the interaction matrix is stable. In other words, the abundance of species seems to not affect the sign of eigenvalues'. This statement has to be interpreted in the limit of large systems (S going to infinity) and in probabilistic terms, i.e., it is true almost surely. In fact, Gibbs et al. show that the probability of matrix M=XA being unstable given that A is stable decays exponentially with S (see Fig. 5 of Gibbs el al.) Indeed, they showed this analytically for uniformly distributed abundances X and constant self-regulation terms, which is precisely the case studied in the present manuscript. Conversely, the probability of M becoming stable given that A is unstable also decays exponentially with S (Fig. 9 of Gibbs et al.). This means that multiplying a random matrix by a diagonal random matrix has no effect on stability in the limit of large S.

Response to Reviewer 3 specific comment 2
--------------------------------------------------------------------------------








Reviewer 3 specific comment 3
--------------------------------------------------------------------------------

This result is seemingly contradictory with the simulations reported in the manuscript under review (Figs. 3 and 4). This is because the author didn't scale the variance of A (sigma^2) with 1/S, indeed he maintained constant variance while increasing S. Increasing S without scaling variances is misleading because eigenvalue distributions are not comparable for different values of S. My feeling is that if he repeats the figures 3 and 4 using this scaling for sigma^2, then the artifactual difference will disappear in the limit of large S, as Gibbs et al found.

Response to Reviewer 3 specific comment 3
--------------------------------------------------------------------------------



Reviewer 3 specific comment 4
--------------------------------------------------------------------------------

The author also states in this manuscript that stability can be largely increased by 'selecting' response rates so that linear stability is forced to increase. This is an obvious technical result but (in the light of the results by Gibbs et al.) it is also obvious that these cases will form a set of null measure in the limit of large systems.

Response to Reviewer 3 specific comment 4
--------------------------------------------------------------------------------















