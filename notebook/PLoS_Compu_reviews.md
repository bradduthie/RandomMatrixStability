Reviewer 1
================================================================================

Reviewer 1 General comments
--------------------------------------------------------------------------------

This paper takes on a timely question: given a random matrix A and diagonal non-negative matrix X, when will M= XA be stable.  It's important b/c this more accurately reflects the dynamics near equilibrium in a system where species abundances are variable (and then are equal to the diagonal elements of X).  The current paper also claims a different interpretation of the stability of a matrix of this form, in connection with variable response rates across species.

Even though this question probably did lie dormant for many years, it has been tackled more recently, and some of the work cited here seems to tackle very similar questions (e.g. Gibbs et al). But those connections/differences aren't really explained here. In short, I'd like to understand if what the author has identified goes beyond what is known from other work about the effect of a diagonal matrix X on the stability of M.

I found the presentation and framing of these results overall very confusing. Here are some comments: 

Response to Reviewer 1 General comments
--------------------------------------------------------------------------------

I am grateful for the helpful comments of Reviewer 1. It is important to emphasise that my focal questions on the stability of complex systems are not restricted to ecological communities (in which all rows and columns of matrix A and M are species densities), but instead address more general questions regarding how the stability of complex systems is affected by varying response rate of system components (with the exception of the specific question that I address concerning ecological community feasibility, of course). System components might be species densities, but they could additionally include a mix of evolving species traits (e.g., Patel et al. 2018) or human social decisions (e.g., management, harvesting, rewilding, land use, etc.; Hui and Richardson 2018). Components might also not include any species, as in the case of purely social systems (e.g., banks; Haldane & May 2011), or physiological or physical systems (e.g., brain networks; Gray & Robinson 2008; 2009). Reviewer 1 is correct that my paper offers a different interpretation of system stability given the diagonal non-negative matrix X, one in which elements X_{i} (\gamma_{i} in the manuscript) are interpreted as the rates at which components respond to perturbation.

The question that I address is distinct from that of Gibbs et al (2018; see also Reviewer 3 comments and responses). I have now revised the manuscript to explain these differences and connections more clearly, both mathematically and conceptually. Mathematically distinct from Gibbs et al. (2018), I am not focused on the effect of a diagonal non-negative matrix X on a stable matrix A as S \to Inf (i.e., as the number of system components becomes arbitrarily large given a stable system). Conceptually distinct, I am not focused specifically on the effect of species abundances on ecosystem stability in a Lotka-Volterra system as species number increases. Rather, I am interested in how the stability of M is affected by X and A (the latter of which is not necessarily assumed stable) at high levels of system complexity, in which the probability of a random complex system being stable becomes rare but non-negligible (under conditions originally identified by May 1972), and for any type of complex system. 


Reviewer 1 specific comment 1
--------------------------------------------------------------------------------

(a) The most obvious interpretation (to me) for gamma_i is as the set of equilibrium abundances for the set of populations, N_i. But  this intepretation is not explained clearly, and indeed would presumably not lead one to restrict the support of gamma to (0,2), or (given what we know about realistic distributions of abundances) to consider cases with half the species abundant, and half very rare.

Response to Reviewer 1 specific comment 1
--------------------------------------------------------------------------------

The use of random matrix theory to investigate the stability of complex systems has perhaps most often been applied to ecological communities, in which the complex system is restricted entirely to interacting species where interactions are described by a square matrix A of dimensions S \times S. This matrix A can then be multiplied by a vector describing the density of each of S species (V) to find the change in species density (V') over some time t: V' = AV (May 1972; 1973). It is therefore understandable why the most obvious interpretation of a vector \gamma_i of length S might be equilibrium species abundances, particularly for those with a background in theoretical ecology; nevertheless, this is an incorrect interpretation of \gamma_i.

In my original manuscript, I was careful to define the vector \gamma as the rate at which the components of a complex system (potentially, but not necessarily, biological species) respond to system perturbation. For the specific question of ecological community feasibility, I also explicitly defined equilibrium abundances as the vector v. I now clarify the use and interpretation of \gamma and v earlier in the manuscript, and explicitly state that \gamma is not interpreted as a vector of species abundances.


Reviewer 1 specific comment 2
--------------------------------------------------------------------------------

(b) If there is another interpretation of gamma_i in terms of evolution of traits, then this is not adequately explained.  Patel et al is cited in support of this, where stability of an eco-evolutionary system is explored using a block matrix of S species and a set of traits that can evolve. But it isn't clear to me in what limit the stability of such a system will reduce to the stability of a matrix of the form M=XA for community matrix A.  Whereas the effect of having a distribution of equilibrium abundances is very clearly of this form.

In summary I would recommend making the connection with species equilibrium abundances clearer, and in parallel including the derivation of gamma_i by considering (presumably extremely fast) trait evolution.

Response to Reviewer 1 specific comment 2
--------------------------------------------------------------------------------

There is another interpretation of \gamma that was explicitly defined in the original manuscript. This interpretation is the rate at which the components of a complex system (note, not just biological species or the values of their traits) respond to system perturbation. I have now revised my manuscript to explain this interpretation more thoroughly, and especially more explicitly with respect to the mathematics; \gamma modifies the matrix A to account for different rates of response inherent to different system components.

Because \gamma is not a vector of equilibrium species abundances, I do not make this connection as requested by Reviewer 1; instead, I have revised my manuscript to make the interpretation of \gamma clearer, and further emphasised that the systems of interest in my manuscript are not restricted to species densities or biological traits.


Reviewer 2
================================================================================

Reviewer 2 General comments
--------------------------------------------------------------------------------

This work it is proposes that accounting for rate variation of nodes in responding to perturbation increases (it is incorrectly stated "drives") the stability in large complex systems.

The paper is well written and presented, and the results are reproducible. However, as I will try to explain, I do not find it meets the criteria for publication in PLOS COMPUTATIONAL BIOLOGY (i.e., Originality, Innovation, High importance to researchers in the field, Significant biological and/or methodological insight, Substantial evidence for its conclusions)

Response to Reviewer 2 General comments
--------------------------------------------------------------------------------

I have changed the title from, "Component response rate variation drives stability in large complex systems" to "Component response rate variation increases the opportunity for stability in large complex systems". I agree with Reviewer 2 that the previous title was technically incorrect.

I am grateful for the comment on the quality of the writing and reproducibility of my results. And while I believe that this manuscript reports original and important work, this helpful review has helped convinced me that it would ultimately fit better in a more multi-disciplinary journal.


Reviewer 2 specific comment 1
--------------------------------------------------------------------------------

1) The Jacobian matrix A (in the notation of the author -  see pag. 8) already incorporates heterogeneity of the response rate, as s_i=\sum_j A_ij is not a priori constant, i.e. different nodes do have different response rates based on their connections. Adding an external heterogeneity therefore, should be justified, having in mind a non-linear model that when linearized, gives a Jacobian that is M_j=\gamma_i A_ij.

Response to Reviewer 2 specific comment 1
--------------------------------------------------------------------------------

This is an interesting and helpful comment, as it is of course true that variation in mean row values will arise when constructing random matrices of finite size. And this random variation could quite reasonably be interpreted as variation in component response rates. However, under such conditions, it is important to also emphasise that the *expected* difference between mean row values (and hence component response rates) is still zero; i.e., differences only arise due to error in i.i.d. random sampling for matrix element values. The effect of this error on mean row value differences decreases with increasing matrix size S. 

In fact, we can calculate the expected standard deviation of mean row values -- an metric of the natural variation in component response rates -- given a randomly generated matrix, A. This value is just the standard error of the row means. Hence for default values of $\sigma = 0.4$ and $C = 1$ for, e.g., $S = 10$, standard devation of mean row values should be as follows:

SD_{mean_row_vals} = \sigma / sqrt{S} = 0.4 / \sqrt{10} = 0.1265

We can do the same for small values of $S = 2$ (SD_{mean_row_vals} = 0.2828) and large values of $S = 50$ (SD_{mean_row_vals} = 0.0566). The R code below recreates this for an $S$; note that the -1 values on the diagonal will cause slight deviations in practice. 

```{r}
SD_mn_row_vals <- function(S = 10, iters = 1000, sigma = 0.4){
    sims  <- rep(x = 0, times = iters);
    for (i in 1:iters){
        A0_dat   <- rnorm(n = S * S, mean = 0, sd = sigma);
        A0       <- matrix(data = A0_dat, nrow = S); 
        diag(A0) <- -1;
        rowmns   <- apply(X = A0, MAR = 1, FUN = mean);
        sims[i]  <- sd(rowmns);
    }
    return(mean(sims));
}
eg_run <- SD_mn_row_vals(S = 50);
```

Importantly, these values are small compared to the equivalent values derived from imposing variation in component response rates a priori, and the inclusion of $Var(\gamma)$ adds to the variation that already is present in A. This added variation is ca 0.9523 for my illustrative example in which half $\gamma$ values equal 1.95 and the other half equal 0.05. For simulations in which $\gamma$ is randomly sampled from 0 to 2, the value is ca 0.577. 

I have edited the text of the manuscript now to make it clear that such variation in row means of A exist prior to the inclusion of $\gamma$, and that such variation could reasonably be interpreted as variation in component response rate. I now introduce $\gamma$ as an a priori increase in the expected difference between component response rates and more carefully interpret the difference between results before and after $\gamma$ is included.


Reviewer 2 specific comment 2
--------------------------------------------------------------------------------

2) The latter case can be realised for a GLV, with \gamma_i=x^*_i. This case corresponds to the case when \gamma_i and A are correlated, and it has been recently fully analysed [ref].

Response to Reviewer 2 specific comment 2
--------------------------------------------------------------------------------

Reviewer 2 here points out the approach in which a vector affecting the community matrix is interpreted as equilibrium species abundances ($x^*_{i}$) in a generalised Lotka-Volterra ecological model, as, e.g., in Gibbs et al (2018). In my Response to Reviewer 3 specific comments 1 and 2, I explain how my manuscript differs from analyses of species equilibrium abundances.


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

Reviewer 2 states that the numerical results both do not produce substantial evidence for the manuscript's conclusions (sentence 1), but also acknowledges an increase in stability ("at best around 2-3%"; sentence 2 -- the attached figures replicating my analysis were much appreciated). Technically this is a contradiction, so I interpret the reviewer to mean that, while stability does in fact increase, the actual magnitude of this increase is not large enough to be interesting. While I agree that the effect here is modest, I would argue that it is an important effect to recognise when considering whether or not complex systems are predicted to be stable, and one that has very broad implications that might affect any complex system (i.e., not just ecological communities).

I disagree that the unstable matrices are typically more unstable when $\gamma \neq 1$. 


Reviewer 3
================================================================================

Reviewer 3 General comments
--------------------------------------------------------------------------------

Reviewer #3: In this contribution, the author studies the variability of component response rates (different time scales) in a general model of population dynamics. He finds that the (asymptotic) stability of the system is largely increased if variability in rates is introduced, compared to the stability of a system where the vector of rates is constant. Although the paper is well written, I have serious concerns about its merits for publication in PLoS Computational Biology.

Response to Reviewer 3 General comments
--------------------------------------------------------------------------------

I am grateful for Reviewer 3's comments. In my revised manuscript, I further emphasise that my model is not restricted to population dynamics, as Reviewer 3 incorrectly interpreted, but encompasses any complex system with interacting components. Most of Reviewer 3's additional concerns are based on a comparison with Gibbs et al. (2018; Phys Rev E. 98:022410. doi: 10.1103/PhysRevE.98.022410. pre-print: https://arxiv.org/pdf/1708.08837.pdf). I have revised my manuscript to further emphasise the difference between my work and Gibbs et al. (2018), and I address more specific concerns below.


Reviewer 3 specific comment 1
--------------------------------------------------------------------------------

First, although the motivation of the author was (in principle) different from that of Gibbs et al., PRE (2018), at the end there is substantial overlap with the latter paper. Gibbs et al. conducted a thorough study of the effect of random population abundances on the stability of community matrices for generalized Lotka-Volterra systems, which reduces to analyzing the stability of the matrix M=XA, where X is a diagonal matrix with random (positive) abundances and A is a (random) interaction matrix. This turns out to be the same problem discussed in the present manuscript, so novelty is substantially reduced (in part because the analysis conducted by Gibbs et al. is deeper and much more extensive).

Response to Reviewer 3 specific comment 1
--------------------------------------------------------------------------------

Reviewer 3 acknowledges that my motivations were different in principle from Gibbs et al. (2018), and I have revised my manuscript to further emphasise these differences while also acknowledging where overlap occurs. My manuscript does not focus on primarily random population abundances (though considers them in sections on feasibility), but rather the reaction rates of system components to system perturbation. While the mathematics of underlying my approach is based on analyses of the stability of the form M = XA, the diagonal matrices are not interpreted as abundances, and my simulations do not address the same theoretical question as in Gibbs et al. (2018). The novelty of my manuscript is conceptual rather than mathematical, as is also the case for Gibbs et al. (2018). Indeed, Ahmadian et al. (2015; Phys. Rev. E Stat. Nonlin. Soft. Matter Phys. 91:012820. doi:10.1103/PhysRevE. 91.012820) provided a framework for investigating eigenvalue densities of random matrices of the even more general form M = B + XAY, which would include Gibbs et al.'s (2018) M = XA. What was unique about Gibbs et al. (2018) was its application of the properties of random matrices of this form to a long-standing question about the relationship between community stability and community feasibility. They showed that given an interaction matrix (A) is stable, the vector of population abundances (X) becomes decreasingly likely to destabilise the resulting community matrix (M) as the community size increases. 

My work applies the same mathematical framework to a different theoretical question, and I have revised my manuscript to make this distinction clear.


Reviewer 3 specific comment 2
--------------------------------------------------------------------------------

Second, the main point of the present manuscript is that stability is largely increased by introducing variability in response rates. Although the author cites the work by Gibbs et al., is surprising that he didn't acknowledge that his results are contradictory with those provided by Gibbs et al. The key point of the latter paper is that population abundances have no effect on stability, as May assumed in his seminal paper. In fact, Gibbs et al. state that 'the community matrix is stable if an only if the interaction matrix is stable. In other words, the abundance of species seems to not affect the sign of eigenvalues'. This statement has to be interpreted in the limit of large systems (S going to infinity) and in probabilistic terms, i.e., it is true almost surely. In fact, Gibbs et al. show that the probability of matrix M=XA being unstable given that A is stable decays exponentially with S (see Fig. 5 of Gibbs el al.) Indeed, they showed this analytically for uniformly distributed abundances X and constant self-regulation terms, which is precisely the case studied in the present manuscript. Conversely, the probability of M becoming stable given that A is unstable also decays exponentially with S (Fig. 9 of Gibbs et al.). This means that multiplying a random matrix by a diagonal random matrix has no effect on stability in the limit of large S.

Response to Reviewer 3 specific comment 2
--------------------------------------------------------------------------------

Reviewer 3 has slightly mis-stated the main point of my manuscript, in part due to some imprecision in the wording of my previous text that Reviewer 1 noticed (particularly regarding the title). The main point of the manuscript is not that "stability is largely increased by introducing variability in response rates", but that *systems in which there is variability in component response rates are more likely to be stable*. I made this point in the second paragraph of my Discussion:

"It is important to emphasise that variation in component response rate is not stabilising per se; that is, adding variation in component response rate to a particular system does not necessarily increase the probability that the system will be stable. Rather, systems that are observed to be stable are more likely to have varying component response rates, and for this variation to be critical to their stability."

I now also emphasise this point in the abstract and clarify the range of S in which I am interested. Given this distinction, there is no contradiction between my results and the work of Gibbs et al. (2018). Three points are particularly important:

(1) Randomly generated matrices (M = XA) that are stable are more likely than not to have variation in X across the range of system sizes (S) that I simulated. This is (to my knowledge) a novel result and fully reprodible using the code underlying my simulations (https://github.com/bradduthie/RandomMatrixStability). 

(2) Using the data from the first table of my Supporting Information, it is possible to confirm that the probability that a stable matrix (A) is destabilised by (X) such that M = XA is unstable decreases with increasing system size (S), consistent with Gibbs et al. (2018) and Reviewer 3's concerns. I now state this point explicitly in the Discussion and show it in a new figure in the Supporting Information. 

(3) At high S, the probability that a system is stable becomes negligible; no such systems were found given S > 32. I now emphasise that my results apply specifically for the upper range of S > 10 for which stability is also non-negligible (i.e., observed at least once in simulations). In this range, systems in which there is variation in X are indeed more stable than in systems in which no such variation exists. Further, this general pattern occurs for all simulated systems at the upper range of S, regardless of system connectance (C), interaction stengths (\sigma), or X and A element distributions (see Supporting Information). In this manuscript, I am not interested in the effect of varying component response rates in systems that are too large and/or complex to be stable (i.e., "S going to infinity"); my focus is entirely restricted to effects in the range where some stability can be reasonably expected. I explain this more carefully now in my Introduction and revise to state that my general conclusions may not apply as S \to \Inf.

I hope these points adequately address Reviewer 3's concerns. There is no contradiction between my results and those of Gibbs et al. (2018), and I have revised my manuscript to more carefully explain why the two works differ.


Reviewer 3 specific comment 3
--------------------------------------------------------------------------------

This result is seemingly contradictory with the simulations reported in the manuscript under review (Figs. 3 and 4). This is because the author didn't scale the variance of A (sigma^2) with 1/S, indeed he maintained constant variance while increasing S. Increasing S without scaling variances is misleading because eigenvalue distributions are not comparable for different values of S. My feeling is that if he repeats the figures 3 and 4 using this scaling for sigma^2, then the artifactual difference will disappear in the limit of large S, as Gibbs et al found.

Response to Reviewer 3 specific comment 3
--------------------------------------------------------------------------------

Reviewer 3 is of course correct that scaling the variance of A (sigma^2) with 1/S makes the difference between \gamma_{i} = 1 and Var(\gamma_{i}) disappear. This was confirmed by simulation during model development, but was not discussed in the manuscript because it was not a focus of mine (this further highlights the difference between my work and that of Gibbs et al. 2018). 

The reason for increasing values of S is to increase the total complexity of the system, and to isolate and ultimately understand the effect of Var(\gamma_{i}) when S approaches a size at which a random system becomes too complex to have a realistic probability of being stable (less than 1 in 1 million). This is the reason for showing increasing system size (S) on the x-axis of Figures 3 and 4. If I were to have scaled the variance of A (\sigma^{2}) with 1/S as Reviewer 3 suggested, then the complexity of the system would not change, defeating the purpose of increasing S. To my knowledge, there is also no a priori reason to expect complex systems to have variances of interaction strengths (\sigma) that scale to their size S, so forcing this to be the case was not justified.

Finally, while eigenvalue distributions are not identical for different values of S without scaling variances, it is not really accurate to say that eigenvalue distributions are "not comparable"; eigenvalue distributions change with S in a predictable way that can be compared when S is changed. It is well-established that the eigenvalues for a random matrix M such as that used in my manuscript, and originally by May (1972; Nature 238:413-414), will be uniformly distributed within a circle on the complex plane within a radius of \sigma\sqrt{SC} < d, where '-d' is the mean of diagonal elements of M (Tao and Tu 2010; Annals of Probability 38:2023-2065). Hence, increasing S specifically increases the radius of this uniform distribution, and does so predictably such that the probability that real eigenvalues are negative (the criteria for system stability) becomes increasingly low with increasing S. Multiplying the random matrix A (which has the uniform distribution with a circle of radius \sigma\sqrt{SC}) by the vector \gamma_{i} in my manuscript resulted in an eigenvalue distribution that was more likely to be stable at the upper range of S. I now clarify the difference between my result and that of Gibbs et al. (2018) in this context in my revised Discussion.


Reviewer 3 specific comment 4
--------------------------------------------------------------------------------

The author also states in this manuscript that stability can be largely increased by 'selecting' response rates so that linear stability is forced to increase. This is an obvious technical result but (in the light of the results by Gibbs et al.) it is also obvious that these cases will form a set of null measure in the limit of large systems.

Response to Reviewer 3 specific comment 4
--------------------------------------------------------------------------------

I showed that manipulating component response rates can, sometimes greatly, increase the probability that a random system will be stable. As with Reviewer 3 specific comment 3 (see above), I now clarify that my focus is on the upper range of S at which large random systems can reasonably be expected to be stable, not at the limit of S \to \Infty (as in Gibbs et al. 2018



Author References Cited
================================================================================

Gibbs, T., Grilli, J., Rogers, T., & Allesina, S. (2018). Effect of population abundances on the stability of large random ecosystems. Physical Review E, 98(2), 022410.

Gray, R. T., & Robinson, P. A. (2009). Stability of random brain networks with excitatory and inhibitory connections. Neurocomputing, 72(7--9), 1849--1858. https://doi.org/10.1016/j.neucom.2008.06.001

Gray, R. T., & Robinson, P. A. (2008). Stability and synchronization of random brain networks with a distribution of connection strengths. Neurocomputing, 71(7--9), 1373--1387. https://doi.org/10.1016/j.neucom.2007.06.002

Hui, C., & Richardson, D. M. (2018). How to invade an ecological network. Trends in Ecology and Ecolution, xx, 1--11. https://doi.org/10.1016/j.tree.2018.11.003

May, R. M. (1972). Will a large complex system be stable? Nature, 238, 413--414.

May, R. M. (1973). Qualitative stability in model ecosystems. Ecology, 54(3), 638–641. https://doi.org/10.2307/1935352

Patel, S., Cortez, M. H., & Schreiber, S. J. (2018). Partitioning the effects of eco-evolutionary feedbacks on community stability. American Naturalist, 191(3), 1--29. https://doi.org/10.1101/104505

Haldane, A. G., & May, R. M. (2011). Systemic risk in banking ecosystems. Nature, 469(7330), 351--355. https://doi.org/10.1038/nature09659




