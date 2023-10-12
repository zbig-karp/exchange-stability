Stability of dyadic exchange: experimental evidence for the impact of
shared group membership
================

This repository contains data and R code necessary to reproduce results
reported in Karpiński et al. (2023). It consists of just three files:

- the files `social-exchange-data.csv` and `social-identity-data.csv`
  contain the data used in the analyses. More information about the
  variables is given below. In both files, columns (variables) are
  separated by semicolons;
- the file `data-analysis.R` is an R script file with R code used in the
  analyses of the data.

The data come form an experiment in which participants were divided into
two categories, based on their preference for the paintings of Kandinsky
over Klee, and then were given 20 opportunities to share valued
resources with one another. After completing all 20 opportunities, the
participants filled in a short questionnaire concerning their
perceptions of their in-group (i.e., the category they had been assigned
to) and their out-group. See the paper by Karpiński et al. (2023) for a
complete description of the experimental setting, protocols, and
instructions.

# The data

## Social exchange process

The data in `social-exchange-data.csv` come from the part of the
experiment in which participants were sharing resources with other
participants over 20 opportunities (or rounds). The data come in a long
format, each observation being an instance of participant $i$ sharing
their resources with participant $j$ in round $t$ of experimental
session $k$. With 18 experimental sessions, 10 participants per session,
20 opportunities for sharing resources given to each participants, and 9
alters to share one’s resources with per round, there is a total of
$18\times 20\times 10\times 9$ = 32,400 unique observations.

The following variables are included in the dataset:

- `sessionid` – session identifier;
- `session_date` – calendar date of the day on which the experiment took
  place;
- `visible` – a binary variable equal to “Yes” for experimental sessions
  in which participants’ social categories were made immediately visible
  to them, and “No” otherwise;
- `split` – a character string corresponding to the way in which
  participants were split between social categories. It assumes the
  following values: “5:5” (meaning that all 10 participants in a session
  were evenly split among the categories), “6:4” (meaning that 6 of the
  10 participants in a session were assigned to one category, while the
  remaining 4 ended up in the other), and “7:3” (meaning that 7 of the
  10 participants in a session were assigned to one category, while the
  remaining 3 ended up in the other);
- `subjectid` – subject’s identifier
- `round` – exchange opportunity, an integer ranging from 1 to 20
- `subject` – subject’s index number in the experimental session, an
  integer ranging from 1 to 10;
- `alter` – alter’s index number in the experimental session, an integer
  ranging from 1 to 10;
- `gave` – a binary variable specific to each subject-alter dyad. It is
  coded 1 if the subject chose to share their resources with the alter
  at a given exchange opportunity. Otherwise, it is coded as 0.
- `received` – a binary variable specific to each subject-alter dyad. It
  is coded 1 if the subject received resources from the alter at a given
  exchange opportunity. Otherwise, it is coded as 0.
- `gave_last` – a lagged version of `gave`: it is coded 1 if the subject
  shared their resources with the alter at the **previous exchange
  opportunity**. Otherwise, it is coded as 0. Note that this variable is
  not defined when `round == 1`.
- `received_last` – a lagged version of `received`: it is coded 1 if the
  subject received resources from the alter at the **previous exchange
  opportunity**. Otherwise, it is coded as 0. Note that this variable is
  not defined when `round == 1`.
- `s_cat` – subject’s social category, coded 1 for the Kandinsky group
  and 2 for the Klee group;
- `a_cat` – alter’s social category, coded 1 for the Kandinsky group and
  2 for the Klee group;
- `dyad_type` – a binary variable distinguishing between “Intra-group”
  (where the subject and alter belong to the same social category) and
  “Inter-group” (where the subject and alter belong to different
  categories) dyads.

## Social identity questionnaire

The dataset in `social-identity-data.csv` contains responses to a
questionnaire that the participants filled in after having completed all
20 exchange opportunities. In the questionnaire, they were asked how
they felt about their in-group and their out-group along 4 dimensions:
belongingness, commonality, closeness, and liking. The responses were
coded on 7 point scale, with higher values indicating more positive
evaluations. The questionnaire is based on the work by Yamagishi and
Kiyonari (2000) and Aksoy (2015).

The data come in a long format, each observation corresponding to an
evaluation along dimension $k$ by individual $i$ in experimental session
$j$. With 18 session, 10 participants per session, and 4 dimensions of
evaluation, the total number of observations is 720.

The dataset includes the following variables:

- `sessionid` – session identifier;
- `session_date` – calendar date of the day on which the experiment took
  place;
- `visible` – a binary variable equal to “Yes” for experimental sessions
  in which participants’ social categories were made immediately visible
  to them, and “No” otherwise;
- `split` – a character string corresponding to the way in which
  participants were split between social categories. It assumes the
  following values: “5:5” (meaning that all 10 participants in a session
  were evenly split among the categories), “6:4” (meaning that 6 of the
  10 participants in a session were assigned to one category, while the
  remaining 4 ended up in the other), and “7:3” (meaning that 7 of the
  10 participants in a session were assigned to one category, while the
  remaining 3 ended up in the other);
- `subjectid` – subject’s identifier
- `subject` – subject’s index number in the experimental session, an
  integer ranging from 1 to 10;
- `s_cat` – subject’s social category, coded 1 for the Kandinsky group
  and 2 for the Klee group;
- `female` – a dummy variable for being female, coded 1 for female and 0
  for male;
- `age` – subject’s age (in years);
- `domain` – evaluation domain for the in-group and the out-group;
- `eval_in` – subjective evaluations of one’s in-group;
- `eval_out` – subjective evaluations of one’s out-group.

The variables `sessionid`, `session_date`, `visible`, `split`,
`subjectid`, `subject`, `s_cat` are common to both datasets.

# Acknowledgement

Karpiński et al. (2023) use data from an experiment designed to study
homophily processes in a small-group social exchange setting. The
experiment is a part of a larger research project on **attraction to
similar others** and **repulsion from dissimilar others** as drivers of
homophily in social relations. The research is supported by a grant from
National Science Centre in Poland (grant number
UMO-2017/25/B/HS6/02543).

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-aksoy2015" class="csl-entry">

Aksoy, Ozan. 2015. “Effects of Heterogeneity and Homophily on
Cooperation.” *Social Psychology Quarterly* 78 (4): 324–44.

</div>

<div id="ref-yamagishi2000" class="csl-entry">

Yamagishi, Toshio, and Toko Kiyonari. 2000. “The Group as the Container
of Generalised Reciprocity.” *Social Psychology Quarterly* 63: 116–32.

</div>

</div>
