#' General Social Survey 2012: Schwartz Values Module
#'
#' The General Social Survey (GSS) collects information from the general public on a wide variety of subjects, including attitudes toward social issues, religion, education, jobs and the economy, government and other institutions, politics, and policy issues. The Schwartz Values module is a part of GSS 2012 questionnaire. This module is a good example of ordinal data, particularly Likert items.
#' 
#' The items are 6-level Likert scale from “Not like me at all” to “Very much like me”. Each question describes a person who adheres to certain values, the respondent is asked to what extent that person is like him.
#' 
#' *Module instruction:* “Now I will briefly describe some people. Please listen to each description and tell me how much each person is or is not like you. Use this card for your answer.”
#' 
#' *Possible responses:*
#' 1.	Not like me at all
#' 2.	Not like me
#' 3.	A little like me
#' 4.	Somewhat like me
#' 5.	Like me
#' 6.	Very much like me
#'
#' @format ## `gss12_values`
#' A tibble with 1,255 rows and 21 columns:
#' \describe{
#'  \item{valorig}{“Doing things in original ways is important”: Thinking up new ideas and being creative is important to (her/him). (S/he) likes to do things in (her/his) own original way.}
#'  \item{valrich}{“Getting rich is important”: It is important to (her/him) to be rich. (S/he) wants to have a lot of money and expensive things.}
#'  \item{valeql}{“Equal opportunity is important”: (S/he) thinks it is important that every person in the world should be treated equally. (S/he) believes everyone should have equal opportunities in life.}
#'  \item{valable}{“Showing abilities is important”: It's important to (her/him) to show (her/his) abilities. (S/he) wants people to admire what (s/he) does.}
#'  \item{valsafe}{“Safety is important”: It is important to (her/him) to live in secure surroundings. (S/he) avoids anything that might endanger (her/his) safety.}
#'  \item{valdiff}{“Doing different things is important”: (S/he) likes surprises and is always looking for new things to do. (S/he) thinks it is important to do lots of different things in life.}
#'  \item{valrule}{“Rules are important”: (S/he) believes that people should do what they're told. (S/he) thinks people should follow rules at all times, even when no-one is watching.}
#'  \item{vallist}{“Listening to different opinions is important”: It is important to (her/him) to listen to people who are different from (her/him). Even when (s/he) disagrees with them, (s/he) still wants to understand them.}
#'  \item{valmod}{“Being modest is important”: It is important to (her/him) to be humble and modest. (S/he) tries not to draw attention to (her/him)self.}
#'  \item{valspl}{“Spoiling oneself is important”: Having a good time is important to (her/him). (S/he) likes to "spoil" (her/him)self.}
#'  \item{valfree}{“Being free and independent is important”: It is important to (her/him) to make (her/his) own decisions about what (s/he) does. (S/he) likes to be free and not depend on others.}
#'  \item{valcare}{“Caring for well-being is important”: It's very important to (her/him) to help the people around (her/him). (S/he) wants to care for their well-being.}
#'  \item{valachv}{“Making achievements is important”: Being very successful is important to (her/him). (S/he) hopes people will recognize (her/his) achievements.}
#'  \item{valdfnd}{“Government's defense of citizens is important”: It is important to (her/him) that the government ensures (her/his) safety against all threats. (S/he) wants the state to be strong so it can defend its citizens.}
#'  \item{valrisk}{“Taking risk is important”: (S/he) looks for adventures and likes to take risks. (S/he) wants to have an exciting life.}
#'  \item{valprpr}{“Doing things properly is important”: It is important to (her/him) always to behave properly. (S/he) wants to avoid doing anything people would say is wrong.}
#'  \item{valrspt}{“Getting respect is important”: It is important to (her/him) to get respect from others. (S/he) wants people to do what (s/he) says.}
#'  \item{valdvot}{“Devotion to close people is important”: It is important to (her/him) to be loyal to (her/his) friends. (S/he) wants to devote (her/him)self to people close to (her/him).}
#'  \item{valeco}{“Ecology or environment is important”: (S/he) strongly believes that people should care for nature. Looking after the environment is important to (her/him).}
#'  \item{valtrdn}{“Tradition is important”: Tradition is important to (her/him). (S/he) tries to follow the customs handed down by (her/his) religion or (her/his) family.}
#'  \item{valfun}{“Having fun is important”: (S/he) seeks every chance (s/he) can to have fun. It is important to (her/him) to do things that give (her/him) pleasure.}
#' }
#' @references Smith, Tom W., Marsden, Peter V., and Hout, Michael. General Social Survey, 2012 Merged Data, Including a Cultural Module [United States]. Inter-university Consortium for Political and Social Research [distributor], 2016-05-26. https://doi.org/10.3886/ICPSR35478.v4
#' @source <https://www.icpsr.umich.edu/web/NADAC/studies/35478>
#' @keywords datasets
#' @keywords category
#' @keywords multivariate
#' @keywords nonparametric
#' @concept survey
#' @concept social science
#' @concept ordinal
#' @concept ordinal data
#' @concept ordinal factor
#' @concept likert
#' @concept likert scale
'gss12_values'