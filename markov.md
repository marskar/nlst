## Markov model

Denote R as LCRAT pre-screening risk, E as an dichotomized image feature (like emphysema yes/no),

and ($\alpha$, $\beta$) as regression coefficients to estimate.

 

The current model *interacts* the image feature and log(pre-screening risk):

 

$$log(P(cancer|R,X)) = \alpha log(R) + \beta E log(R)$$


 

The coefficient of the log is equivalent to the exponent on the argument:

$$log(P(cancer|R,X)) = log(R^\alpha) + log(R^{\beta E})$$

 

Exponentiating both sides gives us the probability of cancer (future risk) given prescreening risk (R) and Image Features (X)

 

$$P(cancer|R,X) = R^\alpha + R^{\beta E}$$

 

$$P(cancer|R,X) = R^{\alpha + \beta E}$$

 

The effect of having emphysema depends on your pre-screening risk R.

What if you don’t interact image features and prescreening risk?  You get the model:

 

$$log(P(cancer|R,X)) = \alpha log(R) + \beta E$$

$$P(cancer|R,X) = R^\alpha \cdot exp(\beta E)$$

 

Here the effect of emphysema is the same, *regardless* of your pre-screening risk R.

Hilary and I liked the interacted model because it seemed easy to write down and explain.  But the uninteracted model is also possible, and has a simplicity of its own in that image features affect risk regardless of pre-screening risk.  Frankly, I don’t recall that any reviewer cared about whether the model was easy to explain.

I think we should do model selection on both possibilities.  Model selection will decide on which image features are in group X (whose effect depends on pre-screening risk) and which are in group W that does not (although I guess some features could be in both).  This implies a hybrid risk model:

 

$$P(cancer|R,X,W) = R^{\alpha + X\beta} \cdot exp(W\gamma)$$.

 

Martin, let’s consider this hybrid model and do some model selection to decide if the LCP score is in X or W, or both.  We may also need to reconsider if the current image features should be in W rather than X.
