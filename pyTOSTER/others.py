def tost_decision(hypothesis, alpha, pvalue, pTOST, mu_text="zero"):
    if hypothesis == "EQU":
        if pvalue <= alpha and pTOST <= alpha:
            combined_outcome = f"NHST: reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: reject null equivalence hypothesis"
        elif pvalue < alpha and pTOST > alpha:
            combined_outcome = f"NHST: reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: don't reject null equivalence hypothesis"
        elif pvalue > alpha and pTOST <= alpha:
            combined_outcome = f"NHST: don't reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: reject null equivalence hypothesis"
        else:
            combined_outcome = f"NHST: don't reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: don't reject null equivalence hypothesis"
    else: # MET
        if pvalue <= alpha and pTOST <= alpha:
            combined_outcome = f"NHST: reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: reject null MET hypothesis"
        elif pvalue < alpha and pTOST > alpha:
            combined_outcome = f"NHST: reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: don't reject null MET hypothesis"
        elif pvalue > alpha and pTOST <= alpha:
            combined_outcome = f"NHST: don't reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: reject null MET hypothesis"
        else:
            combined_outcome = f"NHST: don't reject null significance hypothesis that the effect is equal to {mu_text} \nTOST: don't reject null MET hypothesis"

    return combined_outcome
