double fugacity(P,T,buffer) {
    
    switch buffer
        case: "qfm"
            // Compute gibbs free energy
        case: "iw"
            // Compute gibbs free energy
        otherwise:
    }

    // Philpotts & Ague, Principles of Igneous and Metamorphic
    // Petrology 2nd ed. pg. 262 Eq. 11.41
    fO2 = exp(DeltaG/(R*T));

    return(fO2);
}
