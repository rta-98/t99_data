Hey Viet,

    >This is a bit much to sort through, but I thought extra data wouldn't hurt 

    >Two directories have been created with README.txt in each to describe more about them. Now generally
      speaking:
        > feature_flags_tor_hyb/ contains .csv with minimal changes to the most recent df I sent you, 
          e.g., bo 1.0 to describe a single bond was left intact; hybridization was specified for atoms 
          in pairs and torsions.
        > feature_bo_dashes/  contains .csv with all changes we discussed from our last meeting, e.g., 
          dashs and equals signs to represent single and double bonds for pairs and torsions and 
          torsions with double central bonds, 

    > Lastly, I used this email as a checklist. If I am missing something please don't hesitate to let me know: 
        (x)  Including molecules that only have C, H, F, and O, no other types of atoms.
        (x)  Count the number of Csp3 and Csp2, H, F, and O.
        (x)  Count the number of bonds and specify the bond order and the hybridization of 
                each atom forming that bond. 
        (x)  For example, the number of Csp3-F bond and the number of Csp2-F bond should 
                be 2 different descriptors; Csp3-Csp3 single bond, 
        (x)  Csp2-Csp3 single bond, Csp2-Csp2 single bond should be 3 different descriptors; 
                Csp2-Csp2 single bond and Csp2-Csp2 double bond should be 2 different things.
        (*x)  Count the number of torsions and specify the bond order and the hybridization of 
                the 2 atoms in the middle of that torsion. 
        (*x)  For example, F-Csp3-Csp2-F, F-Csp3-Csp3-F, F-Csp2-Csp2-F should be 3 different torsions; 
                F-Csp2-Csp2-F and F-Csp2=Csp2-F (double bond in the middle) should be 2 different torsions. 
        (**x)  And as I mentioned above, superimposable torsions should not be 2 different torsions.

        *This feature is present in feature_bo_dashes/
        **Each dataset contains a version in which superimposable torsions are combined, and a version
            in which they are not (normalized is the term I use in the README.txt files under each directory.

-Trey
