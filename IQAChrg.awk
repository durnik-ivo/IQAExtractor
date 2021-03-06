#############################################################################
#
# IQAChrg
#
# Script for extracting atomic charges
#
# Tutorial:
#       1) Fill in user edit section
#       2) Run "awk -f IQAChrg.awk" in the terminal
#       3) For output into file run as "awk -f IQAExtractor.awk > file.txt"
#       4) Carefully check results, especially your definition of fragments and sum files
#
# Author: Ivo Durnik
# Contact: durnik@mail.muni.cz
#
#############################################################################
BEGIN {

#############################################################################
# START OF USER EDIT SECTION


# Specify Sum files
#
# AB
ABSumFile = "examples/at.sum"

# A(Complex,Free) / B(Complex,Free)
AFreeSumFile = "examples/t.sum"
BFreeSumFile = "examples/a.sum"

# A(Opt,Free) / B(Opt,Free)
AOptSumFile = "examples/t_opt.sum"
BOptSumFile = "examples/a_opt.sum"

# Specify atom numbers inside ABSumFile
#
# Syntax: "1,2,3" OR "1-3" OR "1-2,3"
fragmentA = "1-15"
fragmentB = "16-30"

# Produce xyz files to check fragment selection? (yes = 1 / no = 0)
#
# If control=1, script will write xyz files extracted from AIMALL sum files 
# (AB.xyz, A.xyz, AFree.xyz, AOpt.xyz, B.xyz, BFree.xyz, BOpt.xyz)
# Complex is splited according to user selection
control = 1;

#############################################################################
# PREPARE ANALYSIS
#############################################################################

# a.u. -> Angstrom
au = 0.529177249;
# a.u. -> kcal/mol
kcal = 627.509608;
# joule -> kcal/mol
J2kcal = 0.000239005736 * 6.02214086 * 10**23
# pi
pi = 3.14159265359;
# permitivity
eps0 = 8.854187817620 * 10**(-12)
# elementary charge
e = 1.602177 * 10**(-19)

# Helper variables
currentFile = "";
i = 0;
fileN = 1;
myNR = 1;

N = 1;

NA = 1;
NB = 1;

NAFree = 1;
NBFree = 1;

NACountFree = 1;
NBCountFree = 1;

NAOpt = 1;
NBOpt = 1;

NACountOpt = 1;
NBCountOpt = 1;

pos = 0;
pos2 = 0;
bool1 = 0;
bool2 = 0;
bool3 = 0;
bool4 = 0;

# Tell awk which files to process (order matters)
ARGV[1] = ABSumFile;
ARGV[2] = AFreeSumFile;
ARGV[3] = BFreeSumFile;
ARGV[4] = AOptSumFile;
ARGV[5] = BOptSumFile;
ARGC = 6;


# Split atom numbers defining molecular fragments into arrays
NsetA = split(fragmentA, setA, ",");
for(i in setA) {
    NSubSet = split(setA[i], subSet, "-");
    if(NSubSet > 1) {
        for(j = subSet[1]; j <= subSet[2]; j++) {
            if(j == subSet[1]){
                delete setA[i];
                NsetA++
            }
            setA[NsetA] = j;
            NsetA++;
        } 
    }
}

NsetB = split(fragmentB, setB, ",");
for(i in setB) {
   NSubSet = split(setB[i], subSet, "-");
    if(NSubSet > 1) {
        for(j = subSet[1]; j <= subSet[2]; j++) {
            if(j == subSet[1]){
                delete setB[i];
                NsetB++
            }
            setB[NsetB] = j;
            NsetB++;
        } 
    }
}

# Split atom numbers for reference states
if(refA != "all") {
    NsetRefA = split(refA, setRefA, ",");
    for(i in setRefA) {
        NSubSet = split(setRefA[i], subSet, "-");
        if(NSubSet > 1) {
            for(j = subSet[1]; j <= subSet[2]; j++) {
                if(j == subSet[1]) {
                    delete setRefA[i];
                    NsetRefA++
                }
                setRefA[NsetRefA] = j;
                NsetRefA++;
            } 
        }
    }
}

if(refB != "all") {
    NsetRefB = split(refB, setRefB, ",");
    for(i in setRefB) {
       NSubSet = split(setRefB[i], subSet, "-");
        if(NSubSet > 1) {
            for(j = subSet[1]; j <= subSet[2]; j++) {
                if(j == subSet[1]) {
                    delete setRefB[i];
                    NsetRefB++
                }
                setRefB[NsetRefB] = j;
                NsetRefB++;
            } 
        }
    }
}

# Check if fragments overlap
for(i in setA) {
    for(j in setB) {
        if( setA[i]==setB[j] ) {
            printf("Fragment selection overlaps!\n");
            printf("A: %s\n", fragmentA);
            printf("B: %s\n", fragmentB);
            exitSwitch = 1;
            exit 1;
        }
    }
}

# Check if control set properly
if( control!=1 && control!=0 ) {
    printf("Control output not set properly!\n");
    printf("1, or 0 expected, but was \"%s\" \n", control);
    exitSwitch = 1;
    exit 1;
}

}

#############################################################################
# COLLECT DATA
#############################################################################
{

#
# Process AB Complex Sum File
#
if(fileN == 1) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:") {
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        AB[N] = $1;
        ABx[N] = $3;
        ABy[N] = $4;
        ABz[N] = $5;
        for(atom in setA) {
            if(setA[atom] == N) {
                A[NA] = $1;
                Ax[NA] = $3;
                Ay[NA] = $4;
                Az[NA] = $5;
                NA++;
                break;
            }
        }
        for(atom in setB) {
            if(setB[atom] == N) {
                B[NB] = $1;
                Bx[NB] = $3;
                By[NB] = $4;
                Bz[NB] = $5;
                NB++;
                break;
            }
        }
        N++;
    }
    if(NF == 0) {
        bool1=0;
    }
    
    # Get atomic charges
    if($0 == "Some Atomic Properties:"){
        bool3=1;
        pos=myNR;
    }
    if( bool3 == 1 && myNR > pos+9){
        posA = 0;
        atomA = -1;
        for( atom in A) {
            if(A[atom] == $1){
                posA = 1;
                atomA = atom;
                break;
            }
        }

        posB = 0;
        atomB = -1;
        for(atom in B){
            if(B[atom] == $1){
                posB = 1;
                atomB = atom;
                break;
            }
        }
        if (posA == 1) {
            AChrg[atomA] = $2;
        }
        if (posB == 1) {
            BChrg[atomB] = $2;
        }
    }
    if( bool3 == 1 && myNR == pos+9+N-1 ) {
        bool3 = 0;
    }

#
# Process A Free Sum File
#
} else if (fileN == 2) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        AFree[NAFree] = $1;
        AFreex[NAFree] = $3;
        AFreey[NAFree] = $4;
        AFreez[NAFree] = $5;
        NAFree++;
    }

    if(NF == 0 ) {
        bool1 = 0;
    }

    # Get Atomic Charges
    if($0 == "Some Atomic Properties:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+9){
        atomA = -1
        for(atom in AFree){
            if(AFree[atom] == $1){
                AFreeChrg[atom] = $2;

            }
        }
    }
    if( bool3 == 1 && myNR == pos+9+NAFree-1 ) {
        bool3 = 0;
    }

   
#
# Process B Free Sum File
#
} else if (fileN == 3) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        BFree[NBFree] = $1;
        BFreex[NBFree] = $3;
        BFreey[NBFree] = $4;
        BFreez[NBFree] = $5;
        NBFree++;
    }
    if(NF == 0) {
        bool1=0;
    }

    # Get Atomic Charges
    if($0 == "Some Atomic Properties:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+9){
        atomB = -1
        for(atom in BFree){
            if( BFree[atom] == $1 ) {
                BFreeChrg[atom] = $2;
            }
        }
    }
    if( bool3 == 1 && myNR == pos+9+NBFree-1 ) {
        bool3 = 0;
    }
#
# Process A(Opt) Sum File
# 
} else if (fileN == 4) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3) {
        AOpt[NAOpt] = $1;
        AOptx[NAOpt] = $3;
        AOpty[NAOpt] = $4;
        AOptz[NAOpt] = $5;
        NAOpt++;
    }
    if(NF == 0 ) {
        bool1 = 0;
    }

    # Get Atomic Charges
    if($0 == "Some Atomic Properties:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+9){
        for(atom in AOpt){
            if(AOpt[atom] == $1){
                atomA = atom;
            }
        }
        AOptChrg[atomA] = $2;
    }
    if( bool3 == 1 && myNR == pos+9+NAOpt-1 ) {
        bool3 = 0;
    }

#
# Process B(Opt) Sum File
# 
} else if (fileN == 5) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }

    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        BOpt[NBOpt] = $1;
        BOptx[NBOpt] = $3;
        BOpty[NBOpt] = $4;
        BOptz[NBOpt] = $5;
        NBOpt++;
    }

    if(NF == 0) {
        bool1=0;
    }

    # Get Atomic Charges
    if($0 == "Some Atomic Properties:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+9){
        for(atom in BOpt){
            if(BOpt[atom] == $1){
                atomB = atom;
            }
        }
        BOptChrg[atomB] = $2;
    }
    if( bool3 == 1 && myNR == pos+9+NBOpt-1 ) {
        bool3 = 0;
    }

#
# In case of error
#
} else {
    print "Error while processing files!"
    exitSwitch = 1;
    exit 1;
}
myNR++;
}

# Change file
ENDFILE {
    fileN++;
    myNR = 1;
}

#############################################################################
# CALCULATE AND PRINT RESULTS
#############################################################################
END {
# Check for error during data colection
if( exitSwitch == 1) {
    exit 1;
}

# Initialize helper variables
# AB
totABChrg = 0.0;

# A(Complex,B)
totAChrg = 0.0;

# B(Complex,A)
totBChrg = 0.0;

# A(Complex,Free)
totAFreeChrg = 0.0;

# B(Complex,Free)
totBFreeChrg = 0.0;

# A(Opt,Free)
totAOptChrg = 0.0;

# B(Opt,Free)
totBOptChrg = 0.0;

# Coulombic interaction
el = 0.0;
elTot = 0.0;
for(i=1; i<NA; i++) {
    for(j=1; j<NB; j++) {
        distX = (Ax[i]-Bx[j])**2;
        distY = (Ay[i]-By[j])**2;
        distZ = (Az[i]-Bz[j])**2;
        dist = sqrt( distX + distY + distZ );
        el = ( 1/(4*pi*eps0) ) * AChrg[i] * e * BChrg[j] * e / ( dist*au*10**(-10) ) * J2kcal; # kcal/mol
        elTot += el;
    }
}

# OUTPUT
#
# Nomenclature
printf(" Nomenclature:\n");
printf(" Fragment(Geometry,Vicinity)\n");
printf(" Geometry - Complex, or Opt\n");
printf(" Vicinity - A, B, or Free\n");
printf("\n");

# Total AB charge check
for(i=1; i<NA; i++) {
    totABChrg += AChrg[i];
}
for(i=1; i<NB; i++) {
    totABChrg += BChrg[i];
}
printf(" Total AB charge: %20.10f\n", totABChrg);
printf("\n");

# Print final output, calculate sums
printf(" Atomic Charges A\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atom        A(Complex,B)        A(Complex,Free)        A(Opt,Free)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NA; i++) {
    printf("%5s %20.10f %20.10f %20.10f\n", A[i], AChrg[i], AFreeChrg[i], AOptChrg[i]);
    totAChrg += AChrg[i];
    totAFreeChrg += AFreeChrg[i];
    totAOptChrg += AOptChrg[i];  
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f\n", totAChrg, totAFreeChrg, totAOptChrg);
printf("\n");
printf("\n");

printf(" Atomic Charges B\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atom        B(Complex,A)        B(Complex,Free)        B(Opt,Free)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NB; i++) {
    printf("%5s %20.10f %20.10f %20.10f\n", B[i], BChrg[i], BFreeChrg[i], BOptChrg[i]);
    totBChrg += BChrg[i];
    totBFreeChrg += BFreeChrg[i];
    totBOptChrg += BOptChrg[i];  
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f\n", totBChrg, totBFreeChrg, totBOptChrg);
printf("\n");
printf("\n");

printf("AB Coulombic Interaction: %20.2f kcal/mol", elTot);
printf("\n");
printf("\n");

# Print xyz files for control of correct selection of FragmentA and FragmentB
if( control == 1 ) {
    printf("\n");
    printf("Control output is switched ON\n");
    printf("Writing:\n")

    # Complex AB
    printf("== AB.xyz\n")
    print N-1 > "AB.xyz";
    print "" >> "AB.xyz";
    for(i=1; i<N; i++) {
        element = AB[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, ABx[i]*au, ABy[i]*au, ABz[i]*au) >> "AB.xyz";
    }

    # Fragment A
    printf("== A.xyz\n")
    print NA-1 > "A.xyz";
    print "" >> "A.xyz";
    for(i=1; i<NA; i++) {
        element = A[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, Ax[i]*au, Ay[i]*au, Az[i]*au) >> "A.xyz";
    }

    # Fragment B
    printf("== B.xyz\n")
    print NB-1 > "B.xyz";
    print "" >> "B.xyz";
    for(i=1; i<NB; i++) {
        element = B[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, Bx[i]*au, By[i]*au, Bz[i]*au) >> "B.xyz";
    }

    # Reference A - Free
    printf("== AFree.xyz\n")
    print NAFree-1 > "AFree.xyz";
    print "" >> "AFree.xyz";
    for(i = 1; i < NAFree; i++) {
        element = AFree[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, AFreex[i]*au, AFreey[i]*au, AFreez[i]*au) >> "AFree.xyz";
    }

    # Reference B - Free
    printf("== BFree.xyz\n")
    print NBFree-1 > "BFree.xyz";
    print "" >> "BFree.xyz";
    for(i = 1; i < NBFree; i++) {
        element = BFree[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, BFreex[i]*au, BFreey[i]*au, BFreez[i]*au) >> "BFree.xyz";
    }

    # Reference A - Opt
    printf("== AOpt.xyz\n")
    print NAOpt-1 > "AOpt.xyz";
    print "" >> "AOpt.xyz";
    for(i = 1; i < NAOpt; i++) {
        element = AOpt[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, AOptx[i]*au, AOpty[i]*au, AOptz[i]*au) >> "AOpt.xyz";
    }

    # Reference B - Opt
    printf("== BOpt.xyz\n")
    print NBOpt-1 > "BOpt.xyz";
    print "" >> "BOpt.xyz";
    for(i = 1; i < NBOpt; i++) {
        element = BOpt[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, BOptx[i]*au, BOpty[i]*au, BOptz[i]*au) >> "BOpt.xyz";
    }

}
}

