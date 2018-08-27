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
ABSumFile = "at.sum"

# A(Complex,Free) / B(Complex,Free)
AFreeSumFile = "t.sum"
BFreeSumFile = "a.sum"

# A(Opt,Free) / B(Opt,Free)
AOptSumFile = "t_opt.sum"
BOptSumFile = "a_opt.sum"

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
        if(refA == "all") {
            AFree[NAFree] = $1;
            AFreex[NAFree] = $3;
            AFreey[NAFree] = $4;
            AFreez[NAFree] = $5;
            NAFree++;
        } else {
            for(atom in setRefA) {
                if(setRefA[atom] == NACountFree) {
                    AFree[NAFree] = $1;
                    AFreex[NAFree] = $3;
                    AFreey[NAFree] = $4;
                    AFreez[NAFree] = $5;
                    NAFree++;
                    break;
                }
            }
        }
        NACountFree++;
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
        if(refB == "all") {
            BFree[NBFree] = $1;
            BFreex[NBFree] = $3;
            BFreey[NBFree] = $4;
            BFreez[NBFree] = $5;
            NBFree++;
        } else {
            for(atom in setRefB) {
                if(setRefB[atom] == NBCountFree) {
                    BFree[NBFree] = $1;
                    BFreex[NBFree] = $3;
                    BFreey[NBFree] = $4;
                    BFreez[NBFree] = $5;
                    NBFree++;
                    break;
                }
            }
        }
        NBCountFree++;
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
        for(atom in AFree){
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
        if(refA == "all") {
            AOpt[NAOpt] = $1;
            AOptx[NAOpt] = $3;
            AOpty[NAOpt] = $4;
            AOptz[NAOpt] = $5;
            NAOpt++;
        } else {
            for(atom in setRefA) {
                if(setRefA[atom] == NACountOpt) {
                    AOpt[NAOpt] = $1;
                    AOptx[NAOpt] = $3;
                    AOpty[NAOpt] = $4;
                    AOptz[NAOpt] = $5;
                    NAOpt++;
                    break;
                }
            }
        }
        NACountOpt++;
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
        if(refB == "all") {
            BOpt[NBOpt] = $1;
            BOptx[NBOpt] = $3;
            BOpty[NBOpt] = $4;
            BOptz[NBOpt] = $5;
            NBOpt++;
        } else {
            for(atom in setRefB) {
                if(setRefB[atom] == NBCountOpt) {
                    BOpt[NBOpt] = $1;
                    BOptx[NBOpt] = $3;
                    BOpty[NBOpt] = $4;
                    BOptz[NBOpt] = $5;
                    NBOpt++;
                    break;
                }
            }
        }
        NBCountOpt++;
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
}

