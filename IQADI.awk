#############################################################################
#
# IQADI
#
# Script for extracting delocaliztion indexes from AIMALL calculation
# Tutorial:
#       1) Fill in user edit section
#       2) Run "awk -f IQADI.awk" in the terminal
#       3) For text output run as "awk -f IQAExtractor.awk > file.txt"
#       4) Carefully check results, especially your selection of fragments
#
# Author: Ivo Durnik
# Contact: durnik@mail.muni.cz
#
#############################################################################
BEGIN {

#############################################################################
# START OF USER EDIT SECTION

# Specify Sum file
sumFile = "at.sum"

# Specify selection for desired fragments
# Syntax: "1,2,3" OR "1-3" OR "1-2,3"
fragmentA = "1-15"
fragmentB = "16-30"

# Produce xyz files to check fragment selection? (yes = 1 / no = 0)
control = 1;


# END OF USER EDIT SECTION
#############################################################################
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
myNR = 1;
N = 1;
NA = 1;
NB = 1;

pos1 = 0;
pos2 = 0;

bool1 = 0;
bool2 = 0;
bool3 = 0;
bool4 = 0;

# Tell awk which files to process (order matters)
ARGV[1] = sumFile;
ARGC=2;

# Split atom numbers defining molecular fragments into arrays
NsetA = split(fragmentA,setA,",");
for(i in setA) {
    NSubSet = split(setA[i],subSet,"-");
    if(NSubSet > 1) {
        for(j=subSet[1]; j<=subSet[2]; j++) {
            if( j == subSet[1] ){
                delete setA[i];
                NsetA++
            }
            setA[NsetA] = j;
            NsetA++;
        } 
    }
}

NsetB = split(fragmentB,setB,",");
for(i in setB) {
   NSubSet = split(setB[i],subSet,"-");
    if(NSubSet > 1) {
        for(j=subSet[1]; j<=subSet[2]; j++) {
            if( j == subSet[1] ){
                delete setB[i];
                NsetB++
            }
            setB[NsetB] = j;
            NsetB++;
        } 
    }
}

# Check if fragments overlap
for(i in setA) {
    for(j in setB) {
        if( setA[i]==setB[j] ) {
            printf("Fragment selection overlaps!\n");
            printf("A: %s\n", fragmentA);
            printf("B: %s\n", fragmentB)
            exitSwitch = 1;
            exit 1;
        }
    }
}

}

#############################################################################
# COLLECT DATA
#############################################################################
{
#
# Check which file is being processed
#
if(currentFile != FILENAME) {
    fileN++;
    myNR = 1
    currentFile = FILENAME;
}

#
# Process Sum File
#
if(fileN == 1) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos1 = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos1+3){
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

    # Get Delocalization Indexes
    if($0 == "Diatomic Electron Pair Contributions and Delocalization Data:"){
        bool2 = 1;
        pos1 = myNR;

        # Initialize variables for first cycle
        ABDI = 0.0;
    }
    if( bool2 == 1 && NF != 0 && myNR > pos1+14){
        posA = 0;
        atomA = -1;
        for( atom in A) {
            if(A[atom] == $1){
                posA = 1;
                atomA = atom;
            } else if(A[atom] == $2) {
                posA = 2;
                atomA = atom;
            }
        }

        posB = 0;
        atomB = -1;
        for(atom in B){
            if(B[atom] == $1){
                posB = 1;
                atomB = atom;
            } else if(B[atom] == $2){
                posB = 2;
                atomB = atom;
            }
        }
        if( (posA==1 && posB==2) || (posB==1 && posA==2) ) {
            ABDI += $4;
        }
    }
    if( NF == 0 ) {
        bool2 = 0;
    }
} else {
    print "Too many files on the input!"
    exitSwitch = 1;
    exit 1;
}
myNR++;
}

#############################################################################
# CALCULATE AND PRINT RESULTS
#############################################################################
END {
# Check for error during data colection
if( exitSwitch == 1) {
    exit 1;
}

# A(Complex)B(Complex)
totABDI = 0.0;

# Print final output, calculate sums
# AB Interaction
printf("\n");
printf("DI between fragment A and fragment B:\n")
printf("%20.10f\n", ABDI);


# Print xyz files for control of correct selection of FragmentA and FragmentB
if( control == 1 ) {
    printf("\n");
    printf("Control output is ON\n");
    printf("Writing:\n")
    # Complex AB
    printf("== AB.xyz\n")
    print N-1 > "AB.xyz";
    print "" >> "AB.xyz";
    for(i=1; i<N; i++){
        element = AB[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, ABx[i]*au, ABy[i]*au, ABz[i]*au) >> "AB.xyz";
    }

    # Fragment A
    printf("== A.xyz\n")
    print NA-1 > "A.xyz";
    print "" >> "A.xyz";
    for(i=1; i<NA; i++){
        element = A[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, Ax[i]*au, Ay[i]*au, Az[i]*au) >> "A.xyz";
    }

    # Fragment B
    printf("== B.xyz\n")
    print NB-1 > "B.xyz";
    print "" >> "B.xyz";
    for(i=1; i<NB; i++){
        element = B[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, Bx[i]*au, By[i]*au, Bz[i]*au) >> "B.xyz";
    }
}
}
